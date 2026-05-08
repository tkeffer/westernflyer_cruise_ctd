"""EOS-80 physical-oceanography processing for SBE 19plus profiles.

All functions take a per-cast DataFrame and a config dict, mutate the
DataFrame in place where convenient, and return it. The pipeline order
is:

    apply_soak_elimination     (flags surface soak rows)
    apply_ctm_correction       (writes cond_final from cond_raw)
    apply_tau_shift            (writes o2_aligned from o2_umol_l)
    apply_loop_edit            (drops pressure-reversal and slow-descent scans)
    apply_physics              (recomputes SP from cond_final, then theta/rho/o2/etc.)
    apply_chl_calibration      (chlorophyll calibration + low-pass filter)
    apply_wild_edit            (two-pass block-statistics spike removal)
    apply_qc_flags             (1 = good, 3 = high-velocity / bad)

`apply_physics` requires `apply_ctm_correction` and `apply_tau_shift` to
have run first; the previous version silently discarded both corrections.

Sample-rate-dependent calculations read `df.attrs['sample_interval_s']`,
which the ingestion module sets from the .cnv `# interval = ...` line
(falling back to the median of `time_elapsed.diff()`).

References:
    UNESCO 1983, "Algorithms for computation of fundamental properties
    of seawater," Tech. Pap. Mar. Sci. 44.
    Fofonoff & Millard 1983 (theta polynomial).
    Garcia & Gordon 1992, "Oxygen solubility in seawater: Better fitting
    equations," L&O 37(6):1307-1312.
"""

import logging

import numpy as np
import pandas as pd

try:
    import gsw  # only used for SP_from_C; everything else is EOS-80 native
    _HAS_GSW = True
except ImportError:
    _HAS_GSW = False

from seabirdscientific.processing import (
    cell_thermal_mass,
    align_ctd,
    FLAG_VALUE,
    loop_edit_pressure,
    low_pass_filter,
    wild_edit,
    MinVelocityType,
)
from seabirdscientific.eos80_processing import (
    density as _sbs_density,
    potential_temperature as _sbs_potential_temperature,
)
from seabirdscientific.conversion import depth_from_pressure


# Default SBE 19plus scan rate, used only when ingestion couldn't find
# the `# interval =` line. Most fixes for the old hard-coded `* 4` look
# this up via `_sample_rate_hz(df)` instead of using this constant.
DEFAULT_SAMPLE_INTERVAL_S = 0.25


def validate_config(config):
    """Validates calibration.csv parameters before the pipeline runs.

    Logs an error message for each invalid value and raises ValueError
    if any are found, so the build fails fast with a clear message rather
    than a cryptic scipy or arithmetic error deep in processing.
    """
    errors = []

    def _check_positive(key, default):
        try:
            v = float(config.get(key, default))
        except (ValueError, TypeError):
            errors.append(f"{key} is not a valid number")
            return
        if v <= 0:
            errors.append(f"{key}={v} must be > 0")

    def _check_non_negative(key, default):
        try:
            v = float(config.get(key, default))
        except (ValueError, TypeError):
            errors.append(f"{key} is not a valid number")
            return
        if v < 0:
            errors.append(f"{key}={v} must be >= 0")

    def _check_min_int(key, default, minimum):
        try:
            v = int(float(config.get(key, default)))
        except (ValueError, TypeError):
            errors.append(f"{key} is not a valid integer")
            return
        if v < minimum:
            errors.append(f"{key}={v} must be >= {minimum}")

    _check_positive('LPF_TIME_CONSTANT',  0.5)
    _check_positive('CTM_ALPHA',         0.04)
    _check_positive('CTM_TAU',           8.0)
    _check_positive('O2_BOOST_RATIO',    1.0)
    _check_positive('QC_VELOCITY',       1.2)
    _check_positive('BIN_SIZE_METERS',   1.0)
    _check_positive('LOOP_MIN_VELOCITY', 0.25)
    _check_positive('MIN_DEPTH_FLOOR',   10.0)
    _check_non_negative('ALIGN_OXY_SHIFT', 20.0)
    _check_min_int('SOAK_WINDOW_SIZE',   20, minimum=1)
    _check_positive('WILD_STD_PASS1',    2.0)
    _check_positive('WILD_STD_PASS2',    20.0)
    _check_min_int('WILD_SCANS_PER_BLOCK', 100, minimum=1)
    _check_non_negative('WILD_DISTANCE_TO_MEAN', 0.0)

    if errors:
        for msg in errors:
            logging.error(f"calibration.csv validation failed: {msg}")
        raise ValueError(
            f"Invalid calibration.csv — {len(errors)} error(s): {errors}"
        )

    logging.info("calibration.csv validation passed.")


def _sample_rate_hz(df):
    """Returns the per-second sample rate for `df` (Hz).

    Uses the ingestion-stamped `df.attrs['sample_interval_s']` if
    present; otherwise falls back to the elapsed-time median, then to
    the SBE 19plus default of 4 Hz.
    """
    interval_s = df.attrs.get('sample_interval_s')
    if interval_s is None or interval_s <= 0:
        if 'time_elapsed' in df.columns and len(df) > 1:
            diffs = df['time_elapsed'].diff().dropna()
            if not diffs.empty:
                interval_s = float(diffs.median())
    if interval_s is None or interval_s <= 0:
        interval_s = DEFAULT_SAMPLE_INTERVAL_S
    return 1.0 / interval_s


def calculate_depth_eos80(p, lat):
    """Depth (m) from pressure (dbar) and latitude (deg).

    Delegates to the Sea-Bird Scientific toolkit's depth_from_pressure,
    which implements the same UNESCO 1983 algorithm.
    """
    return depth_from_pressure(np.asarray(p, dtype=float), float(np.asarray(lat).flat[0]))


def calculate_theta_eos80(s, t, p, pr=0.0):
    """Potential temperature (Fofonoff & Millard 1983) in IPTS-68 deg C.

    Delegates to the Sea-Bird Scientific toolkit's potential_temperature,
    which implements the identical Runge-Kutta polynomial. The toolkit
    requires pr as an array, so a scalar reference pressure (default 0.0)
    is broadcast to match the shape of p.
    """
    p = np.asarray(p, dtype=float)
    pr_arr = np.full_like(p, float(pr))
    return _sbs_potential_temperature(
        np.asarray(s, dtype=float),
        np.asarray(t, dtype=float),
        p,
        pr_arr,
    )


def calculate_density_eos80(s, t, p):
    """In-situ density (kg/m^3) from PSS-78 salinity, IPTS-68 temp,
    and pressure (dbar).

    Delegates to the Sea-Bird Scientific toolkit's density function,
    which implements the identical UNESCO 1983 polynomial.
    """
    return _sbs_density(
        np.asarray(s, dtype=float),
        np.asarray(t, dtype=float),
        np.asarray(p, dtype=float),
    )


def apply_soak_elimination(df, config):
    """Flags surface-soak rows where the package is still warming up.

    A row is "in soak" until we observe `SOAK_WINDOW_SIZE` consecutive
    samples below `MIN_SAFE_PRESSURE` dbar AND descending faster than
    `LOOP_MIN_VELOCITY` m/s. If no such window exists in the cast,
    every row stays flagged and a warning is logged (the dashboard's
    default `is_soak=0` filter would then hide the cast - that's
    intentional, the cast is suspect).
    """
    soak_depth = float(config.get('MIN_DEPTH_FLOOR', 5.0))
    min_velocity = float(config.get('LOOP_MIN_VELOCITY', 0.25))
    window = int(float(config.get('SOAK_WINDOW_SIZE', 20)))

    rate_hz = _sample_rate_hz(df)
    # dp/dt in dbar/s == m/s for seawater near surface
    dp_dt = df['pres_raw'].diff() * rate_hz
    descent_mask = (df['pres_raw'] > soak_depth) & (dp_dt > min_velocity)

    # First index where the trailing `window` rows are all descending
    rolling_all = descent_mask.rolling(window=window, min_periods=window).sum().eq(window)
    df['is_soak'] = True
    if rolling_all.any():
        # `rolling_all` is True at the *end* of the window, so the soak
        # ends at (idx - window + 1)
        first_end = rolling_all.idxmax()
        first_pos = df.index.get_loc(first_end) - window + 1
        first_pos = max(first_pos, 0)
        df.iloc[first_pos:, df.columns.get_loc('is_soak')] = False
    else:
        logging.warning(
            f"Soak window not found ({window} consecutive descent samples "
            f"@ >{min_velocity} m/s below {soak_depth} dbar); cast will be "
            "flagged entirely as soak."
        )
    return df


def apply_ctm_correction(df, config):
    """Cell-thermal-mass conductivity correction (writes `cond_final`).

    Matches the SBE Data Processing workflow:
        1. Low-Pass Filter on cond_raw (Sea-Bird Scientific toolkit)
        2. Cell Thermal Mass correction (Sea-Bird Scientific toolkit)

    The Sea-Bird Low-Pass Filter is a zero-phase first-order Butterworth:
        a = 1 / (1 + 2*tc/dt),  b = a * (1 - 2*tc/dt)
        filtfilt([a, a], [1, b], data)

    The toolkit's `cell_thermal_mass` implements the recursive SBE App-Note 22
    filter (Morison et al. 1994):

        a = 2α / (dt/τ + 2)
        b = 1 − 2a/α
        ctm[n] = C_smooth[n] + lfilter([a, 0], [1, b], dc_dt × dT)

    where dc_dt = 0.1 × (1 + 0.006 × (T − 20)) and dT = ΔT per scan.

    Downstream `apply_physics` recomputes Practical Salinity from
    `cond_final` so the correction propagates to SP, density, and theta.
    """
    alpha = float(config.get('CTM_ALPHA', 0.04))
    tau   = float(config.get('CTM_TAU', 8.0))
    dt    = 1.0 / _sample_rate_hz(df)
    tc    = float(config.get('LPF_TIME_CONSTANT', 0.5))

    # Step 1: Sea-Bird Low-Pass Filter on conductivity (matches SBE Data Processing).
    cond_smooth = low_pass_filter(df['cond_raw'].to_numpy(), tc, dt)

    # Step 2: Recursive cell-thermal-mass correction (Sea-Bird Scientific).
    cond_ctm = cell_thermal_mass(
        df['temp_raw'].to_numpy(),
        cond_smooth,
        amplitude=alpha,
        time_constant=tau,
        sample_interval=dt,
    )

    df['cond_final'] = pd.Series(cond_ctm, index=df.index)
    return df


def apply_loop_edit(df, config):
    """Remove scans flagged as pressure reversals or slow descent (Loop Edit).

    Delegates to the Sea-Bird Scientific toolkit's `loop_edit_pressure`,
    which replicates the SBE Data Processing Loop Edit algorithm:

    - Computes depth from pressure and latitude (via GSW).
    - Flags scans where vertical velocity < LOOP_MIN_VELOCITY (m/s).
    - Flags scans that fail the monotonic depth test (pressure reversals).

    The toolkit returns a ``cast`` array:  -1 = downcast, 1 = upcast,
    0 = bad.  We discard scans flagged as bad (cast == 0).  main.py
    trims the result to the downcast only (up to max pressure), so both
    downcast and upcast rows must survive this step.

    ``remove_surface_soak=False`` because soak removal is handled
    upstream by ``apply_soak_elimination``.
    """
    min_velocity = float(config.get('LOOP_MIN_VELOCITY', 0.25))
    dt = 1.0 / _sample_rate_hz(df)

    pres  = df['pres_raw'].to_numpy()
    n     = len(pres)
    lat   = float(df['lat'].iloc[0]) if ('lat' in df.columns and not pd.isna(df['lat'].iloc[0])) else 23.0
    flag  = np.zeros(n, dtype=float)

    cast = loop_edit_pressure(
        pressure=pres,
        latitude=lat,
        flag=flag,
        sample_interval=dt,
        min_velocity_type=MinVelocityType.FIXED,
        min_velocity=min_velocity,
        window_size=0.0,           # unused for FIXED type
        mean_speed_percent=0.0,    # unused for FIXED type
        remove_surface_soak=False, # handled by apply_soak_elimination
        min_soak_depth=0.0,
        max_soak_depth=10.0,
        use_deck_pressure_offset=False,
        exclude_flags=False,
    )

    good_mask = cast != 0
    n_bad = int((~good_mask).sum())

    if n_bad > 0:
        logging.info(
            "Loop Edit (SBE toolkit): removing %d/%d scans (%.1f%%).",
            n_bad, n, 100.0 * n_bad / n
        )

    df = df[good_mask].reset_index(drop=True)
    return df


def apply_tau_shift(df, config):
    """Aligns oxygen samples with T/C by advancing o2_umol_l in time.

    The SBE43 oxygen sensor responds more slowly than the T/C sensors.
    This step advances the O2 signal by ALIGN_OXY_SHIFT scans so that
    the oxygen value is aligned with the corresponding T/C measurement.

    Uses the Sea-Bird Scientific toolkit's `align_ctd`, which performs
    interpolation-based time shifting (more accurate than integer-sample
    rolling).  A positive offset_s pulls future values forward in time,
    which is equivalent to ``df['o2_umol_l'].shift(-shift_steps)``.

    FLAG_VALUE (-9.99e-29) is returned for edge samples that fall outside
    the valid range; these are replaced with NaN and then forward/back
    filled so downstream physics sees no gaps.

    Writes to ``o2_aligned`` so ``apply_physics`` can apply the
    density-normalization step on the already-aligned values.
    """
    if 'o2_umol_l' not in df.columns:
        logging.warning("No o2_umol_l column found; skipping oxygen tau shift.")
        return df

    shift_steps = int(float(config.get('ALIGN_OXY_SHIFT', 20.0)))
    dt = 1.0 / _sample_rate_hz(df)
    offset_s = shift_steps * dt    # positive = advance O2 forward in time

    o2_shifted = align_ctd(
        df['o2_umol_l'].to_numpy(),
        offset=offset_s,
        sample_interval=dt,
    )

    # Replace FLAG_VALUE edge sentinels with NaN, then fill
    o2_series = pd.Series(o2_shifted, index=df.index)
    o2_series = o2_series.where(o2_series != FLAG_VALUE).ffill().bfill()
    df['o2_aligned'] = o2_series
    return df


def apply_physics(df, config):
    """Computes SP, theta, density, derived oxygen, pH, chl, and smooths.

    Requires `apply_ctm_correction` and `apply_tau_shift` to have run.
    Recomputes Practical Salinity from the CTM-corrected conductivity
    via `gsw.SP_from_C` (cond in mS/cm = S/m * 10). Falls back to the
    SBE-precomputed `sal_raw` only if `gsw` is unavailable, with a
    warning.
    """
    t68_conv = float(config.get('T68_CONVERSION', 1.00024))
    tc = float(config.get('LPF_TIME_CONSTANT', 0.5))
    dt = 1.0 / _sample_rate_hz(df)
    sal_offset = float(config.get('SAL_OFFSET', 0.0))

    df['t68'] = df['temp_raw'] * t68_conv
    df['in_situ_temp'] = df['temp_raw']  # ITS-90; what the dashboard plots

    # --- Salinity from corrected conductivity ---
    if 'cond_final' in df.columns and _HAS_GSW:
        # gsw.SP_from_C wants conductivity in mS/cm; SBE c0S/m is S/m, so *10.
        # Temperature input is ITS-90, pressure is dbar.
        #
        # Matched filter: cond_final has already been low-pass filtered
        # (tc=LPF_TIME_CONSTANT) in apply_ctm_correction.  temp_raw receives
        # the same filter here before entering PSS-78 so T and C have identical
        # frequency response going into salinity derivation.
        temp_for_sp = low_pass_filter(df['temp_raw'].to_numpy(), tc, dt)

        # Align CTD: advance temperature to match the C cell sampling point.
        # On the SBE 19plus the T sensor leads the C cell during descent, so
        # at any scan C is sampling water that T already measured. SeaSoft's
        # Align CTD module corrects this by advancing T by ALIGN_TEMP_SHIFT
        # scans (AlignCTD.psa value = 0.5 s = 2 scans at 4 Hz).
        # Uses the same toolkit align_ctd function as the O2 tau shift.
        temp_shift = int(float(config.get('ALIGN_TEMP_SHIFT', 0)))
        if temp_shift > 0:
            offset_s = temp_shift * dt
            temp_for_sp = align_ctd(temp_for_sp, offset_s, dt)

        df['SP'] = gsw.SP_from_C(df['cond_final'].to_numpy() * 10.0,
                                  temp_for_sp,
                                  df['pres_raw'].to_numpy()) + sal_offset
    else:
        if not _HAS_GSW:
            logging.warning("gsw not available; falling back to SBE pre-computed sal_raw.")
        df['SP'] = df['sal_raw'] + sal_offset

    # Theta and density (EOS-80 expects IPTS-68 temp)
    df['theta'] = calculate_theta_eos80(df['SP'].to_numpy(),
                                        df['t68'].to_numpy(),
                                        df['pres_raw'].to_numpy())
    df['rho'] = calculate_density_eos80(df['SP'].to_numpy(),
                                        df['t68'].to_numpy(),
                                        df['pres_raw'].to_numpy())

    # Oxygen: convert tau-aligned umol/L to umol/kg using in-situ density.
    # If neither o2_aligned nor o2_umol_l exists (no O2 sensor), fill with NaN
    # rather than raising a KeyError.
    if 'o2_aligned' in df.columns:
        o2_source = df['o2_aligned']
    elif 'o2_umol_l' in df.columns:
        o2_source = df['o2_umol_l']
    else:
        logging.warning("No oxygen column found; o2_final will be NaN.")
        o2_source = np.nan
    df['o2_final'] = (o2_source * float(config.get('O2_BOOST_RATIO', 1.0))) / (df['rho'] / 1000.0)
    if 'ph_raw' in df.columns:
        df['ph_final'] = df['ph_raw'] + float(config.get('PH_DRIFT', 0.0))
    else:
        df['ph_final'] = np.nan
    df['chl_final'] = df['chl_raw'] if 'chl_raw' in df.columns else np.nan

    # Low-pass filter sensor-direct columns only (Sea-Bird toolkit LPF).
    # SP, rho, and theta are derived from inputs that are already filtered
    # (cond_final from apply_ctm_correction; temp_for_sp from line 387 above),
    # so a second LPF here would over-smooth them and wash out the CTM
    # correction — making SP match SeaSoft "Salinity 1st" instead of "2nd".
    # SeaSoft's Derive module does not re-filter its outputs.
    # in_situ_temp is temp_raw (unfiltered), so it gets one LPF here.
    # o2_final comes from a tau-shifted raw sensor signal and benefits from LPF.
    # interpolate/bfill/ffill first so edge NaNs don't propagate into the filter.
    targets = ['o2_final', 'in_situ_temp']
    for col in targets:
        if col not in df.columns:
            continue
        df[col] = df[col].interpolate(method='linear', limit_direction='both').bfill().ffill()
        if df[col].isna().all():
            continue  # nothing to smooth
        df[col] = low_pass_filter(df[col].to_numpy(), tc, dt)
    return df



def apply_chl_calibration(df, config):
    """Chlorophyll-specific calibration and low-pass filtering.

    Matches SBE Data Processing: applies calibration (slope/offset) then
    a single Sea-Bird Low-Pass Filter with LPF_TIME_CONSTANT (seconds).
    Uses the same time constant as all other derived variables.

    The 1-metre bin averaging in main.py provides additional low-pass
    filtering downstream.
    """
    slope = float(config.get('CHL_SLOPE', 1.0))
    offset = float(config.get('CHL_OFFSET', 0.0))
    tc = float(config.get('LPF_TIME_CONSTANT', 0.5))
    dt = 1.0 / _sample_rate_hz(df)

    if 'chl_raw' not in df.columns:
        return df

    df['chl_final'] = (df['chl_raw'] * slope) + offset
    df['chl_final'] = df['chl_final'].interpolate(method='linear', limit_direction='both').bfill().ffill()
    df['chl_final'] = low_pass_filter(df['chl_final'].to_numpy(), tc, dt)
    return df


def apply_wild_edit(df, config):
    """Two-pass block-statistics spike removal (Sea-Bird Wild Edit).

    Replicates SeaSoft's Wild Edit module: each variable is divided into
    blocks of WILD_SCANS_PER_BLOCK scans.  Within each block, values
    farther than WILD_STD_PASS1 standard deviations from the block mean
    are flagged on the first pass, then a second pass flags values farther
    than WILD_STD_PASS2 standard deviations from the remaining mean.
    Flagged samples are replaced with FLAG_VALUE (-9.99e-29).

    main.py replaces FLAG_VALUE with NaN before the pandas groupby so
    that outlier scans do not contaminate bin averages.

    Applied after apply_surgical_chl so spikes are identified in the
    calibrated, filtered variables.  applied before apply_qc_flags so
    velocity-flagged scans are not counted as spikes.
    """
    std1 = float(config.get('WILD_STD_PASS1', 2.0))
    std2 = float(config.get('WILD_STD_PASS2', 20.0))
    block = int(float(config.get('WILD_SCANS_PER_BLOCK', 100)))
    dist = float(config.get('WILD_DISTANCE_TO_MEAN', 0.0))

    targets = ['SP', 'theta', 'rho', 'o2_final', 'ph_final', 'chl_final', 'in_situ_temp']
    n = len(df)
    base_flags = np.zeros(n, dtype=float)

    for col in targets:
        if col not in df.columns:
            continue
        data = df[col].to_numpy(dtype=float, copy=True)
        # Wild edit expects no NaNs; temporarily replace NaN with FLAG_VALUE
        nan_mask = np.isnan(data)
        data[nan_mask] = FLAG_VALUE
        result = wild_edit(
            data,
            base_flags.copy(),
            std_pass_1=std1,
            std_pass_2=std2,
            scans_per_block=block,
            distance_to_mean=dist,
            exclude_bad_flags=False,
        )
        # Restore original NaNs (wild_edit may not have preserved them)
        result[nan_mask] = np.nan
        n_flagged = int(np.sum(result == FLAG_VALUE))
        if n_flagged > 0:
            logging.debug("Wild Edit: %s flagged %d/%d scans.", col, n_flagged, n)
        df[col] = result

    return df


def apply_qc_flags(df, config):
    """Velocity-based QC: 1 = good, 3 = high-velocity / bad."""
    threshold = float(config.get('QC_VELOCITY', 1.2))  # m/s; matches calibration.csv default
    rate_hz = _sample_rate_hz(df)
    df['qc_flag'] = 1
    mask_bad_velocity = (df['pres_raw'].diff().abs() * rate_hz) >= threshold
    df.loc[mask_bad_velocity, 'qc_flag'] = 3
    return df
