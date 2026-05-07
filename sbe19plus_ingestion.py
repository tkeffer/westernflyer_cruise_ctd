"""Sea-Bird SBE 19plus .hex file ingestion using the Sea-Bird Scientific toolkit.

Reads raw .hex files, converts to engineering units using calibration
coefficients parsed from an .xmlcon file, and attaches station metadata
from the cruise log.

Conversion pipeline per cast:
    read_hex_file         raw A/D counts per scan
    convert_temperature   ITS-90 °C (use_mv_r=True, required for SBE 19plus V2)
    convert_pressure      dbar (temperature compensation voltage from header)
    convert_conductivity  S/m (input is raw_int/256 from read_hex_file)
    SP_from_C             initial practical salinity for O2 calibration
    convert_sbe43_oxygen  ml/L (with tau and hysteresis corrections) → µmol/L × 44.660
    convert_sbe18_ph      pH (total scale)
    convert_eco           µg/L CHL
    attach metadata + timestamps

The returned DataFrame is identical in schema to the former .cnv version
so that eos80_processing.py and main.py require no changes.
"""

import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import gsw
import numpy as np
import pandas as pd

import seabirdscientific.cal_coefficients as cc
from seabirdscientific.conversion import (
    convert_temperature,
    convert_pressure,
    convert_conductivity,
    convert_sbe43_oxygen,
    convert_sbe18_ph,
    convert_eco,
)
from seabirdscientific.instrument_data import (
    InstrumentType,
    Sensors,
    read_hex_file,
)

# 1 mL/L O2 = 44.660 µmol/L  (standard SBE conversion factor)
_O2_MLPERL_TO_UMOLPERL = 44.660

# Sensors enabled on the Western Flyer CTD in profiling mode.
# Order must follow the Sensors enum ordering used by read_SBE19plus_format_0.
_PROFILING_SENSORS = [
    Sensors.Temperature,   # channel 0 — 6 hex chars, raw A/D counts
    Sensors.Conductivity,  # channel 1 — 6 hex chars, raw_int/256
    Sensors.Pressure,      # channel 2 — 6+4 hex chars (pressure + T comp)
    Sensors.ExtVolt0,      # volt 0 — SBE43 dissolved oxygen
    Sensors.ExtVolt1,      # volt 1 — SBE18 pH
    Sensors.ExtVolt2,      # volt 2 — WetLabs ECO CHL-a fluorometer
]


# ---------------------------------------------------------------------------
# xmlcon calibration parser
# ---------------------------------------------------------------------------

def parse_xmlcon(xmlcon_path):
    """Parse sensor calibration coefficients from a Sea-Bird .xmlcon file.

    Reads the Jan 2026 calibration for CTD S/N 8289.  Returns a dict
    with keys 'temp', 'cond', 'pres', 'o2', 'ph', 'chl', each holding
    the appropriate seabirdscientific.cal_coefficients dataclass.

    Raises ValueError if any required sensor block is absent.
    """
    xmlcon_path = Path(xmlcon_path)
    tree = ET.parse(xmlcon_path)
    root = tree.getroot()
    coefs = {}

    # --- Temperature (Steinhart-Hart, ITS-90) ---
    for sensor in root.iter('TemperatureSensor'):
        coefs['temp'] = cc.TemperatureCoefficients(
            a0=float(sensor.find('A0').text),
            a1=float(sensor.find('A1').text),
            a2=float(sensor.find('A2').text),
            a3=float(sensor.find('A3').text),
        )
        break   # only one T sensor

    # --- Conductivity (G/H/I/J, equation="1") ---
    for sensor in root.iter('ConductivitySensor'):
        eq1 = None
        for c in sensor.findall('Coefficients'):
            if c.get('equation') == '1':
                eq1 = c
                break
        if eq1 is None:
            raise ValueError(f"{xmlcon_path.name}: ConductivitySensor has no equation='1' block.")
        coefs['cond'] = cc.ConductivityCoefficients(
            g=float(eq1.find('G').text),
            h=float(eq1.find('H').text),
            i=float(eq1.find('I').text),
            j=float(eq1.find('J').text),
            cpcor=float(eq1.find('CPcor').text),
            ctcor=float(eq1.find('CTcor').text),
            wbotc=float(eq1.find('WBOTC').text),
        )
        break

    # --- Pressure (Strain-gauge, PA0/PA1/PA2 + PTCA/PTCB/PTEMPA) ---
    for sensor in root.iter('PressureSensor'):
        coefs['pres'] = cc.PressureCoefficients(
            pa0=float(sensor.find('PA0').text),
            pa1=float(sensor.find('PA1').text),
            pa2=float(sensor.find('PA2').text),
            ptca0=float(sensor.find('PTCA0').text),
            ptca1=float(sensor.find('PTCA1').text),
            ptca2=float(sensor.find('PTCA2').text),
            ptcb0=float(sensor.find('PTCB0').text),
            ptcb1=float(sensor.find('PTCB1').text),
            ptcb2=float(sensor.find('PTCB2').text),
            ptempa0=float(sensor.find('PTEMPA0').text),
            ptempa1=float(sensor.find('PTEMPA1').text),
            ptempa2=float(sensor.find('PTEMPA2').text),
        )
        break

    # --- Oxygen SBE43 (equation="1", 2007 Sea-Bird form) ---
    for sensor in root.iter('OxygenSensor'):
        eq1 = None
        for c in sensor.findall('CalibrationCoefficients'):
            if c.get('equation') == '1':
                eq1 = c
                break
        if eq1 is None:
            raise ValueError(f"{xmlcon_path.name}: OxygenSensor has no equation='1' block.")
        coefs['o2'] = cc.Oxygen43Coefficients(
            soc=float(eq1.find('Soc').text),
            v_offset=float(eq1.find('offset').text),
            tau_20=float(eq1.find('Tau20').text),
            a=float(eq1.find('A').text),
            b=float(eq1.find('B').text),
            c=float(eq1.find('C').text),
            e=float(eq1.find('E').text),
            d0=float(eq1.find('D0').text),
            d1=float(eq1.find('D1').text),
            d2=float(eq1.find('D2').text),
            h1=float(eq1.find('H1').text),
            h2=float(eq1.find('H2').text),
            h3=float(eq1.find('H3').text),
        )
        break

    # --- pH SBE18 (slope / offset) ---
    for sensor in root.iter('pH_Sensor'):
        coefs['ph'] = cc.PH18Coefficients(
            slope=float(sensor.find('Slope').text),
            offset=float(sensor.find('Offset').text),
        )
        break

    # --- Chlorophyll WetLabs ECO (scale factor / Vblank) ---
    for sensor in root.iter('FluoroWetlabECO_AFL_FL_Sensor'):
        coefs['chl'] = cc.ECOCoefficients(
            slope=float(sensor.find('ScaleFactor').text),
            offset=float(sensor.find('Vblank').text),
        )
        break

    missing = [k for k in ('temp', 'cond', 'pres', 'o2', 'ph', 'chl') if k not in coefs]
    if missing:
        raise ValueError(f"parse_xmlcon: {xmlcon_path.name} missing sensors: {missing}")

    logging.info(
        "Loaded calibration from %s: T/C/P/O2/pH/CHL", xmlcon_path.name
    )
    return coefs


# ---------------------------------------------------------------------------
# Hex header parser
# ---------------------------------------------------------------------------

def parse_hex_header(hex_path):
    """Parse the ASCII header block of an SBE 19plus .hex file.

    Returns
    -------
    start_time : pd.Timestamp or None
        Cast start time from the ``* cast`` line.  None if not found.
    sample_interval_s : float
        Sample interval in seconds (default 0.25 for SBE 19plus profiling).

    Notes
    -----
    The voltage board calibration lines ("* volt N: offset = X, slope = Y")
    are intentionally NOT parsed here.  The Sea-Bird toolkit conversion
    functions (convert_sbe43_oxygen, convert_sbe18_ph, convert_eco) expect
    the raw A/D voltages produced by read_hex_file.  The sensor calibration
    coefficients in the .xmlcon (Soc/v_offset for O2, Slope/Offset for pH,
    ScaleFactor/Vblank for CHL) were derived against raw A/D voltages, so
    applying the board correction as a separate multiplicative step would
    introduce ~2-3 pH unit error and incorrect O2/CHL values.
    SeaSoft's datcnv incorporates the board correction internally within
    each sensor equation — it is not a standalone pre-processing step.
    """
    hex_path = Path(hex_path)
    start_time = None
    sample_interval_s = 0.25     # SBE 19plus profiling-mode default (4 Hz)

    with open(hex_path, 'r', errors='replace') as f:
        for line in f:
            line_s = line.rstrip()

            # Cast timestamp: "* cast  N  DD Mon YYYY HH:MM:SS samples ..."
            if line_s.startswith('* cast') and start_time is None:
                m = re.search(
                    r'(\d{1,2}\s+\w{3}\s+\d{4}\s+\d{2}:\d{2}:\d{2})',
                    line_s,
                )
                if m:
                    start_time = pd.to_datetime(
                        m.group(1), format='%d %b %Y %H:%M:%S', errors='coerce'
                    )

            # Sample interval (various possible formats)
            if re.search(r'sample.interval', line_s, re.I):
                m = re.search(r':\s*([\d.]+)', line_s)
                if m:
                    try:
                        sample_interval_s = float(m.group(1))
                    except ValueError:
                        pass

            if line_s.startswith('*END*'):
                break

    return start_time, sample_interval_s


# ---------------------------------------------------------------------------
# Shared helpers (unchanged from .cnv version)
# ---------------------------------------------------------------------------

def dm_to_decimal(dm_str):
    """Converts DD-MM.MM[NSEW] compass format to decimal degrees.

    Returns None and logs a warning on unparseable input rather than
    raising, so a single bad cruise-log row doesn't kill the whole build.
    """
    if dm_str is None:
        return None
    if isinstance(dm_str, (int, float)):
        return float(dm_str)
    try:
        match = re.match(r"(\d+)-(\d+\.?\d*)([NSEW])", str(dm_str).strip())
        if not match:
            logging.warning(
                "Coordinate '%s' did not match DD-MM.MM[NSEW]; returning None.", dm_str
            )
            return None
        deg, mins, direction = match.groups()
        decimal = float(deg) + (float(mins) / 60.0)
        if direction in ('S', 'W'):
            decimal *= -1
        return decimal
    except (ValueError, AttributeError, TypeError) as exc:
        logging.warning("Failed to parse coordinate '%s': %s", dm_str, exc)
        return None


def load_config(file_path):
    """Loads cruise processing parameters from a TOML config file.

    Sections (e.g. [processing], [wild_edit]) are flattened into a single
    dict so all existing ``config.get('KEY', default)`` calls work unchanged.
    Duplicate keys across sections raise rather than silently overwriting.
    """
    import tomllib
    file_path = Path(file_path)
    logging.info("Loading configuration: %s", file_path)
    with open(file_path, 'rb') as f:
        raw = tomllib.load(f)
    flat = {}
    for section, value in raw.items():
        if isinstance(value, dict):
            for k, v in value.items():
                if k in flat:
                    raise ValueError(f"Duplicate key '{k}' in {file_path}")
                flat[k] = str(v)
        else:
            flat[section] = str(value)
    return flat


def load_cruise_log(file_path):
    """Loads cruise metadata CSV keyed by file_key (hex stem)."""
    logging.info("Loading cruise log: %s", file_path)
    df = pd.read_csv(file_path)
    if df['file_key'].duplicated().any():
        dups = df.loc[df['file_key'].duplicated(), 'file_key'].tolist()
        raise ValueError(f"Duplicate file_key entries in {file_path}: {dups}")
    return df.set_index('file_key').to_dict('index')


# ---------------------------------------------------------------------------
# Core ingestion function
# ---------------------------------------------------------------------------

def ingest_sbe_hex(file_path, cruise_log, config, sensor_coefs):
    """Convert one SBE 19plus V2 .hex file to a calibrated engineering-units DataFrame.

    Parameters
    ----------
    file_path : str or Path
        Path to the raw .hex file.
    cruise_log : dict
        Keyed by file_key (hex stem).  Each value has lat, lon,
        station_name, sb_cast, wf_cast, start_utc.
    config : dict
        calibration.csv parameters (CRUISE_ID and processing knobs).
    sensor_coefs : dict
        Output of ``parse_xmlcon()`` — the Jan 2026 coefficient objects.

    Returns
    -------
    pd.DataFrame
        One row per scan.  Columns:

        temp_raw    — temperature, ITS-90 °C
        cond_raw    — conductivity, S/m
        pres_raw    — pressure, dbar
        time_elapsed — elapsed seconds from cast start
        o2_umol_l  — dissolved oxygen, µmol/L
        ph_raw      — pH (total scale)
        chl_raw     — chlorophyll-a fluorescence, µg/L
        lat, lon    — decimal degrees
        station_name, sb_cast, wf_cast, station_id
        time_iso    — UTC timestamp (pd.Timestamp)
        cruise_id

        ``df.attrs['sample_interval_s']`` carries the per-cast sample
        interval (seconds) for use by downstream rate-dependent functions.
    """
    file_path = Path(file_path)
    file_stem = file_path.stem

    # --- Cruise-log lookup with cruise-median fallback ---
    valid_lats = [dm_to_decimal(v.get('lat')) for v in cruise_log.values()]
    valid_lons = [dm_to_decimal(v.get('lon')) for v in cruise_log.values()]
    valid_lats = [x for x in valid_lats if x is not None]
    valid_lons = [x for x in valid_lons if x is not None]
    fallback_lat = float(pd.Series(valid_lats).median()) if valid_lats else float('nan')
    fallback_lon = float(pd.Series(valid_lons).median()) if valid_lons else float('nan')

    default_meta = {
        'lat': fallback_lat, 'lon': fallback_lon, 'station_name': "Unknown",
        'sb_cast': 0, 'wf_cast': 0, 'start_utc': "2025-01-01 00:00:00",
    }
    meta = cruise_log.get(file_stem, default_meta)
    if file_stem not in cruise_log:
        logging.warning(
            "No cruise_log entry for %s; using cruise median coords "
            "(%.4f, %.4f).", file_stem, fallback_lat, fallback_lon
        )

    lat = dm_to_decimal(meta.get('lat', default_meta['lat']))
    lon = dm_to_decimal(meta.get('lon', default_meta['lon']))
    lat = lat if lat is not None else fallback_lat
    lon = lon if lon is not None else fallback_lon

    # --- Parse hex header ---
    start_time, sample_interval_s = parse_hex_header(file_path)

    if start_time is None or pd.isnull(start_time):
        start_time = pd.to_datetime(meta.get('start_utc', default_meta['start_utc']))
        logging.warning(
            "%s: no cast start_time in hex header; using cruise_log start_utc.", file_path.name
        )

    # --- Read raw hex data via Sea-Bird Scientific toolkit ---
    raw = read_hex_file(
        file_path,
        InstrumentType.SBE19Plus,
        _PROFILING_SENSORS,
    )

    if raw is None or raw.empty:
        raise ValueError(f"{file_path.name}: read_hex_file returned no data.")

    n = len(raw)

    # --- Temperature (ITS-90 °C) ---
    # use_mv_r=True: SBE 19plus V2 stores T as mV-converted-to-resistance counts
    temp_c = convert_temperature(
        raw['temperature'].to_numpy(),
        sensor_coefs['temp'],
        use_mv_r=True,
    )

    # --- Pressure (dbar) ---
    # 'temperature compensation' column from read_hex_file is already in volts
    # (instrument_data divides by COUNTS_TO_VOLTS = 13107)
    pres_dbar = convert_pressure(
        raw['pressure'].to_numpy(),
        raw['temperature compensation'].to_numpy(),
        sensor_coefs['pres'],
        units='dbar',
    )

    # --- Conductivity (S/m) ---
    # 'conductivity' column from read_hex_file is raw_int/256;
    # convert_conductivity further divides by 1000 to get kHz frequency
    cond_sm = convert_conductivity(
        raw['conductivity'].to_numpy(),
        temp_c,
        pres_dbar,
        sensor_coefs['cond'],
    )

    # --- Raw A/D voltages (no separate board calibration) ---
    # read_hex_file returns volt N columns already in volts (raw_int / 13107).
    # The toolkit conversion functions expect these raw A/D voltages directly.
    # The sensor calibration coefficients (Soc/v_offset for O2, Slope/Offset
    # for pH, ScaleFactor/Vblank for CHL) were derived against raw A/D volts,
    # so NO separate board-calibration step is applied here.
    v0 = raw['volt 0'].to_numpy()   # SBE43 dissolved O2
    v1 = raw['volt 1'].to_numpy()   # SBE18 pH
    v2 = raw['volt 2'].to_numpy()   # WetLabs ECO CHL

    # --- Initial practical salinity (needed for O2 conversion) ---
    # gsw.SP_from_C expects conductivity in mS/cm = S/m × 10
    sp_init = gsw.SP_from_C(cond_sm * 10.0, temp_c, pres_dbar)

    # --- Dissolved oxygen (µmol/L) ---
    # Tau correction removes the electrochemical lag of the SBE 43 membrane
    # (matches SeaSoft datcnv_ox_tau_correction = yes).
    # Hysteresis correction removes pressure-cycle hysteresis in the oxygen
    # signal (matches SeaSoft datcnv_ox_hysteresis_correction = yes).
    # Both require sample_interval_s, which is read from the hex header.
    # convert_sbe43_oxygen returns ml/L; multiply by 44.660 for µmol/L.
    o2_mll = convert_sbe43_oxygen(
        v0, temp_c, pres_dbar, sp_init, sensor_coefs['o2'],
        apply_tau_correction=True,
        apply_hysteresis_correction=True,
        sample_interval=sample_interval_s,
    )
    o2_umol_l = o2_mll * _O2_MLPERL_TO_UMOLPERL

    # --- pH (total scale) ---
    ph_raw = convert_sbe18_ph(v1, temp_c, sensor_coefs['ph'])

    # --- Chlorophyll-a fluorescence (µg/L) ---
    # convert_eco: slope × (raw − offset) = ScaleFactor × (V − Vblank)
    chl_raw = convert_eco(v2, sensor_coefs['chl'])

    # --- Assemble DataFrame ---
    df = pd.DataFrame({
        'temp_raw':   temp_c,
        'cond_raw':   cond_sm,
        'pres_raw':   pres_dbar,
        'o2_umol_l':  o2_umol_l,
        'ph_raw':     ph_raw,
        'chl_raw':    chl_raw,
    })

    # Uniform elapsed-time grid (scan 0 = cast start)
    df['time_elapsed'] = np.arange(n, dtype=float) * sample_interval_s

    # ISO timestamps
    df['time_iso'] = start_time + pd.to_timedelta(df['time_elapsed'], unit='s')

    # Station / cast metadata from cruise log
    df['lat'] = lat
    df['lon'] = lon
    df['station_name'] = meta.get('station_name', 'Unknown')
    df['sb_cast'] = meta.get('sb_cast', 0)
    df['wf_cast'] = meta.get('wf_cast', 0)
    df['station_id'] = (
        f"{df['station_name'].iloc[0]}_Cast{df['wf_cast'].iloc[0]}"
    )
    df['cruise_id'] = config.get('CRUISE_ID', 'unknown')

    # Stash sample interval for rate-dependent downstream functions
    df.attrs['sample_interval_s'] = sample_interval_s

    logging.info(
        "Ingested %s: %d scans, T_mean=%.2f°C, P_max=%.1f dbar, dt=%.3fs",
        file_path.name, n, float(np.nanmean(temp_c)),
        float(np.nanmax(pres_dbar)), sample_interval_s,
    )
    return df
