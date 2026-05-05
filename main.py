"""Western Flyer CTD build orchestrator.

Usage:
    python main.py <cruise_id> [--bin-size 1.0] [--xmlcon PATH]
                               [--hex-dir PATH] [--db PATH] [--verbose]

The cruise_id maps to ``cruises/<cruise_id>/`` and is the authoritative
source of truth for the cruise.  Expected layout:

    cruises/<cruise_id>/
        calibration.csv       processing parameters (CTM, bin size, etc.)
        calibration.xmlcon    sensor calibration coefficients (from SeaSoft)
        cruise_log.csv        per-cast metadata (lat, lon, station, etc.)
        hex/                  raw .hex files — one per cast

Processing pipeline (per cast):
    read_hex_file             raw A/D counts via Sea-Bird Scientific toolkit
    convert_temperature       ITS-90 °C  (use_mv_r=True for SBE 19plus V2)
    convert_pressure          dbar
    convert_conductivity      S/m
    convert_sbe43_oxygen      ml/L → µmol/L
    convert_sbe18_ph          pH
    convert_eco               µg/L CHL
    apply_soak_elimination    flag surface-soak rows
    apply_ctm_correction      cell thermal mass (Sea-Bird Scientific)
    apply_tau_shift           O2 lag alignment (Sea-Bird Scientific)
    apply_loop_edit           pressure-reversal removal (Sea-Bird Scientific)
    apply_physics             SP, theta, rho, o2_final
    apply_chl_calibration     CHL calibration + filter
    apply_wild_edit           two-pass block-statistics spike removal
    apply_qc_flags            velocity-based QC

The pipeline writes into a single DuckDB file:
    processed/wf_ctd_eos80.duckdb
        ctd_data        - one row per (station, cast, depth bin)
        build_metadata  - one row per build, with SHA-256 of .xmlcon

Per-cruise records are replaced (not appended) inside a single
transaction, so a failed insert can't leave the database empty.
"""

import argparse
import datetime as dt
import hashlib
import logging
import pathlib
import subprocess
import sys

import duckdb
import numpy as np
import pandas as pd

import eos80_processing
import sbe19plus_ingestion
from sbe19plus_ingestion import parse_xmlcon, ingest_sbe_hex
from seabirdscientific.processing import FLAG_VALUE, get_downcast

BASE_DIR = pathlib.Path(__file__).resolve().parent
DEFAULT_DB_PATH = BASE_DIR / "processed" / "wf_ctd_eos80.duckdb"
LOG_DIR = BASE_DIR / "logs"


def setup_logging(cruise_id, verbose=False):
    """Sets up logging with a per-build, timestamped log file."""
    LOG_DIR.mkdir(exist_ok=True)
    timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = LOG_DIR / f"wf_build_{cruise_id}_{timestamp}.log"
    level = logging.DEBUG if verbose else logging.INFO
    root = logging.getLogger()
    for h in list(root.handlers):
        root.removeHandler(h)
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path, mode='w'),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return log_path


def file_sha256(path):
    """Returns the hex SHA-256 of a file (used for calibration provenance)."""
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def git_commit():
    """Returns the current git commit hash, or 'unknown' if not a repo."""
    try:
        out = subprocess.check_output(
            ['git', '-C', str(BASE_DIR), 'rev-parse', 'HEAD'],
            stderr=subprocess.DEVNULL,
        )
        return out.decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return 'unknown'


def process_one_cast(file_path, cruise_log, config, cruise_id, output_res, sensor_coefs):
    """Runs the full per-cast pipeline and returns the binned DataFrame."""
    df = ingest_sbe_hex(file_path, cruise_log=cruise_log, config=config,
                        sensor_coefs=sensor_coefs)

    df = eos80_processing.apply_soak_elimination(df, config)
    df = eos80_processing.apply_ctm_correction(df, config)
    df = eos80_processing.apply_tau_shift(df, config)
    df = eos80_processing.apply_loop_edit(df, config)
    df = eos80_processing.apply_physics(df, config)
    df = eos80_processing.apply_chl_calibration(df, config)
    df = eos80_processing.apply_wild_edit(df, config)
    df = eos80_processing.apply_qc_flags(df, config)

    df = get_downcast(df, 'pres_raw')

    # Replace Wild Edit sentinel values with NaN so they don't contaminate
    # bin averages in the groupby below.
    df = df.replace(FLAG_VALUE, np.nan)

    df['dbar_bin'] = (df['pres_raw'] / output_res).round() * output_res
    df['depth_m'] = eos80_processing.calculate_depth_eos80(
        df['dbar_bin'].to_numpy(), df['lat'].to_numpy()
    )

    agg_config = {
        'time_iso': 'min',
        'rho': 'mean', 'SP': 'mean', 'theta': 'mean',
        'o2_final': 'mean', 'ph_final': 'mean', 'chl_final': 'mean',
        'in_situ_temp': 'mean',
        'lat': 'first', 'lon': 'first',
        'qc_flag': 'max',
    }
    agg_config = {k: v for k, v in agg_config.items() if k in df.columns}

    # is_soak is a groupby key so soak and non-soak scans bin into separate
    # rows.  This prevents soak artifacts from contaminating non-soak bin
    # averages while keeping soak data in DuckDB for the dashboard toggle.
    groupby_cols = ['station_id', 'station_name', 'sb_cast', 'wf_cast', 'dbar_bin', 'depth_m']
    if 'is_soak' in df.columns:
        groupby_cols.append('is_soak')

    df_binned = (df.groupby(groupby_cols).agg(agg_config).reset_index())

    df_binned['cruise_id'] = cruise_id
    df_binned['filename'] = file_path.name
    return df_binned


def commit_to_duckdb(full_df, cruise_id, db_path, calibration_path, build_started):
    """Replaces the cruise's rows atomically and records build provenance."""
    db_path.parent.mkdir(exist_ok=True)
    cal_hash = file_sha256(calibration_path)
    commit = git_commit()

    with duckdb.connect(str(db_path)) as con:
        con.register('df_view', full_df)
        con.execute("CREATE TABLE IF NOT EXISTS ctd_data AS SELECT * FROM df_view WHERE 1=0")
        con.execute("""
            CREATE TABLE IF NOT EXISTS build_metadata (
                cruise_id        VARCHAR,
                build_started_ts TIMESTAMP,
                build_finished_ts TIMESTAMP,
                calibration_sha256 VARCHAR,
                git_commit       VARCHAR,
                rows_written     BIGINT
            )
        """)
        con.execute("BEGIN TRANSACTION")
        try:
            con.execute("DELETE FROM ctd_data WHERE cruise_id = ?", [cruise_id])
            con.execute("INSERT INTO ctd_data SELECT * FROM df_view")
            con.execute("DELETE FROM build_metadata WHERE cruise_id = ?", [cruise_id])
            con.execute(
                "INSERT INTO build_metadata VALUES (?, ?, ?, ?, ?, ?)",
                [cruise_id, build_started, dt.datetime.now(), cal_hash, commit, len(full_df)],
            )
            con.execute("COMMIT")
        except Exception:
            con.execute("ROLLBACK")
            raise


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Build the Western Flyer CTD DuckDB for one cruise.")
    parser.add_argument("cruise_id", help="Cruise directory name under cruises/, e.g. baja2025")
    parser.add_argument("--bin-size", type=float, default=None,
                        help="Override BIN_SIZE_METERS from calibration.csv (m)")
    parser.add_argument("--db", type=pathlib.Path, default=DEFAULT_DB_PATH,
                        help=f"DuckDB output path (default: {DEFAULT_DB_PATH})")
    parser.add_argument("--xmlcon", type=pathlib.Path, default=None,
                        help="Path to .xmlcon calibration file "
                             "(default: cruises/<cruise_id>/calibration.xmlcon)")
    parser.add_argument("--hex-dir", type=pathlib.Path, default=None,
                        help="Directory containing raw .hex files "
                             "(default: cruises/<cruise_id>/hex/)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable DEBUG-level logging")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    cruise_id = args.cruise_id
    cruise_dir = BASE_DIR / "cruises" / cruise_id

    (BASE_DIR / "processed").mkdir(exist_ok=True)
    log_path = setup_logging(cruise_id, verbose=args.verbose)

    if not cruise_dir.exists():
        logging.error(f"Cruise directory '{cruise_dir}' not found.")
        sys.exit(1)

    build_started = dt.datetime.now()
    logging.info(f"--- STARTING PIPELINE FOR {cruise_id.upper()} ---")
    logging.info(f"Build log: {log_path}")

    calibration_path = cruise_dir / 'calibration.csv'
    config = sbe19plus_ingestion.load_config_csv(calibration_path)
    cruise_log = sbe19plus_ingestion.load_cruise_log(cruise_dir / 'cruise_log.csv')

    config['CRUISE_ID'] = cruise_id
    eos80_processing.validate_config(config)

    output_res = float(args.bin_size if args.bin_size is not None
                       else config.get('BIN_SIZE_METERS', 1.0))

    # Resolve .xmlcon path (default: cruises/<cruise_id>/calibration.xmlcon)
    xmlcon_path = args.xmlcon if args.xmlcon else cruise_dir / 'calibration.xmlcon'
    if not xmlcon_path.exists():
        logging.error(f"Calibration .xmlcon not found: {xmlcon_path}")
        sys.exit(1)
    sensor_coefs = parse_xmlcon(xmlcon_path)

    # Resolve hex directory (default: cruises/<cruise_id>/hex/)
    hex_dir = args.hex_dir if args.hex_dir else cruise_dir / 'hex'
    if not hex_dir.exists():
        logging.error(f"Hex directory not found: {hex_dir}")
        sys.exit(1)

    files = sorted(hex_dir.glob("*.hex"))
    if not files:
        logging.warning(f"No .hex files found in {hex_dir}")

    processed_data = []
    failures = []
    for file_path in files:
        logging.info(f"Processing: {file_path.name}")
        try:
            df_binned = process_one_cast(
                file_path, cruise_log, config, cruise_id, output_res, sensor_coefs
            )
            processed_data.append(df_binned)
        except Exception as exc:
            logging.exception(f"Failed to process {file_path.name}: {exc}")
            failures.append(file_path.name)
            continue

    if not processed_data:
        logging.warning("No data processed; database not updated.")
        if failures:
            logging.warning(f"All {len(failures)} casts failed: {failures}")
        sys.exit(1 if failures else 0)

    full_df = pd.concat(processed_data, ignore_index=True)
    logging.info(f"Committing {len(full_df)} binned records to {args.db}")
    commit_to_duckdb(full_df, cruise_id, args.db, xmlcon_path, build_started)

    if failures:
        logging.warning(f"Pipeline finished with {len(failures)} failed casts: {failures}")
    else:
        logging.info("Pipeline finished successfully.")


if __name__ == "__main__":
    main()
