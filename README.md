# Western Flyer CTD Data

This repository contains the processing pipeline and visualization suite for Sea-Bird CTD data collected aboard the R/V Western Flyer.

## Feature Summary

The Western Flyer CTD codebase is a modular, automated processing pipeline designed for Sea-Bird SBE 19plus V2 CTD instrumentation. It provides a standardized framework for raw hex data conversion, physical correction, and interactive visualization.

* **Raw .hex Processing:** The pipeline starts from raw instrument `.hex` files rather than SeaSoft-processed `.cnv` files. Calibration coefficients are read directly from the instrument's `.xmlcon` configuration file, giving you full provenance over every conversion step.
* **Sea-Bird Scientific Toolkit:** All sensor conversion and correction routines delegate to the official [Sea-Bird Scientific Python toolkit](https://github.com/Sea-Bird-Scientific/seabirdscientific). This includes temperature/pressure/conductivity conversion, O2 tau and hysteresis corrections, cell-thermal-mass correction, Sea-Bird Low-Pass Filter, O2 tau-shift alignment, loop edit, Wild Edit spike removal, and downcast extraction.
* **Modular Cruise Architecture:** Each expedition is isolated under `cruises/<cruise_id>/`. Adding a new cruise requires only creating a directory with the raw `.hex` files and metadata; the directory name is the cruise ID.
* **Sidecar Data Integration:** The pipeline joins raw CTD measurements with `config.toml` (processing parameters), `calibration.xmlcon` (sensor coefficients), and `stations.csv` (station metadata) automatically.
* **Automated EOS-80 Physics Pipeline:** Core oceanographic correction logic runs in this fixed order: soak elimination → cell-thermal-mass (CTM) conductivity correction → O2 tau-shift alignment → loop edit → EOS-80 physics (Practical Salinity from CTM-corrected conductivity via PSS-78, then potential temperature and density) → chlorophyll calibration (`apply_chl_calibration`) → Wild Edit spike removal → velocity-based QC → downcast extraction → 1 m bin averaging.
* **Atomic, Idempotent Storage:** Processed data is managed via DuckDB. Per-cruise records are replaced atomically inside a single transaction. A `build_metadata` table records the SHA-256 of the `.xmlcon` file and the git commit for every build, so every row in `ctd_data` is traceable to the exact calibration that produced it.
* **Per-Build Audit Logging:** Each invocation writes a fresh, timestamped log file (`logs/wf_build_<cruise>_<timestamp>.log`).
* **Interactive Dashboard Suite:** A built-in visualization suite powered by the HoloViz ecosystem (Panel, HoloViews, Bokeh) renders interactive dashboards. Supports multi-variable vertical profiling, T-S analysis, vertical section plots with isopycnal (σθ) contours, AOU, MLD/stability, the metabolic index Φ, geolocation, and tabular export.

## Root Directory Notes

The project root can be named anything (e.g., `westernflyer_cruise_ctd`).

To avoid Windows permission issues, place the root folder inside your Windows user directory:
`C:\Users\<your_username>\westernflyer_cruise_ctd`

## Python Installation Requirements

This project requires **Python 3.12**.

### Windows
When installing, the installer will ask: "Add Python to PATH?" **Select NO.**

* **Typical path:** `C:\Users\<your_username>\AppData\Local\Programs\Python\Python312\`

### macOS (Homebrew)
    brew install python@3.12

### Linux (Ubuntu / Debian)
    sudo add-apt-repository ppa:deadsnakes/ppa
    sudo apt update
    sudo apt install python3.12 python3.12-venv

## Project Structure

*`logs/` and `processed/` are created automatically on first run.*

    .
    |   ctd_holoviews.py
    |   eos80_processing.py
    |   main.py
    |   pyproject.toml
    |   sbe19plus_ingestion.py
    +---cruises
    |   \---baja2025
    |       |   config.toml        processing parameters
    |       |   calibration.xmlcon     sensor calibration (from SeaSoft)
    |       |   stations.csv         station metadata
    |       |
    |       \---hex
    |               20250416_cast1.hex
    |               20250416_cast2.hex
    |               ... (etc)
    +---logs
    +---processed
    \---scripts
        +---linux_mac
        |       ctd_build.sh              (generic; takes <cruise_id>)
        |       ctd_build_baja2025.sh     (wrapper around ctd_build.sh baja2025)
        |       ctd_dashboard.sh
        |
        \---windows
                ctd_build.bat             (generic; takes <cruise_id>)
                ctd_build_baja2025.bat    (wrapper around ctd_build.bat baja2025)
                ctd_dashboard.bat

## File Descriptions

* **main.py** — Orchestrator. Parses CLI args, loads the `.xmlcon` calibration, runs the per-cast pipeline, writes the DuckDB atomically, and stamps build provenance.
* **sbe19plus_ingestion.py** — Hex ingestion and calibration. Reads raw `.hex` files via `read_hex_file`, applies temperature/pressure/conductivity/O2 (with tau and hysteresis corrections)/pH/CHL conversions from the Sea-Bird Scientific toolkit, and attaches cruise-log metadata.
* **eos80_processing.py** — EOS-80 physics pipeline. CTM correction, O2 tau-shift, loop edit, Sea-Bird Low-Pass Filter, Wild Edit spike removal, and depth/theta/density calculations all delegate to the Sea-Bird Scientific toolkit. Salinity is derived from CTM-corrected conductivity via `gsw.SP_from_C`.
* **ctd_holoviews.py** — Read-only Panel dashboard. Reactive widgets update all tabs when cruise or station selection changes.
* **pyproject.toml** — Package definition. Enables `pip install -e .` to install the project and all dependencies in one step.
* **config.toml** — Processing constants per cruise (CTM coefficients, bin size, QC thresholds, Wild Edit parameters). See schema below.
* **calibration.xmlcon** — Sensor calibration file exported from SeaSoft. Contains T/C/P/O2/pH/CHL coefficients. This is the Jan 2026 calibration for CTD S/N 8289. Replace with the current xmlcon before each cruise.
* **stations.csv** — Per-cast metadata: station name, lat/lon, cast numbers, start time.

## Prerequisites

Python 3.12 is required (see above). The **Sea-Bird Scientific Python Toolkit** is declared as a dependency in `pyproject.toml` and is installed automatically from PyPI when you run `pip install -e .` — no manual download required.

## Getting Started

    git clone https://github.com/joeacarlisle/westernflyer_cruise_ctd.git
    cd westernflyer_cruise_ctd

Then follow the **Prerequisites** and **Setup** sections below to install Python 3.12 and the project dependencies.

## Setup: Virtual Environment and Installation

*Run from the project root directory.*

### Windows (Command Prompt)

    "C:\Users\<your_username>\AppData\Local\Programs\Python\Python312\python.exe" -m venv .venv
    .venv\Scripts\activate
    pip install -e .

### macOS / Linux

    python3.12 -m venv .venv
    source .venv/bin/activate
    pip install -e .

`pip install -e .` reads `pyproject.toml` and installs all dependencies in one step, including the Sea-Bird Scientific toolkit pulled from PyPI.

The `-e` (editable) flag means source changes are picked up immediately without reinstalling.

## Preparing a New Cruise

1. Create a cruise directory: `cruises/<cruise_id>/`
2. Place raw `.hex` files in `cruises/<cruise_id>/hex/`
3. Copy the sensor calibration file from SeaSoft to `cruises/<cruise_id>/calibration.xmlcon`
4. Create `cruises/<cruise_id>/config.toml` (see schema below)
5. Create `cruises/<cruise_id>/stations.csv` (see schema below)

### stations.csv schema

| Column | Description |
|--------|-------------|
| `file_key` | Hex filename stem (e.g. `20250416_cast1` for `20250416_cast1.hex`) |
| `station_name` | Station identifier (e.g. `STA01`) |
| `lat` | Latitude in `DD-MM.MMN` format (e.g. `23-30.00N`) or decimal degrees |
| `lon` | Longitude in `DD-MM.MMW` format (e.g. `110-15.00W`) or decimal degrees |
| `sb_cast` | Sea-Bird instrument cast number |
| `wf_cast` | Western Flyer expedition cast number |
| `start_utc` | Cast start time fallback if not in hex header (`YYYY-MM-DD HH:MM:SS`) |

## Running the Pipeline

**Activate the virtual environment first.** Required every time you open a new terminal.

### Windows (Command Prompt)
    cd C:\Users\<your_username>\westernflyer_cruise_ctd
    .venv\Scripts\activate

You should see `(.venv)` at the start of your prompt.

### macOS / Linux
    cd /path/to/westernflyer_cruise_ctd
    source .venv/bin/activate

---

Run the build with the cruise ID as the argument:

* **Windows:** `scripts\windows\ctd_build.bat <cruise_id>`
* **Linux / macOS:** `./scripts/linux_mac/ctd_build.sh <cruise_id>`

Alternatively, call `main.py` directly for full control. The following are example invocations — run only the one that fits your needs:

    # Basic build
    python main.py baja2025

    # Override bin size or enable verbose logging
    python main.py baja2025 --bin-size 0.5 --verbose

    # Write output to a custom DuckDB path
    python main.py baja2025 --db /path/to/custom.duckdb

    # Override the xmlcon or hex directory locations
    python main.py baja2025 --xmlcon /path/to/custom.xmlcon
    python main.py baja2025 --hex-dir /path/to/hex_files/

### CLI arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `cruise_id` | *(required)* | Cruise directory name under `cruises/` |
| `--bin-size` | From `config.toml` | Vertical bin size in meters |
| `--xmlcon` | `cruises/<id>/calibration.xmlcon` | Path to sensor calibration file |
| `--hex-dir` | `cruises/<id>/hex/` | Directory containing raw `.hex` files |
| `--db` | `processed/wf_ctd_eos80.duckdb` | DuckDB output path |
| `--verbose` / `-v` | Off | Enable DEBUG-level logging |

*On macOS / Linux, make scripts executable once:*

    chmod +x scripts/linux_mac/ctd_build.sh
    chmod +x scripts/linux_mac/ctd_build_baja2025.sh
    chmod +x scripts/linux_mac/ctd_dashboard.sh

## Running the Dashboard

* **Windows:** `scripts\windows\ctd_dashboard.bat`
* **Linux / macOS:** `./scripts/linux_mac/ctd_dashboard.sh`

## config.toml Schema

`config.toml` uses TOML format with named sections (e.g. `[processing]`, `[wild_edit]`). The cruise ID is **not** stored here — it comes from the directory name.

### Calibration Offsets

| Parameter | Description |
|-----------|-------------|
| `CTD_SERIAL` | Instrument serial number (informational; not read by pipeline) |
| `SAL_OFFSET` | Additive offset applied to recomputed Practical Salinity (PSU). Use to match bottle salinity samples. |
| `PH_DRIFT` | Additive offset for pH. Use to correct for sensor drift between calibrations. |
| `O2_BOOST_RATIO` | Multiplicative gain applied to oxygen before density normalization. Use to match Winkler titration values. |
| `CHL_SLOPE` | Multiplicative secondary gain for chlorophyll (default 1.0). The xmlcon ScaleFactor is applied during ingestion; use this only when a post-deployment fluorometer calibration against extracted chlorophyll samples is available. |
| `CHL_OFFSET` | Additive offset for chlorophyll (µg/L). |

### SBE Processing Constants

| Parameter | Description |
|-----------|-------------|
| `T68_CONVERSION` | Scale factor converting ITS-90 temperature to IPTS-68 (default 1.00024). EOS-80 equations require IPTS-68. |
| `CTM_ALPHA` | Cell-thermal-mass amplitude α (default 0.04). SBE 19plus default. Passed to `cell_thermal_mass()`. |
| `CTM_TAU` | Cell-thermal-mass time constant τ in seconds (default 8.0). SBE 19plus default. Passed to `cell_thermal_mass()`. |
| `ALIGN_OXY_SHIFT` | Samples to advance oxygen data to correct for SBE 43 electrochemical lag. At 4 Hz, 20 samples = 5 s ≈ 5 m at 1 m/s. Passed to `align_ctd()`. |
| `LPF_TIME_CONSTANT` | Sea-Bird Low-Pass Filter time constant in seconds. Applied to conductivity before CTM, to temperature before PSS-78, and to all derived physics variables. Default 0.5 s (matches SeaSoft default for 4 Hz profiling). Passed to `low_pass_filter()`. |

### Wild Edit Parameters

| Parameter | Description |
|-----------|-------------|
| `WILD_STD_PASS1` | First-pass spike rejection threshold (standard deviations from block mean). Default 2.0. |
| `WILD_STD_PASS2` | Second-pass spike rejection threshold (standard deviations). Default 20.0. |
| `WILD_SCANS_PER_BLOCK` | Number of scans per statistics block. Default 100 (25 s at 4 Hz). |
| `WILD_DISTANCE_TO_MEAN` | Minimum absolute distance from block mean before a scan can be flagged (0 = no floor). Default 0.0. |

### Soak Detection Parameters

| Parameter | Description |
|-----------|-------------|
| `MIN_DEPTH_FLOOR` | Minimum depth (m) before soak detection begins. Default 10.0. |
| `LOOP_MIN_VELOCITY` | Minimum descent velocity (m/s) for soak detection and loop edit. Default 0.25. |
| `SOAK_WINDOW_SIZE` | Consecutive samples that must satisfy descent criteria to end the soak period. Default 20. |

### QC and Resolution Thresholds

| Parameter | Description |
|-----------|-------------|
| `QC_VELOCITY` | Descent speed (m/s) above which a scan is flagged QC=3. Default 1.2. |
| `BIN_SIZE_METERS` | Vertical bin resolution in the database (m). Default 1.0. Override with `--bin-size`. |

## Contact

For questions regarding this data archive or the Western Flyer Foundation, please contact **joeacarlisle@gmail.com**
