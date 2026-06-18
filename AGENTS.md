# Project Components and Architecture

This document describes the functional components that comprise the Western Flyer CTD processing pipeline and provides guidelines for AI agents interacting with this repository.

## Functional Components

The codebase is organized into modular components, each responsible for a distinct stage of the data lifecycle.

### 1. Ingestion Component (`sbe19plus_ingestion.py`)
**Role:** Raw Data Interface & Calibration
- **Responsibilities:**
    - Parses Sea-Bird SBE 19plus V2 `.hex` files (raw A/D counts).
    - Extracts sensor calibration coefficients from `.xmlcon` files.
    - Loads cruise-specific station metadata from `stations.csv`.
    - Converts raw counts to initial engineering units (T, C, P, O2, pH, CHL) using the official Sea-Bird Scientific Python toolkit.
- **Key Interface:** `ingest_sbe_hex()`

### 2. Physics Component (`eos80_processing.py`)
**Role:** Oceanographic Correction & Property Derivation
- **Responsibilities:**
    - Performs physical corrections: Soak elimination, Cell Thermal Mass (CTM) correction, O2 tau-shift alignment, and Loop Edit (pressure reversal removal).
    - Derives EOS-80 physical properties: Practical Salinity (PSS-78), potential temperature ($\theta$), and in-situ density ($\rho$).
    - Applies secondary calibrations and quality control (QC) flags based on descent velocity and statistical spike removal (Wild Edit).
- **Key Interface:** `apply_physics()` and the `apply_*` pipeline functions.

### 3. Orchestrator Component (`main.py`)
**Role:** Workflow Manager & Persistence
- **Responsibilities:**
    - Orchestrates the end-to-end pipeline for a specific cruise.
    - Manages CLI arguments, logging, and directory resolution.
    - Ensures atomic and idempotent database commits to DuckDB.
    - Records build provenance (git commit hash, calibration file SHA-256, and timestamps) for auditability.
- **Key Interface:** `main()` entry point and `process_one_cast()`.

### 4. Analysis Component (`ctd_holoviews.py`)
**Role:** Interactive Visualization & Export
- **Responsibilities:**
    - Serves a Panel-based dashboard for interactive data exploration.
    - Renders vertical profiles, T-S analysis, and vertical sections with isopycnal contours.
    - Calculates real-time analytics: Apparent Oxygen Utilization (AOU), Mixed Layer Depth (MLD), and the Metabolic Index ($\Phi$).
    - Facilitates tabular data export for external analysis.
- **Key Interface:** `dashboard.servable()`

### 5. Infrastructure Component (Scripts)
**Role:** Environment-Specific Automation
- **Linux/macOS (`scripts/linux_mac/`):** Bash scripts for building the database and launching the dashboard.
- **Windows (`scripts/windows/`):** Batch files for equivalent operations in a CMD environment.

---

## AI Agent Interaction Guidelines

When an AI agent (like Junie or others) interacts with this repository, it should adhere to the following operational patterns:

### Adding a New Cruise
1. **Directory Structure:** Create `cruises/<cruise_id>/`.
2. **Raw Data:** Place `.hex` files in `cruises/<cruise_id>/hex/`.
3. **Calibration:** Ensure `calibration.xmlcon` and `config.toml` are present in the cruise root.
4. **Metadata:** Populate `stations.csv` with mappings from file keys to geographic coordinates and cast numbers.

### Modifying the Pipeline
- **Physics Changes:** All oceanographic logic should reside in `eos80_processing.py`. Use the `seabirdscientific` toolkit wherever possible to maintain parity with official Sea-Bird standards.
- **Ingestion Changes:** Add new sensors or change A/D conversion logic in `sbe19plus_ingestion.py`.
- **Database Schema:** If adding columns, update the aggregation configuration in `main.py` and the query logic in `ctd_holoviews.py`.

### Execution and Verification
- **Run the Build:** Use `python main.py <cruise_id>` to process data.
- **Logs:** Monitor `logs/wf_build_<cruise_id>_<timestamp>.log` for detailed processing diagnostics.
- **Database:** Verify output in `processed/wf_ctd_eos80.duckdb`. The `build_metadata` table should be checked to confirm successful persistence and provenance.

### Dependencies
- This project strictly requires **Python 3.12**.
- Other dependencies are listed in `pyproject.toml`.
