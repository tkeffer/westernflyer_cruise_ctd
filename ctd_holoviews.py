import pandas as pd
import numpy as np
import duckdb
import holoviews as hv
import panel as pn
import geoviews as gv
import geoviews.tile_sources as gvts
import gsw
import cartopy.crs as ccrs
import pathlib
import re
from io import BytesIO

# 1. ENGINE INITIALIZATION
pn.extension('tabulator')
hv.extension('bokeh')
gv.extension('bokeh')

# 2. DATABASE CONNECTION
BASE_DIR = pathlib.Path(__file__).resolve().parent
DB_PATH = BASE_DIR / "processed" / "wf_ctd_eos80.duckdb"
TABLE_NAME = "ctd_data"

if not DB_PATH.exists():
    raise FileNotFoundError(f"Database not found at {DB_PATH}. Run main.py first.")

con = duckdb.connect(str(DB_PATH), read_only=True)

# Fetch Metadata
raw_cruises = con.execute(f"SELECT DISTINCT cruise_id FROM {TABLE_NAME}").df()['cruise_id'].tolist()

# Initialize default stations for the first cruise
default_cruise = raw_cruises[0]
raw_stations = con.execute(f"SELECT DISTINCT station_id FROM {TABLE_NAME} WHERE cruise_id = ?", (default_cruise,)).df()['station_id'].tolist()
stations = sorted(raw_stations, key=lambda x: int(re.findall(r'\d+', x)[-1]) if re.findall(r'\d+', x) else 0)

# 3. GLOBAL WIDGETS
cruise_select = pn.widgets.Select(name='Cruise ID', options=raw_cruises, value=default_cruise)
station_select = pn.widgets.Select(name='Station ID', options=stations, value=stations[0] if stations else None)
depth_slider = pn.widgets.RangeSlider(name='Depth Range (m)', start=0, end=1000, value=(0, 600), step=1.0)
qc_checkbox = pn.widgets.Checkbox(name='Filter QC (Flag < 3)', value=False)
soak_toggle = pn.widgets.Checkbox(name='Show Soak Data (Surface)', value=False)

# Callback to update stations when cruise changes
def update_stations(event):
    new_cruise = event.new
    new_stations = con.execute(f"SELECT DISTINCT station_id FROM {TABLE_NAME} WHERE cruise_id = ?", (new_cruise,)).df()['station_id'].tolist()
    sorted_stations = sorted(new_stations, key=lambda x: int(re.findall(r'\d+', x)[-1]) if re.findall(r'\d+', x) else 0)
    station_select.options = sorted_stations
    if sorted_stations:
        station_select.value = sorted_stations[0]

cruise_select.param.watch(update_stations, 'value')

# 4. DATA LOGIC
@pn.depends(cruise_select, station_select, depth_slider, qc_checkbox, soak_toggle)
def get_clean_df(target_cruise, target_id, z_range, filter_qc, show_soak):
    if not target_id: return pd.DataFrame()
    qc_clause = "AND qc_flag < 3" if filter_qc else ""
    soak_clause = "" if show_soak else "AND is_soak = 0"
    query = f"""
        SELECT * FROM {TABLE_NAME}
        WHERE cruise_id = ?
        AND station_id = ?
        AND dbar_bin BETWEEN ? AND ?
        {soak_clause}
        {qc_clause} ORDER BY dbar_bin ASC
    """
    df = con.execute(query, (target_cruise, target_id, z_range[0], z_range[1])).df()
    if df.empty or 'theta' not in df.columns: return df

    # EOS-80 Analytics
    df['sigma'] = df['rho'] - 1000
    # theta is stored as IPTS-68 (used by EOS-80 equations); gsw functions
    # expect ITS-90 potential temperature, so convert here for derived quantities.
    theta_90 = df['theta'] / 1.00024
    df['sat_o2'] = gsw.O2sol_SP_pt(df['SP'], theta_90)
    df['AOU'] = df['sat_o2'] - df['o2_final']

    # Metabolic Index (Phi) — uses ITS-90 potential temperature in Arrhenius term
    k, Eo = 8.617e-5, 0.45
    df['phi'] = (df['o2_final'] / np.exp(-Eo / (k * (theta_90 + 273.15)))) / 1e6

    # CHL cannot be negative — truncate to 0 (Gaussian smoothing can push
    # near-zero values slightly below zero at the base of the DCM).
    if 'chl_final' in df.columns:
        df['chl_final'] = df['chl_final'].clip(lower=0.0)

    return df

def download_csv():
    df = get_clean_df(cruise_select.value, station_select.value, depth_slider.value, qc_checkbox.value, soak_toggle.value)
    sio = BytesIO()
    df.to_csv(sio, index=False); sio.seek(0)
    return sio

csv_button = pn.widgets.FileDownload(callback=download_csv, filename='CTD_Export.csv', label='Export CSV', button_type='primary', sizing_mode='stretch_width')

# 5. TAB FUNCTIONS
def _cruise_summary_impl(target_cruise):
    query = f"""
        SELECT station_id as "Station ID", wf_cast as "WF #", sb_cast as "SB #",
        MIN(lat) as "Lat", MIN(lon) as "Lon", MIN(time_iso)::TIMESTAMP as "Start Time",
        CAST(MAX(depth_m) AS DECIMAL(10,1)) as "Max Depth (m)",
        ROUND((SUM(CASE WHEN qc_flag = 1 THEN 1 ELSE 0 END) * 100.0) / COUNT(*), 1) as "Health %"
        FROM {TABLE_NAME} WHERE cruise_id = ?
        GROUP BY station_id, wf_cast, sb_cast ORDER BY "Start Time" ASC
    """
    df_sum = con.execute(query, (target_cruise,)).df()
    config = {'columns': [{'field': 'Health %', 'formatter': 'progress', 'formatterParams': {'color': '#28a745', 'legend': True}}]}
    table = pn.widgets.Tabulator(df_sum, theme='midnight', show_index=False, sizing_mode='stretch_both', configuration=config)
    return pn.Column("# Cruise Summary", table, sizing_mode='stretch_both')

view_cruise_summary = pn.Column(
    pn.bind(_cruise_summary_impl, cruise_select),
    sizing_mode='stretch_both'
)

@pn.depends(cruise_select, station_select, depth_slider, qc_checkbox, soak_toggle)
def view_profiles(target_cruise, target_id, z_range, filter_qc, show_soak):
    df = get_clean_df(target_cruise, target_id, z_range, filter_qc, show_soak)
    if df.empty: return pn.pane.Alert("Data Pending...")
    v_opts = dict(invert_yaxis=True, height=550, show_grid=True, xaxis='top', tools=['hover'], xticks=3, padding=0.05,
                  fontsize={'labels': '8pt', 'xticks': '7pt', 'yticks': '7pt', 'legend': '7pt'})

    p1 = hv.Curve(df, 'theta', 'depth_m', label='Pot. Temp').opts(**v_opts, color='blue', width=175, xlabel='theta (°C)')
    p6 = hv.Curve(df, 'in_situ_temp', 'depth_m', label='In-Situ Temp').opts(**v_opts, color='purple', width=140, yaxis=None, xlabel='T (°C)')
    p2 = hv.Curve(df, 'SP', 'depth_m', label='Prac. Sal').opts(**v_opts, color='red', width=140, yaxis=None, xlabel='SP (PSU)')
    p3 = hv.Curve(df, 'o2_final', 'depth_m', label='O2').opts(**v_opts, color='black', width=155, yaxis=None, xlabel='O2 (µmol/kg)')
    p4 = hv.Curve(df, 'ph_final', 'depth_m', label='pH').opts(**v_opts, color='orange', width=125, yaxis=None, xlabel='pH')
    p5 = hv.Curve(df, 'chl_final', 'depth_m', label='Chl').opts(**v_opts, color='green', width=125, yaxis=None, xlabel='Chl (mg/m³)')

    return (p1 + p6 + p2 + p3 + p4 + p5).cols(6).opts(shared_axes=True, merge_tools=True)

@pn.depends(cruise_select, station_select, depth_slider, qc_checkbox, soak_toggle)
def view_ts_analysis(target_cruise, target_id, z_range, filter_qc, show_soak):
    df = get_clean_df(target_cruise, target_id, z_range, filter_qc, show_soak)
    if df.empty: return pn.pane.Alert("Data Null")
    return hv.Points(df, ['SP', 'theta'], ['depth_m', 'sigma']).opts(
        color='depth_m', cmap='Viridis_r', width=600, height=500, colorbar=True, title="T-S Analysis (EOS-80)", tools=['hover']
    )

@pn.depends(cruise_select, station_select, depth_slider, qc_checkbox, soak_toggle)
def view_aou(target_cruise, target_id, z_range, filter_qc, show_soak):
    df = get_clean_df(target_cruise, target_id, z_range, filter_qc, show_soak)
    if df.empty: return pn.pane.Alert("No Data")
    opts = dict(invert_yaxis=True, height=550, width=600, show_grid=True, tools=['hover'])
    sat_l = hv.Curve(df, 'sat_o2', 'depth_m', label='Sat. Cap').opts(**opts, color='black', line_dash='dashed')
    o2_l = hv.Curve(df, 'o2_final', 'depth_m', label='Observed').opts(**opts, color='cyan')
    fill = hv.Area(df, ('sat_o2', 'o2_final'), 'depth_m', label='AOU').opts(**opts, color='orange', alpha=0.3)
    return (fill * sat_l * o2_l).opts(title="Apparent Oxygen Utilization (AOU)")

@pn.depends(cruise_select, station_select, depth_slider, qc_checkbox, soak_toggle)
def view_stability(target_cruise, target_id, z_range, filter_qc, show_soak):
    df = get_clean_df(target_cruise, target_id, z_range, filter_qc, show_soak)
    if df.empty: return pn.pane.Alert("No Data")
    surf = df['sigma'].iloc[0]
    opts = dict(invert_yaxis=True, height=500, width=400, tools=['hover'])
    density_curve = hv.Curve(df, 'sigma', 'depth_m', label='Density').opts(**opts, color='blue')
    # Guard against fully-mixed profiles: only draw MLD line if threshold is crossed
    over_threshold = (df['sigma'] - surf) > 0.03
    if over_threshold.any():
        mld_v = df.loc[over_threshold.idxmax(), 'depth_m']
        density_panel = density_curve * hv.HLine(mld_v).opts(color='red')
    else:
        density_panel = density_curve
    return pn.Row(density_panel,
                  hv.Area(df, 'o2_final', 'depth_m', label='Oxygen Concentration').opts(**opts, color='magenta', alpha=0.2))

@pn.depends(cruise_select, station_select, depth_slider, qc_checkbox, soak_toggle)
def view_metabolic_index(target_cruise, target_id, z_range, filter_qc, show_soak):
    df = get_clean_df(target_cruise, target_id, z_range, filter_qc, show_soak)
    if df.empty: return pn.pane.Alert("No Data")
    return (hv.Curve(df, 'phi', 'depth_m', label='Φ').opts(color='#e67e22', invert_yaxis=True, height=550, width=500, tools=['hover']) * hv.VLine(1.0).opts(color='red', line_dash='dashed'))

@pn.depends(cruise_select, station_select)
def view_map_geolocation(target_cruise, target_id):
    cruise_coords = con.execute(f"SELECT DISTINCT station_id, lat, lon FROM {TABLE_NAME} WHERE cruise_id = ?", (target_cruise,)).df()
    pts = gv.Points(cruise_coords, ['lon', 'lat'], vdims=['station_id'], crs=ccrs.PlateCarree()).opts(size=8, color='#f1c40f', alpha=0.6, tools=['hover'])
    sel = gv.Points(cruise_coords[cruise_coords['station_id'] == target_id], ['lon', 'lat'], crs=ccrs.PlateCarree()).opts(size=18, color='red', marker='circle', line_color='white')
    return (gvts.EsriOceanBase * gvts.EsriOceanReference * pts * sel).opts(width=900, height=600, title="Geolocation")

def _tabular_impl(target_cruise, target_id, z_range, filter_qc, show_soak):
    df = get_clean_df(target_cruise, target_id, z_range, filter_qc, show_soak)
    if df.empty:
        return pn.pane.Alert("No data for the current selection.")
    return pn.widgets.Tabulator(df, pagination='remote', page_size=15, theme='midnight', sizing_mode='stretch_both')

view_tabular_data = pn.Column(
    pn.bind(_tabular_impl, cruise_select, station_select, depth_slider, qc_checkbox, soak_toggle),
    sizing_mode='stretch_both'
)

# 6. ASSEMBLY
tabs = pn.Tabs(
    ("Vertical Profiles", view_profiles),
    ("T-S Analysis", view_ts_analysis),
    ("Oxygen Utilization (AOU)", view_aou),
    ("Stability & MLD", view_stability),
    ("Metabolic Index", view_metabolic_index),
    ("Geolocation", view_map_geolocation),
    ("Tabular Data", view_tabular_data),
    ("Cruise Summary", view_cruise_summary),
    dynamic=True, active=0
)

dashboard = pn.template.FastListTemplate(
    title="Western Flyer - CTD Data",
    sidebar=[cruise_select, station_select, depth_slider, qc_checkbox, soak_toggle, pn.pane.Markdown("---"), csv_button],
    main=[tabs], accent_base_color="#00f2ff", header_background="#1a1a1a"
)

dashboard.servable()
