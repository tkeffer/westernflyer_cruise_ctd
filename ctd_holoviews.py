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
from scipy.interpolate import griddata
from holoviews.operation import contours as hv_contours
from bokeh.models import HoverTool, ColumnDataSource, LinearAxis, Range1d
from bokeh.plotting import figure

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
                  line_width=2.5, fontsize={'labels': '8pt', 'xticks': '7pt', 'yticks': '7pt', 'legend': '7pt'})

    p1 = hv.Curve(df, 'theta', 'depth_m', label='Pot. Temp').opts(**v_opts, color='cyan', width=175, xlabel='theta (°C)')
    p6 = hv.Curve(df, 'in_situ_temp', 'depth_m', label='In-Situ Temp').opts(**v_opts, color='blue', width=140, yaxis=None, xlabel='T (°C)')
    p2 = hv.Curve(df, 'SP', 'depth_m', label='Prac. Sal').opts(**v_opts, color='red', width=140, yaxis=None, xlabel='SP (PSU)')
    p3 = hv.Curve(df, 'o2_final', 'depth_m', label='O2').opts(**v_opts, color='black', width=155, yaxis=None, xlabel='O2 (µmol/kg)')
    p4 = hv.Curve(df, 'ph_final', 'depth_m', label='pH').opts(**v_opts, color='goldenrod', width=125, yaxis=None, xlabel='pH')
    p5 = hv.Curve(df, 'chl_final', 'depth_m', label='Chl').opts(**v_opts, color='green', width=125, yaxis=None, xlabel='Chl (mg/m³)')

    return (p1 + p6 + p2 + p3 + p4 + p5).cols(6).opts(shared_axes=True, merge_tools=True)

@pn.depends(cruise_select, station_select, depth_slider, qc_checkbox, soak_toggle)
def view_multi_profile(target_cruise, target_id, z_range, filter_qc, show_soak):
    df = get_clean_df(target_cruise, target_id, z_range, filter_qc, show_soak)
    if df.empty: return pn.pane.Alert("Data Pending...")

    source = ColumnDataSource(df)

    def get_range(col):
        m, M = df[col].min(), df[col].max()
        if pd.isna(m) or pd.isna(M):
            return Range1d(0, 1)
        if m == M:
            return Range1d(m - 1, m + 1)
        padding = (M - m) * 0.05
        return Range1d(m - padding, M + padding)

    # Base figure
    p = figure(height=700, width=800,
               y_axis_label='Depth (m)',
               y_range=Range1d(df['depth_m'].max(), df['depth_m'].min()),
               toolbar_location="right",
               tools="pan,wheel_zoom,box_zoom,reset,save")

    p.xaxis.visible = False # Hide default xaxis

    # Define extra ranges
    p.extra_x_ranges = {
        "ph": get_range('ph_final'),
        "o2": get_range('o2_final'),
        "sal": get_range('SP'),
        "temp": get_range('in_situ_temp'),
        "chl": get_range('chl_final'),
    }

    # Add lines
    p.line('ph_final', 'depth_m', source=source, x_range_name="ph", color='goldenrod', line_width=2.5, legend_label='pH')
    p.line('o2_final', 'depth_m', source=source, x_range_name="o2", color='black', line_width=2.5, legend_label='Oxygen')
    p.line('SP', 'depth_m', source=source, x_range_name="sal", color='red', line_width=2.5, legend_label='Salinity')
    p.line('in_situ_temp', 'depth_m', source=source, x_range_name="temp", color='blue', line_width=2.5, legend_label='In-Situ Temp')
    p.line('chl_final', 'depth_m', source=source, x_range_name="chl", color='green', line_width=2.5, legend_label='Chlorophyll')

    # Top axes (from top to bottom: pH, Oxygen, Salinity)
    # Bokeh stacks axes outward from the plot.
    # To get pH (top), Oxygen (middle), Salinity (bottom/closest to plot) on top of the plot:
    # 1. Salinity (added first, closest)
    # 2. Oxygen
    # 3. pH (added last, furthest)
    p.add_layout(LinearAxis(x_range_name="sal", axis_label="Salinity (PSU)",
                            axis_label_text_color="red", axis_line_color="red",
                            major_label_text_color="red", major_tick_line_color="red"), 'above')
    p.add_layout(LinearAxis(x_range_name="o2", axis_label="Oxygen (µmol/kg)",
                            axis_label_text_color="black", axis_line_color="black",
                            major_label_text_color="black", major_tick_line_color="black"), 'above')
    p.add_layout(LinearAxis(x_range_name="ph", axis_label="pH",
                            axis_label_text_color="goldenrod", axis_line_color="goldenrod",
                            major_label_text_color="goldenrod", major_tick_line_color="goldenrod"), 'above')

    # Bottom axes (from top to bottom: Temperature, Chlorophyll)
    # "top to bottom" for bottom axes means Temperature is closer to the plot than Chlorophyll.
    # 1. Temperature (added first, closest)
    # 2. Chlorophyll (added last, furthest)
    p.add_layout(LinearAxis(x_range_name="temp", axis_label="Temperature (°C)",
                            axis_label_text_color="blue", axis_line_color="blue",
                            major_label_text_color="blue", major_tick_line_color="blue"), 'below')
    p.add_layout(LinearAxis(x_range_name="chl", axis_label="Chlorophyll (mg/m³)",
                            axis_label_text_color="green", axis_line_color="green",
                            major_label_text_color="green", major_tick_line_color="green"), 'below')

    hover = HoverTool(tooltips=[
        ("Depth", "@depth_m{0.1} m"),
        ("pH", "@ph_final{0.00}"),
        ("Oxygen", "@o2_final{0.1}"),
        ("Salinity", "@SP{0.00}"),
        ("Temp", "@in_situ_temp{0.2}"),
        ("Chl", "@chl_final{0.2}"),
    ], mode='vline')
    p.add_tools(hover)

    p.legend.location = "top_left"
    p.legend.click_policy = "hide"
    p.legend.label_text_font_size = "8pt"
    p.legend.background_fill_alpha = 0.5

    return p

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

# ── Section plot helpers ───────────────────────────────────────────────────────
def _haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance between two coordinate pairs in kilometres."""
    R = 6371.0
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = (np.sin(dlat / 2) ** 2
         + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon / 2) ** 2)
    return 2 * R * np.arcsin(np.sqrt(a))

_SECTION_VARS = {
    'Pot. Temperature (°C)': ('theta',        'RdBu_r',  '°C'),
    'In-Situ Temp (°C)':     ('in_situ_temp', 'RdBu_r',  '°C'),
    'Practical Salinity':    ('SP',           'viridis', 'PSU'),
    'Density σ (kg/m³)':     ('sigma',        'viridis', 'kg/m³'),
    'Oxygen (µmol/kg)':      ('o2_final',     'plasma',  'µmol/kg'),
    'pH':                    ('ph_final',     'RdYlBu',  ''),
    'Chlorophyll (mg/m³)':   ('chl_final',    'Greens',  'mg/m³'),
}

section_var_select = pn.widgets.Select(
    name='Section Variable',
    options=list(_SECTION_VARS.keys()),
    value='Pot. Temperature (°C)',
    width=230,
)

def _section_impl(target_cruise, z_range, filter_qc, section_var_label):
    col, cmap, unit = _SECTION_VARS[section_var_label]
    qc_clause = "AND qc_flag < 3" if filter_qc else ""
    query = f"""
        SELECT station_id, station_name, lat, lon, time_iso, depth_m, dbar_bin,
               theta, in_situ_temp, SP, rho, o2_final, ph_final, chl_final
        FROM {TABLE_NAME}
        WHERE cruise_id = ?
        AND dbar_bin BETWEEN ? AND ?
        AND is_soak = 0
        {qc_clause}
        ORDER BY station_id, dbar_bin ASC
    """
    df = con.execute(query, (target_cruise, z_range[0], z_range[1])).df()
    if df.empty:
        return pn.pane.Alert("No data available for section plot.")

    # Always compute sigma — needed for isopycnal contours regardless of variable
    df['sigma'] = df['rho'] - 1000

    if col not in df.columns:
        return pn.pane.Alert(f"Column '{col}' not found in data.")

    # Order stations chronologically along the cruise track
    station_order = (
        df.groupby('station_id')['time_iso'].min()
        .sort_values().index.tolist()
    )

    # Cumulative along-track distance (km) for each station
    sta_pos = df.groupby('station_id')[['lat', 'lon']].first().loc[station_order]
    lats, lons = sta_pos['lat'].values, sta_pos['lon'].values
    cum_dist = np.zeros(len(station_order))
    for i in range(1, len(station_order)):
        cum_dist[i] = cum_dist[i - 1] + _haversine_km(
            lats[i - 1], lons[i - 1], lats[i], lons[i]
        )
    dist_map = dict(zip(station_order, cum_dist))
    df['dist_km'] = df['station_id'].map(dist_map)

    # Base (x, y) grid — shared by main variable and sigma
    x_all = df['dist_km'].values.astype(float)
    y_all = df['depth_m'].values.astype(float)
    base_mask = ~(np.isnan(x_all) | np.isnan(y_all))

    # ── Main variable ─────────────────────────────────────────────────────────
    z = df[col].values.astype(float)
    mask = base_mask & ~np.isnan(z)
    xm, ym, zm = x_all[mask], y_all[mask], z[mask]

    if len(xm) < 10:
        return pn.pane.Alert("Insufficient data for section interpolation.")

    # Regular 300 × 200 grid
    xi = np.linspace(xm.min(), xm.max(), 300)
    yi = np.linspace(ym.min(), ym.max(), 200)
    Xi, Yi = np.meshgrid(xi, yi)
    Zi = griddata((xm, ym), zm, (Xi, Yi), method='linear')

    vmin = np.nanpercentile(zm, 2)
    vmax = np.nanpercentile(zm, 98)
    label = section_var_label + (f' [{unit}]' if unit else '')

    img = hv.Image(
        (xi, yi, Zi),
        kdims=['Distance Along Track (km)', 'Depth (m)'],
        vdims=[label],
    ).opts(
        cmap=cmap, colorbar=True, clim=(vmin, vmax),
        width=900, height=500, invert_yaxis=True,
        tools=['hover'],
        title=f"{section_var_label} Section — {target_cruise}",
        fontsize={'title': '10pt', 'labels': '9pt', 'xticks': '8pt', 'yticks': '8pt'},
    )

    # ── Isopycnal contours (σθ, every 0.5 kg/m³) ─────────────────────────────
    s_vals = df['sigma'].values.astype(float)
    s_mask = base_mask & ~np.isnan(s_vals)
    Zi_sigma = griddata(
        (x_all[s_mask], y_all[s_mask]), s_vals[s_mask], (Xi, Yi), method='linear'
    )
    s_clean = Zi_sigma[~np.isnan(Zi_sigma)]
    if len(s_clean) > 0:
        lvl_min = np.ceil(s_clean.min() * 2) / 2
        lvl_max = np.floor(s_clean.max() * 2) / 2
        levels = np.arange(lvl_min, lvl_max + 0.5, 0.5).tolist()
    else:
        levels = 10

    sigma_img = hv.Image(
        (xi, yi, Zi_sigma),
        kdims=['Distance Along Track (km)', 'Depth (m)'],
        vdims=['sigma'],
    )
    iso = hv_contours(sigma_img, levels=levels).opts(
        line_color='black', line_width=0.9, show_legend=False
    )

    # ── Station marker lines ──────────────────────────────────────────────────
    sta_lines = hv.Overlay([
        hv.VLine(d).opts(color='white', line_width=0.8, line_dash='dashed', alpha=0.5)
        for d in cum_dist
    ])

    # ── Station hover markers (triangle at surface, shows name on hover) ────────
    sta_names = df.groupby('station_id')['station_name'].first()
    marker_depth = yi.min() + (yi.max() - yi.min()) * 0.02
    sta_hover = hv.Points(
        pd.DataFrame({
            'Distance Along Track (km)': cum_dist,
            'Depth (m)':                 np.full(len(station_order), marker_depth),
            'Station':                   [sta_names.get(s, s) for s in station_order],
        }),
        kdims=['Distance Along Track (km)', 'Depth (m)'],
        vdims=['Station'],
    ).opts(
        color='white', size=7, marker='triangle', show_legend=False,
        tools=[HoverTool(tooltips=[
            ('Station',                   '@{Station}'),
            ('Distance Along Track (km)', '@{Distance Along Track (km)}'),
            ('Depth (m)',                 '@{Depth (m)}'),
        ])],
    )

    return (img * iso * sta_lines * sta_hover).opts(show_grid=False)

view_section = pn.Column(
    pn.Row(section_var_select, margin=(8, 0, 4, 10)),
    pn.bind(_section_impl, cruise_select, depth_slider, qc_checkbox, section_var_select),
    sizing_mode='stretch_both',
)

# ── Tabular data ───────────────────────────────────────────────────────────────
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
    ("Cruise Summary", view_cruise_summary),
    ("Geolocation", view_map_geolocation),
    ("Vertical Section", view_section),
    ("Vertical Profiles", view_profiles),
    ("Multi-Axis Profile", view_multi_profile),
    ("T-S Analysis", view_ts_analysis),
    ("Stability & MLD", view_stability),
    ("Oxygen Utilization (AOU)", view_aou),
    ("Metabolic Index", view_metabolic_index),
    ("Tabular Data", view_tabular_data),
    dynamic=True, active=3
)

dashboard = pn.template.FastListTemplate(
    title="Western Flyer - CTD Data",
    sidebar=[cruise_select, station_select, depth_slider, qc_checkbox, soak_toggle, pn.pane.Markdown("---"), csv_button],
    main=[tabs], accent_base_color="#00f2ff", header_background="#1a1a1a"
)

dashboard.servable()
