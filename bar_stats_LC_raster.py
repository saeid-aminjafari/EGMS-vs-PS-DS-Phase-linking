import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.patches as mpatches
from osgeo import gdal
# Set global default font size for tick labels and axis labels
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['axes.labelsize'] = 14  # Optional: also make axis titles bigger

# =========================================================
# CONFIG
# =========================================================
raster_name = "NMD2023bas_v0_2_LO"   # paletted raster
psds_layer_name = "PSDS_LOA_coh70_SW99TM_statsLC_NMD"
egms_layer_name = "egms_LOA_coh70_SW99TM_statsLC_NMD"

LC_FIELD = "LC_nmd2023v021"   # sampled raster LC code
RMSE_FIELD = "rmse"           # both layers use rmse

# =========================================================
# LOAD POINT LAYERS
# =========================================================
def load_points(layer_name):
    layer = QgsProject.instance().mapLayersByName(layer_name)[0]
    records = []
    fields = [f.name() for f in layer.fields()]
    for feat in layer.getFeatures():
        rec = dict(zip(fields, feat.attributes()))
        records.append(rec)
    return pd.DataFrame(records)

ps = load_points(psds_layer_name)
eg = load_points(egms_layer_name)

# force numeric
ps[LC_FIELD] = pd.to_numeric(ps[LC_FIELD], errors="coerce")
eg[LC_FIELD] = pd.to_numeric(eg[LC_FIELD], errors="coerce")
ps[RMSE_FIELD] = pd.to_numeric(ps[RMSE_FIELD], errors="coerce")
eg[RMSE_FIELD] = pd.to_numeric(eg[RMSE_FIELD], errors="coerce")

# =========================================================
# GET LC CLASS NAMES FROM RASTER SYMBOLOGY
# =========================================================
def get_lc_names_from_raster(raster_name):
    raster = QgsProject.instance().mapLayersByName(raster_name)[0]
    renderer = raster.renderer()

    lc_map = {}

    # Paletted symbology
    if hasattr(renderer, "classes"):
        for c in renderer.classes():
            # c.value is a float, c.label is already a string
            try:
                code = int(c.value)
                label = c.label   # not c.label()
                lc_map[code] = label
            except:
                pass

    # Additional fallback: legend items
    try:
        for (sym, lbl) in renderer.legendSymbologyItems():
            parts = lbl.split(" ")
            if parts[0].isdigit():
                lc_map[int(parts[0])] = lbl
    except:
        pass

    return lc_map

LC_NAMES = get_lc_names_from_raster(raster_name)
# =========================================================
#COMPUTE AREA PER LAND COVER CLASS FROM RASTER
# =========================================================
def compute_raster_area(layer_name):
    layer = QgsProject.instance().mapLayersByName(layer_name)[0]
    path = layer.dataProvider().dataSourceUri().split("|")[0]

    ds = gdal.Open(path)
    band = ds.GetRasterBand(1)

    # Pixel size from geotransform
    gt = ds.GetGeoTransform()
    px = abs(gt[1])
    py = abs(gt[5])
    pixel_area = px * py  # m²

    # Read entire raster safely
    arr = band.ReadAsArray()

    area_map = {}
    unique, counts = np.unique(arr, return_counts=True)

    for val, cnt in zip(unique, counts):
        code = int(val)
        area_map[code] = cnt * pixel_area

    return area_map

area_lc = compute_raster_area(raster_name)
# Remove LC = 0 (unknown/outside corridor)
if 0 in area_lc:
    del area_lc[0]

ps = ps[ps[LC_FIELD] != 0]
eg = eg[eg[LC_FIELD] != 0]
# =========================================================
# DETERMINE LANDCOVER CLASSES TO USE
# =========================================================
classes_with_points = sorted(set(ps[LC_FIELD].dropna().unique()) |
                             set(eg[LC_FIELD].dropna().unique()))
# =========================================================
# COUNT POINTS PER LC (percentage)
# =========================================================
ps_counts = ps.groupby(LC_FIELD)[RMSE_FIELD].count()
eg_counts = eg.groupby(LC_FIELD)[RMSE_FIELD].count()

total_ps = ps_counts.sum()
total_eg = eg_counts.sum()

# sort by area (descending)
classes_with_points = sorted(classes_with_points, key=lambda lc: area_lc.get(int(lc), 0), reverse=True)

percent_ps = [(ps_counts.get(lc, 0) / total_ps * 100) for lc in classes_with_points]
percent_eg = [(eg_counts.get(lc, 0) / total_eg * 100) for lc in classes_with_points]

# area in km²
area_km2 = [area_lc.get(int(lc), 0) / 1e6 for lc in classes_with_points]
# =========================================================
# BOXPLOT DATA
# =========================================================
box_ps = [ps[ps[LC_FIELD] == lc][RMSE_FIELD].dropna() for lc in classes_with_points]
box_eg = [eg[eg[LC_FIELD] == lc][RMSE_FIELD].dropna() for lc in classes_with_points]

# =========================================================
# PLOTS
# =========================================================
fig, (ax_box, ax_bar) = plt.subplots(1, 2, figsize=(12, 3))
fig.patch.set_facecolor('none')
# ---------------------------------------------------------
# RMSE BOXPLOTS (PSDS vs EGMS)
# ---------------------------------------------------------
positions_ps = np.arange(len(classes_with_points)) - 0.2
positions_eg = np.arange(len(classes_with_points)) + 0.2

ax_box.boxplot(box_ps, positions=positions_ps, widths=0.35,
               patch_artist=True, boxprops=dict(facecolor="orange", alpha=1),
               medianprops=dict(color="black"))

ax_box.boxplot(box_eg, positions=positions_eg, widths=0.35,
               patch_artist=True, boxprops=dict(facecolor="blue", alpha=1),
               medianprops=dict(color="black"))
ax_box.set_facecolor('none') # Remove panel background
ax_box.set_title("e) Lodose (RMSE)")
ax_box.set_ylabel("RMSE (mm)")
ax_box.set_xticks(np.arange(len(classes_with_points)))
ax_box.set_xticklabels(classes_with_points, rotation=45, ha="right")
ax_box.grid(axis='y', alpha=0.3)
psds_patch = mpatches.Patch(facecolor="orange", alpha=1, label="PSDS")
egms_patch = mpatches.Patch(facecolor="blue", alpha=1, label="EGMS")
ax_box.legend(handles=[psds_patch, egms_patch])


# ---------------------------------------------------------
# BARPLOT: PERCENT SHARE OF LC CLASSES
# ---------------------------------------------------------
x = np.arange(len(classes_with_points))
bw = 0.35

bars_ps = ax_bar.bar(x - bw/2, percent_ps, bw, label="PSDS", color="orange", alpha=1)
bars_eg = ax_bar.bar(x + bw/2, percent_eg, bw, label="EGMS", color="blue", alpha=1)
ax_bar.set_facecolor('none') # Remove panel background
ax_bar.set_title("f) Lodose (Point share)")
ax_bar.set_ylabel("% of all points")
ax_bar.set_xticks(x)
ax_bar.set_xticklabels(classes_with_points, rotation=45, ha="right")
ax_bar.grid(axis='y', alpha=0.3)
ax_bar.legend(["PSDS", "EGMS"])
# --- Annotate bars with point counts ---
for rect, lc in zip(bars_ps, classes_with_points):
    height = rect.get_height()
    N = ps_counts.get(lc, 0)
    ax_bar.text(rect.get_x() + rect.get_width()/2, height + 0.5, f"{N}",
                ha='center', va='bottom', fontsize=8)
for rect, lc in zip(bars_eg, classes_with_points):
    height = rect.get_height()
    N = eg_counts.get(lc, 0)
    ax_bar.text(rect.get_x() + rect.get_width()/2, height + 0.5, f"{N}",
                ha='center', va='bottom', fontsize=8)
# ---------------------------------------------------------
# OVERLAY: LAND-COVER AREA AS A SECONDARY AXIS
# ---------------------------------------------------------
ax_area = ax_bar.twinx()
ax_area.plot(x, area_km2, color="black", marker="o", linestyle="-",
             linewidth=1.2, label="area (km²)")

ax_area.set_ylabel("Area (km²)")

ax_area.grid(False)
# ---------------------------------------------------------
# TEXT LEGEND (LC code → LC class name)
# ---------------------------------------------------------
#legend_text = "\n".join([
#    f"{lc}: {LC_NAMES.get(int(lc), 'Unknown')}"
#    for lc in classes_with_points
#])

#fig.text(
#    0.4, 0.3, legend_text,
#    va='center', ha='left',
#    fontsize=9, family='monospace'
#)

#plt.subplots_adjust(right=0.78)
plt.tight_layout() #rect=[0,0,0.78,1]
plt.show()