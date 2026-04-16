import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.patches as mpatches

# ---------------------------------------------------------
# CONFIG
# ---------------------------------------------------------
egms_layer_name = "egms_sthlm_coh60_SW99TM_cor50_statsLC"
psds_layer_name = "PSDS_sthlm_coh60_SW99TM_cor50_statsLC"

CLC = "CLC_sthlm_cor50"
LC_FIELD = "Code_18"
AREA_FIELD = "area_m2"

# PSDS fields
PSDS_RMSE = "rmse_median"
PSDS_COH = "coherence_median"
PSDS_COUNT = "coherence_count"

# EGMS fields
EGMS_RMSE = "rmse_median"
EGMS_COH = "temporal_coherence_median"
EGMS_COUNT = "rmse_count"

# ---------------------------------------------------------
# LOAD LAYERS
# ---------------------------------------------------------
def load_stats(layer_name):
    layer = QgsProject.instance().mapLayersByName(layer_name)[0]
    records = []
    fields = [f.name() for f in layer.fields()]
    for feat in layer.getFeatures():
        rec = dict(zip(fields, feat.attributes()))
        records.append(rec)
    return pd.DataFrame(records)

psds = load_stats(psds_layer_name)
egms = load_stats(egms_layer_name)

# ---------------------------------------------------------
# CLEAN NUMERICS
# ---------------------------------------------------------
numeric_fields = [
    PSDS_RMSE, PSDS_COH, PSDS_COUNT, AREA_FIELD,
    EGMS_RMSE, EGMS_COH, EGMS_COUNT
]

for f in numeric_fields:
    if f in psds.columns:
        psds[f] = pd.to_numeric(psds[f], errors="coerce")
    if f in egms.columns:
        egms[f] = pd.to_numeric(egms[f], errors="coerce")

# Remove water course class and forest

#remove_classes = ["511", "312"]

#psds = psds[~psds[LC_FIELD].isin(remove_classes)]
#egms = egms[~egms[LC_FIELD].isin(remove_classes)]

# ---------------------------------------------------------
# COMPUTE DENSITY (sum count / sum area)
# ---------------------------------------------------------
density_ps = (
    psds.groupby(LC_FIELD)[PSDS_COUNT].sum() /
    psds.groupby(LC_FIELD)[AREA_FIELD].sum()
) * 1e6

density_eg = (
    egms.groupby(LC_FIELD)[EGMS_COUNT].sum() /
    egms.groupby(LC_FIELD)[AREA_FIELD].sum()
) * 1e6

# Order LC classes by PSDS density
lc_order = density_ps.sort_values(ascending=False).index.tolist()

# ---------------------------------------------------------
# PREPARE BOX PLOT DATA
# ---------------------------------------------------------
box_rmse_ps = [psds[psds[LC_FIELD] == lc][PSDS_RMSE].dropna() for lc in lc_order]
box_rmse_eg = [egms[egms[LC_FIELD] == lc][EGMS_RMSE].dropna() for lc in lc_order]

box_coh_ps  = [psds[psds[LC_FIELD] == lc][PSDS_COH].dropna() for lc in lc_order]
box_coh_eg  = [egms[egms[LC_FIELD] == lc][EGMS_COH].dropna() for lc in lc_order]

# ---------------------------------------------------------
# Get land cover codes' names'
# ---------------------------------------------------------
def get_lc_names_from_layer(layer_name, lc_field):
    """
    Reads category labels from QGIS symbology.
    Returns dict: {code: label}.
    """
    layer = QgsProject.instance().mapLayersByName(layer_name)[0]
    renderer = layer.renderer()

    lc_map = {}

    # Works for categorized symbology
    if hasattr(renderer, "categories"):
        for cat in renderer.categories():
            code = cat.value()
            label = cat.label()
            lc_map[str(code)] = label
            lc_map[int(code)] = label  # allow int lookup too

    return lc_map
LC_NAMES = get_lc_names_from_layer(CLC, LC_FIELD)
lc_labels = [
    LC_NAMES.get(code) or 
    LC_NAMES.get(str(code)) or 
    LC_NAMES.get(int(code)) or 
    str(code)
    for code in lc_order
]
# ---------------------------------------------------------
# PLOTTING
# ---------------------------------------------------------
fig, axes = plt.subplots(3, 1, figsize=(18, 20))
ax_rmse, ax_coh, ax_den = axes
ps_patch = mpatches.Patch(facecolor="tab:blue",  alpha=0.5, label="PSDS")
eg_patch = mpatches.Patch(facecolor="tab:orange", alpha=0.5, label="EGMS")
# ======================================
# RMSE BOXPLOTS
# ======================================
positions_ps = np.arange(len(lc_order)) - 0.2
positions_eg = np.arange(len(lc_order)) + 0.2

ax_rmse.boxplot(
    box_rmse_ps, positions=positions_ps, widths=0.35,
    patch_artist=True, boxprops=dict(facecolor="tab:blue", alpha=0.5),
    medianprops=dict(color="black")
)

ax_rmse.boxplot(
    box_rmse_eg, positions=positions_eg, widths=0.35,
    patch_artist=True, boxprops=dict(facecolor="tab:orange", alpha=0.5),
    medianprops=dict(color="black")
)

ax_rmse.set_title("RMSE Distribution per Landcover")
ax_rmse.set_ylabel("RMSE (mm)")
ax_rmse.set_xticks(np.arange(len(lc_order)))
ax_rmse.set_xticklabels(lc_order, rotation=45, ha="right")
ax_rmse.grid(axis='y', alpha=0.3)
ax_rmse.legend(handles=[ps_patch, eg_patch], loc="upper right")

# ======================================
# COHERENCE BOXPLOTS
# ======================================
ax_coh.boxplot(
    box_coh_ps, positions=positions_ps, widths=0.35,
    patch_artist=True, boxprops=dict(facecolor="tab:blue", alpha=0.5),
    medianprops=dict(color="black")
)

ax_coh.boxplot(
    box_coh_eg, positions=positions_eg, widths=0.35,
    patch_artist=True, boxprops=dict(facecolor="tab:orange", alpha=0.5),
    medianprops=dict(color="black")
)

ax_coh.set_title("Coherence Distribution per Landcover")
ax_coh.set_ylabel("Coherence")
ax_coh.set_ylim(0, 1.0)
ax_coh.set_xticks(np.arange(len(lc_order)))
ax_coh.set_xticklabels(lc_order, rotation=45, ha="right")
ax_coh.grid(axis='y', alpha=0.3)
ax_coh.legend(handles=[ps_patch, eg_patch], loc="upper right")

# ======================================
# DENSITY BAR PLOT (KEEP AS BAR)
# ======================================
x = np.arange(len(lc_order))
bar_width = 0.35

ax_den.bar(x - bar_width/2, density_ps.loc[lc_order], bar_width, label="PSDS", color="tab:blue")
ax_den.bar(x + bar_width/2, density_eg.loc[lc_order], bar_width, label="EGMS", color="tab:orange")

ax_den.set_title("Point Density per Landcover (log scale)")
ax_den.set_ylabel("Points per km²")
ax_den.set_yscale("log")
ax_den.set_xticks(x)
ax_den.set_xticklabels(lc_order, rotation=45, ha="right")
ax_den.grid(axis='y', alpha=0.3)
ax_den.legend()

# ---------------------------------------------------------
# LANDCOVER LEGEND (code → name)
# ---------------------------------------------------------
legend_text = "\n".join([
    f"{code}: {LC_NAMES.get(str(code), 'Unknown')}"
    for code in lc_order
])

# Add side text box
fig.text(
    0.2, 0.45, legend_text,
    va='center', ha='left',
    fontsize=9, family='monospace'
)

plt.subplots_adjust(right=0.78)   # ← creates space on the right for the legend
plt.tight_layout(rect=[0,0,0.78,1])  
plt.show()