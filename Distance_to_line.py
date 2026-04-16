import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
import sys
# Set global default font size for tick labels and axis labels
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['axes.labelsize'] = 14  # Optional: also make axis titles bigger
# ---------------------------------------------------------
# CONFIG – CHANGE NAMES / FIELDS IF NEEDED
# ---------------------------------------------------------
DIST_FIELD = "distance"

# Coherence fields
#COH_PSDS_FIELD = "coherence"
#COH_EGMS_FIELD = "temporal_coherence"
RMSE_PSDS_FIELD = "rmse"
RMSE_EGMS_FIELD = "rmse"

# DEM fields
DEM_PSDS_FIELD = "dem"
DEM_EGMS_FIELD = "height_wgs84"

layers_cfg = [
    ("PSDS_LOA_coh70_SW99TM_jl_dist", "PSDS", "rail"),
    ("PSDS_LOA_coh70_SW99TM_vl_dist", "PSDS", "road"),
    ("egms_LOA_coh70_SW99TM_jl_dist", "EGMS", "rail"),
    ("egms_LOA_coh70_SW99TM_vl_dist", "EGMS", "road"),
]

# Distance bands
bands = [
    ("0–5 m",    0.0, 5.0),
    ("5–10 m",   5.0, 10.0),
    ("10–15 m",  10.0, 15.0),
    ("15–20 m",  15.0, 20.0),
]

# ---------------------------------------------------------
# ROBUST MEDIAN SMOOTHING FUNCTION
# ---------------------------------------------------------
def smooth_median(x, y, window=0.20):
    x = np.asarray(x)
    y = np.asarray(y)

    # Guard: if no points -> return empty
    if x.size == 0:
        return x, y

    order = np.argsort(x)
    xs = x[order]
    ys = y[order]

    n = len(xs)
    half = int(max(3, n * window / 2))
    y_smooth = np.zeros(n)

    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half)
        y_smooth[i] = np.median(ys[lo:hi])

    return xs, y_smooth

# ---------------------------------------------------------
# LOAD SHORTEST-LINE LAYER WITH RMSE & DEM VALUES
# ---------------------------------------------------------
def load_shortest_layer(layer_name, product, infra):
    layer = QgsProject.instance().mapLayersByName(layer_name)[0]

    records = []
    fields = [f.name() for f in layer.fields()]

    for feat in layer.getFeatures():
        attrs = dict(zip(fields, feat.attributes()))
        d = attrs.get(DIST_FIELD)

        if product == "PSDS":
            rmse = attrs.get(RMSE_PSDS_FIELD)
            dem = attrs.get(DEM_PSDS_FIELD)
        else:
            rmse = attrs.get(RMSE_EGMS_FIELD)
            dem = attrs.get(DEM_EGMS_FIELD)

        records.append({
            "product": product,
            "infra": infra,
            "dist": d,
            "rmse": rmse,
            "dem": dem
        })

    df = pd.DataFrame(records)

    df["dist"] = pd.to_numeric(df["dist"], errors="coerce")
    df["rmse"]  = pd.to_numeric(df["rmse"],  errors="coerce")
    df["dem"]  = pd.to_numeric(df["dem"],  errors="coerce")

    return df.dropna(subset=["dist"])

# ---------------------------------------------------------
# LOAD ALL
# ---------------------------------------------------------
dfs = []
for layer_name, prod, infra in layers_cfg:
    dfs.append(load_shortest_layer(layer_name, prod, infra))

df = pd.concat(dfs, ignore_index=True)

# Only distances 0–20 m
df_corr = df[df["dist"] <= 20.0].copy()

# ---------------------------------------------------------
# COUNTS / FRACTIONS PER BAND
# ---------------------------------------------------------
rows = []
for (prod, infra), sub in df_corr.groupby(["product", "infra"]):
    total = len(sub)
    if total == 0:
        continue

    for label, lo, hi in bands:
        count = ((sub["dist"] >= lo) & (sub["dist"] < hi)).sum()
        frac = count / total * 100.0
        rows.append({
            "product": prod,
            "infra": infra,
            "band": label,
            "count": count,
            "fraction": frac
        })

stats_df = pd.DataFrame(rows)

print("=== RAW COUNTS PER BAND ===")
print(stats_df.pivot_table(index=["product","infra"], columns="band",
                           values="count", fill_value=0))

print("\n=== FRACTION (%) PER BAND ===")
print(stats_df.pivot_table(index=["product","infra"], columns="band",
                           values="fraction", fill_value=0).round(2))

# =========================================================
# COMBINED PLOT: LEFT = STACKED BAR, RIGHT = BOXPLOTS
# =========================================================
fig_comb, (ax1, ax_box) = plt.subplots(
    1, 2, figsize=(12, 3.2), layout="constrained")
fig_comb.patch.set_facecolor('none') # No figure background
# ---------------------------------------------------------
# LEFT PANEL (ax1)– STACKED BAR PLOT
# ---------------------------------------------------------
ax1.set_facecolor('none') # No panel background
combo_order = [("PSDS", "rail"), ("EGMS", "rail"),
               ("PSDS", "road"), ("EGMS", "road")]
combo_labels = ["Rail-PSDS", "Rail-EGMS", "Road-PSDS", "Road-EGMS"]
plot_labels = combo_labels
frac_matrix = []

# compute total counts per combo to annotate later
total_counts = {
    label: len(df_corr[(df_corr["product"] == p) & (df_corr["infra"] == i)])
    for (p, i), label in zip(combo_order, combo_labels)
}

for band_label, _, _ in bands:
    row_vals = []
    for prod, infra in combo_order:
        sub = stats_df[(stats_df["product"] == prod) &
                       (stats_df["infra"] == infra) &
                       (stats_df["band"] == band_label)]
        row_vals.append(sub["fraction"].iloc[0] if len(sub)==1 else 0.0)
    frac_matrix.append(row_vals)

frac_matrix = np.array(frac_matrix)

x = np.arange(len(combo_order))
bottom = np.zeros(len(combo_order))
colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a"]

for i, (band_label, _, _) in enumerate(bands):
    ax1.bar(x, frac_matrix[i], bottom=bottom, color=colors[i], label=band_label)
    bottom += frac_matrix[i]

# annotations: show N for each group
for idx, label in enumerate(plot_labels):
    count = total_counts[label]
    ax1.text(idx, 102, f"N={count}", ha='center', va='bottom', fontsize=14)

ax1.set_xticks(x)
ax1.set_xticklabels(plot_labels, rotation=10, ha="right")
ax1.set_ylabel("% of points")
ax1.set_title("e) Lodose (Point count)")
ax1.set_ylim(0, 110)
ax1.grid(axis="y", alpha=0.3)
ax1.legend(title="", loc='lower left')


# ---------------------------------------------------------
# Right PANEL – BOXPLOTS FOR RAIL + ROAD
# (both shown in one panel)
# ---------------------------------------------------------
ax_box.set_facecolor('none') # No panel background

rail_df = df_corr[df_corr["infra"] == "rail"]
road_df = df_corr[df_corr["infra"] == "road"]

rail_psds = rail_df[rail_df["product"]=="PSDS"]["dist"]
rail_egms = rail_df[rail_df["product"]=="EGMS"]["dist"]
road_psds = road_df[road_df["product"]=="PSDS"]["dist"]
road_egms = road_df[road_df["product"]=="EGMS"]["dist"]

# Combined plot structure:
# Groups: Rail-PSDS, Rail-EGMS, Road-PSDS, Road-EGMS
box_data = [rail_psds, rail_egms, road_psds, road_egms]
labels   = ["Rail-PSDS", "Rail-EGMS", "Road-PSDS", "Road-EGMS"]
counts   = [len(rail_psds), len(rail_egms), len(road_psds), len(road_egms)]

bp = ax_box.boxplot(
    box_data,
    labels=labels,
    patch_artist=True,
    boxprops=dict(facecolor="lightgray",alpha=0.7),
    medianprops=dict(color="black")
)

ax_box.set_title("f) Lodose (Distance)")
ax_box.set_ylabel("Distance (m)")
ax_box.grid(axis="y", alpha=0.3)
ax_box.set_xticklabels(labels, rotation=10, ha="right")

# Annotate N below each box
#y_min = ax_box.get_ylim()[0]
#for i, N in enumerate(counts):
#    ax_box.text(i+1, y_min - 0.5, f"{N}", ha="center", va="top", fontsize=14)

# Adjust bottom to fit annotations
#ax_box.set_ylim(y_min - 1, ax_box.get_ylim()[1])

plt.tight_layout()
plt.show()


# ---------------------------------------------------------
# SIMPLE CORRELATION FUNCTION
# ---------------------------------------------------------
def corr(x, y, label):
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    keep = ~(x.isna() | y.isna())
    x = x[keep]
    y = y[keep]
    if len(x) < 5:
        print(f"{label}: insufficient data")
        return
    pear = pearsonr(x, y)[0]
    spear = spearmanr(x, y)[0]
    print(f"\n{label}:")
    print(f"Pearson r  = {pear:.3f}")
    print(f"Spearman ρ = {spear:.3f}")

# Correlations (optional, console only)
corr(rail_df[rail_df["product"]=="PSDS"]["dist"],
     rail_df[rail_df["product"]=="PSDS"]["rmse"],
     "Rail PSDS: dist vs rmse")
corr(rail_df[rail_df["product"]=="EGMS"]["dist"],
     rail_df[rail_df["product"]=="EGMS"]["rmse"],
     "Rail EGMS: dist vs rmse")
corr(road_df[road_df["product"]=="PSDS"]["dist"],
     road_df[road_df["product"]=="PSDS"]["rmse"],
     "Road PSDS: dist vs rmse")
corr(road_df[road_df["product"]=="EGMS"]["dist"],
     road_df[road_df["product"]=="EGMS"]["rmse"],
     "Road EGMS: dist vs rmse")

corr(rail_df[rail_df["product"]=="PSDS"]["dem"],
     rail_df[rail_df["product"]=="PSDS"]["dist"],
     "Rail PSDS: DEM vs dist")
corr(rail_df[rail_df["product"]=="EGMS"]["dem"],
     rail_df[rail_df["product"]=="EGMS"]["dist"],
     "Rail EGMS: DEM vs dist")
corr(road_df[road_df["product"]=="PSDS"]["dem"],
     road_df[road_df["product"]=="PSDS"]["dist"],
     "Road PSDS: DEM vs dist")
corr(road_df[road_df["product"]=="EGMS"]["dem"],
     road_df[road_df["product"]=="EGMS"]["dist"],
     "Road EGMS: DEM vs dist")

# =========================================================
# Two figures: COMBINED 2-PANEL PLOT each
#  2x: 1×2: (RMSE, DEM) × (Rail, Road)
# =========================================================
color_PSDS = "#003399"  # Dark Blue
color_egms = "#CC5500" # Dark Orange (Burnt Orange)
fig_rmse, axes_rmse = plt.subplots(1, 2, figsize=(10, 3), sharey=True, layout="constrained")
fig_dem, axes_dem = plt.subplots(1, 2, figsize=(10, 3), sharey=True)
fig_rmse.patch.set_facecolor('none')
fig_dem.patch.set_facecolor('none')

def plot_scatter_with_trend(ax, x_psds, y_psds, x_egms, y_egms,
                            xlabel, ylabel, title):
    # Scatter
    ax.scatter(x_psds, y_psds, s=8, alpha=0.3, edgecolors="none", color=color_PSDS, label="PSDS")
    ax.scatter(x_egms, y_egms, s=8, alpha=0.3, edgecolors="none", color=color_egms, label="EGMS")

    # Trend lines (guard against empty)
    if len(x_psds) > 5:
        xs, ys = smooth_median(x_psds, y_psds)
        ax.plot(xs, ys, color=color_PSDS, lw=2.5, linestyle="--", zorder=5)
    if len(x_egms) > 5:
        xs, ys = smooth_median(x_egms, y_egms)
        ax.plot(xs, ys, color=color_egms, lw=2.5, linestyle="--", zorder=5)
    # Correlation annotations
    def safe_corr(x, y):
        x = pd.to_numeric(x, errors="coerce")
        y = pd.to_numeric(y, errors="coerce")
        mask = ~(x.isna() | y.isna())
        x = x[mask]
        y = y[mask]
        if len(x) < 5:
            return np.nan, np.nan
        return pearsonr(x, y)[0], spearmanr(x, y)[0]

    r_ps, rho_ps = safe_corr(x_psds, y_psds)
    r_eg, rho_eg = safe_corr(x_egms, y_egms)

    text = (f"PSDS: r={r_ps:.2f}, ρ={rho_ps:.2f}\n"
            f"EGMS: r={r_eg:.2f}, ρ={rho_eg:.2f}")

    ax.text(0.02, 0.98, text,
            transform=ax.transAxes,
            ha="left", va="top",
            fontsize=12, color="black",
            bbox=dict(facecolor="white", alpha=0.6, edgecolor="none"))

    ax.set_title(title)
    ax.set_facecolor('none')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.4)
    ax.legend(loc='upper right')


# ---------------------------------------------------------
# fig_rmse-LEFT: Rail – Distance vs RMSE
# ---------------------------------------------------------
ps = rail_df[rail_df["product"]=="PSDS"]
eg = rail_df[rail_df["product"]=="EGMS"]
plot_scatter_with_trend(
    axes_rmse[0],
    ps["dist"], ps["rmse"],
    eg["dist"], eg["rmse"],
    xlabel="Distance to Rail (m)",
    ylabel="RMSE (mm)",
    title="e) Lodose (Rail)"
)

# ---------------------------------------------------------
# fig_rmse-RIGHT: Road – Distance vs RMSE
# ---------------------------------------------------------
ps = road_df[road_df["product"]=="PSDS"]
eg = road_df[road_df["product"]=="EGMS"]
plot_scatter_with_trend(
    axes_rmse[1],
    ps["dist"], ps["rmse"],
    eg["dist"], eg["rmse"],
    xlabel="Distance to Road(m)",
    ylabel="RMSE (mm)",
    title="f) Lodose (Road)"
)


# ---------------------------------------------------------
# fig_dem-LEFT: Rail – Distance vs DEM
# (FIXED: distance on x-axis, DEM on y-axis)
# ---------------------------------------------------------
ps = rail_df[rail_df["product"]=="PSDS"]
eg = rail_df[rail_df["product"]=="EGMS"]
plot_scatter_with_trend(
    axes_dem[0],
    ps["dist"], ps["dem"],     # <-- FIXED
    eg["dist"], eg["dem"],     # <-- FIXED
    xlabel="Distance to Rail(m)",
    ylabel="DEM (m)",
    title="e) Lodose (Rail)"
)

# ---------------------------------------------------------
# fig_dem-RIGHT: Road – Distance vs DEM
# (FIXED: distance on x-axis, DEM on y-axis)
# ---------------------------------------------------------
ps = road_df[road_df["product"]=="PSDS"]
eg = road_df[road_df["product"]=="EGMS"]
plot_scatter_with_trend(
    axes_dem[1],
    ps["dist"], ps["dem"],     # <-- FIXED
    eg["dist"], eg["dem"],     # <-- FIXED
    xlabel="Distance to Road(m)",
    ylabel="DEM (m)",
    title="f) Lodose (Road)"
)

plt.tight_layout()
plt.show()

