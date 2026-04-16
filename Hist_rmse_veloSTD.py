import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# ---------------------------------------------------------
# CONFIG
# ---------------------------------------------------------
egms_layer = "egms_kir_coh70_SW99TM_statGRD100"
psds_layer = "PSDS_kir_coh70_SW99TM_statGRD100"

# Field names
EGMS_RMSE = "rmse_median"
EGMS_velo  = "mean_veloc_median" #or stddev

PSDS_RMSE = "rmse_median"
PSDS_velo = "velocity_median"

# ---------------------------------------------------------
# HELPER: Load fields for whole layer
# ---------------------------------------------------------
def load_metrics(layer_name, rmse_field, velo_field):
    layer = QgsProject.instance().mapLayersByName(layer_name)[0]

    records = []
    fields = [f.name() for f in layer.fields()]

    for feat in layer.getFeatures():
        rec = dict(zip(fields, feat.attributes()))
        records.append(rec)

    df = pd.DataFrame(records)

    # numeric conversion (keep NaN)
    df[rmse_field] = pd.to_numeric(df[rmse_field], errors="coerce")
    df[velo_field]  = pd.to_numeric(df[velo_field],  errors="coerce")

    # remove NaN for analysis
    rmse = df[rmse_field].dropna().to_numpy()
    velo = df[velo_field].dropna().to_numpy()

    return rmse, velo

# ---------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------
egms_rmse, egms_velo = load_metrics(egms_layer, EGMS_RMSE, EGMS_velo)
psds_rmse, psds_velo = load_metrics(psds_layer, PSDS_RMSE, PSDS_velo)

# ---------------------------------------------------------
# HELPER: Compute CDF
# ---------------------------------------------------------
def compute_cdf(arr):
    arr_sorted = np.sort(arr)
    cdf_vals = np.arange(1, len(arr_sorted)+1) / len(arr_sorted)
    return arr_sorted, cdf_vals

# EGMS CDFs
egms_rmse_x, egms_rmse_y = compute_cdf(egms_rmse)
egms_velo_x,  egms_velo_y  = compute_cdf(egms_velo)

# PSDS CDFs
psds_rmse_x, psds_rmse_y = compute_cdf(psds_rmse)
psds_velo_x,  psds_velo_y  = compute_cdf(psds_velo)

# ---------------------------------------------------------
# PLOT
# ---------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 3))
fig.patch.set_facecolor('none')
# ==== PANEL 1: RMSE CDF ====
ax1.plot(egms_rmse_x, egms_rmse_y, label="EGMS RMSE", color="tab:blue")
ax1.plot(psds_rmse_x, psds_rmse_y, label="PSDS RMSE", color="tab:orange")
ax1.set_facecolor('none') # Remove panel background
ax1.set_title("a) Kiruna (RMSE)")
ax1.set_xlabel("RMSE (mm)")
ax1.set_ylabel("Cumulative Probability")
ax1.grid(True, linestyle="--", alpha=0.4)
ax1.legend()

# ==== PANEL 2: Velocity CDF ====
ax2.plot(egms_velo_x, egms_velo_y, label="EGMS velocity", color="tab:blue")
ax2.plot(psds_velo_x, psds_velo_y, label="PSDS velocity", color="tab:orange")
ax2.set_facecolor('none') # Remove panel background
ax2.set_title("b) Kiruna (Velocity)")
ax2.set_xlabel("Velocity (mm/yr)")
ax2.set_ylabel("Cumulative Probability")
ax2.grid(True, linestyle="--", alpha=0.4)
ax2.legend()

# Set a common font size for all tick labels
tick_font_size = 14 

ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)

plt.tight_layout()
plt.show()