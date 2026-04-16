import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# ---------------------------------------------------------
# CONFIG
# ---------------------------------------------------------
egms_layer = "egms_sthlm_coh60_SW99TM_cor50_stats100"
psds_layer = "PSDS_sthlm_coh60_SW99TM_cor50_stats100"

# EGMS fields
EGMS_RMSE = "rmse_median"
EGMS_velo  = "mean_velocity_median"
EGMS_W    = "mean_velocity_count"     # coverage weight

# PSDS fields
PSDS_RMSE = "rmse_median"
PSDS_velo  = "velocity_median"
PSDS_W    = "velocity_count"

# ---------------------------------------------------------
# HELPER: Build weighted sample array
# ---------------------------------------------------------
def load_weighted_values(layer_name, value_field, weight_field):

    layer = QgsProject.instance().mapLayersByName(layer_name)[0]

    # extract table
    records = []
    fields = [f.name() for f in layer.fields()]
    for feat in layer.getFeatures():
        rec = dict(zip(fields, feat.attributes()))
        records.append(rec)

    df = pd.DataFrame(records)

    # convert numeric (keep NaN)
    df[value_field]  = pd.to_numeric(df[value_field], errors="coerce")
    df[weight_field] = pd.to_numeric(df[weight_field], errors="coerce")

    # remove invalid weights (<= 0)
    df = df[df[weight_field] > 0]

    # drop NaN values
    df = df.dropna(subset=[value_field, weight_field])

    # build weighted sample:
    # repeat each value according to its coverage weight
    values = df[value_field].to_numpy()
    weights = df[weight_field].to_numpy().astype(int)

    weighted_samples = np.repeat(values, weights)

    return weighted_samples

# ---------------------------------------------------------
# LOAD WEIGHTED DATA
# ---------------------------------------------------------
egms_rmse_w = load_weighted_values(egms_layer, EGMS_RMSE, EGMS_W)
egms_velo_w  = load_weighted_values(egms_layer, EGMS_velo,  EGMS_W)

psds_rmse_w = load_weighted_values(psds_layer, PSDS_RMSE, PSDS_W)
psds_velo_w  = load_weighted_values(psds_layer, PSDS_velo,  PSDS_W)

# ---------------------------------------------------------
# HELPER: Compute CDF
# ---------------------------------------------------------
def compute_cdf(arr):
    arr_sorted = np.sort(arr)
    cdf_vals = np.arange(1, len(arr_sorted)+1) / len(arr_sorted)
    return arr_sorted, cdf_vals

# Compute CDFs
egms_rmse_x, egms_rmse_y = compute_cdf(egms_rmse_w)
psds_rmse_x, psds_rmse_y = compute_cdf(psds_rmse_w)

egms_velo_x,  egms_velo_y  = compute_cdf(egms_velo_w)
psds_velo_x,  psds_velo_y  = compute_cdf(psds_velo_w)

# ---------------------------------------------------------
# PLOT
# ---------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 3))
fig.patch.set_facecolor('none')
# ==== PANEL 1: WEIGHTED RMSE CDF ====
ax1.plot(egms_rmse_x, egms_rmse_y, label="EGMS weighted RMSE", color="tab:blue")
ax1.plot(psds_rmse_x, psds_rmse_y, label="PSDS weighted RMSE", color="tab:orange")
ax1.set_facecolor('none') # Remove panel background
ax1.set_title("i) Stockholm (point-count-weighted RMSE)")
ax1.set_xlabel("RMSE (mm)")
ax1.set_ylabel("Cumulative Probability")
ax1.grid(True, linestyle="--", alpha=0.4)
ax1.legend()

# ==== PANEL 2: WEIGHTED VELOCITY  CDF ====
ax2.plot(egms_velo_x, egms_velo_y, label="EGMS weighted velocity", color="tab:blue")
ax2.plot(psds_velo_x, psds_velo_y, label="PSDS weighted velocity", color="tab:orange")
ax2.set_facecolor('none') # Remove panel background
ax2.set_title("j) Stockholm (point-count-weighted velocity)")
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