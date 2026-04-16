import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# ---------------------------------------------------------
# CONFIG
# ---------------------------------------------------------
egms_layer = "egms_LOA_coh70_SW99TM_statsSL100"
psds_layer = "PSDS_LO_coh70_SW99TM_statsSL100"
target_id = 1

# ---------------------------------------------------------
# HELPER: Load + prepare dataset
# ---------------------------------------------------------
def load_dataset(layer_name, vel_field, coh_field, count_field):
    layer = QgsProject.instance().mapLayersByName(layer_name)[0]
    
    records = []
    fields = [f.name() for f in layer.fields()]
    for feat in layer.getFeatures():
        rec = dict(zip(fields, feat.attributes()))
        records.append(rec)
    df = pd.DataFrame(records)

    # filter
    df = df[df["id"] == target_id].copy()

    # numeric conversion (keep NaN)
    df["order"]      = pd.to_numeric(df["order"], errors="coerce")
    df[vel_field]    = pd.to_numeric(df[vel_field], errors="coerce")
    df[coh_field]    = pd.to_numeric(df[coh_field], errors="coerce")
    df[count_field]  = pd.to_numeric(df[count_field], errors="coerce")

    df = df.sort_values("order")

    # extract
    chain = df["order"]
    vel   = df[vel_field]
    coh   = df[coh_field]
    npts  = df[count_field]

    # baseline shift at order = 0
    if 0 in df["order"].values:
        v0 = df.loc[df["order"] == 0, vel_field].iloc[0]
        vel = vel - v0

    return chain, vel, coh, npts

# ---------------------------------------------------------
# LOAD BOTH PRODUCTS
# ---------------------------------------------------------
chain_e, vel_e, coh_e, n_e = load_dataset(
    egms_layer, "mean_velocity_median", "rmse_median", "rmse_count"
)

'''chain_e, vel_e, coh_e, n_e = load_dataset(
    egms_layer, "velocity_median", "rmse_median", "rmse_count"
)'''

chain_p, vel_p, coh_p, n_p = load_dataset(
    psds_layer, "velocity_median", "rmse_median", "rmse_count"
)

# ---------------------------------------------------------
# DETERMINE SHARED AXIS LIMITS
# ---------------------------------------------------------
# shared left axis (velocity)
vmin = np.nanmin([vel_e.min(), vel_p.min()])
vmax = np.nanmax([vel_e.max(), vel_p.max()])
# small padding
vpad = 0.05 * (vmax - vmin)
shared_vel_limits = (vmin - vpad, vmax + vpad)

# shared right axis (rmse)
cmin = np.nanmin([coh_e.min(), coh_p.min()])
cmax = np.nanmax([coh_e.max(), coh_p.max()])
cpad = 0.05 * (cmax - cmin)
shared_coh_limits = (cmin - cpad, cmax + cpad)
# ---------------------------------------------------------
# UPDATED COLOR PALETTE (Darker for visibility)
# ---------------------------------------------------------
color_vel = "#003399"  # Dark Blue
color_rmse = "#CC5500" # Dark Orange (Burnt Orange)
color_pts = "#555555"  # Dark Gray/Dim Gray
# ---------------------------------------------------------
# PLOT TWO PANELS
# ---------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 3), sharex=True)
fig.patch.set_facecolor('none')
# ==== PANEL 1: EGMS ====
ax1.set_facecolor('none') # Remove panel background
ax1.set_title("i) EGMS-Stockholm")
ax1.set_xlabel("Chainage index (100 m segments)")
ax1.set_ylabel("Velocity (mm/yr)", color=color_vel)
ax1.plot(chain_e, vel_e, color=color_vel, linestyle="None", marker="o")
ax1.tick_params(axis="y", labelcolor=color_vel)
ax1.grid(True, linestyle="--", alpha=0.4)
ax1.set_ylim(shared_vel_limits)

# Right axis
ax1b = ax1.twinx()
ax1b.set_ylabel("RMSE", color=color_rmse)
ax1b.plot(chain_e, coh_e, color=color_rmse, linestyle="None", marker="s") #, linestyle="--"
ax1b.tick_params(axis="y", labelcolor=color_rmse)
ax1b.set_ylim(shared_coh_limits)

# Bar plot for point count
ax1c = ax1.twinx()
ax1c.set_ylabel("Point Count", color=color_pts)
ax1c.bar(chain_e, n_e, width=0.8, color=color_pts, alpha=0.25)
# offset the axis to the right so labels don't overlap
ax1c.spines["right"].set_position(("outward", 50))
ax1c.tick_params(axis="y", labelcolor=color_pts)

# ==== PANEL 2: PSDS ====
ax2.set_facecolor('none') # Remove panel background
ax2.set_title("j) PSDS-Stockholm")
ax2.set_xlabel("Chainage index (100 m segments)")
ax2.set_ylabel("Velocity (mm/yr)", color=color_vel)
ax2.plot(chain_p, vel_p, color=color_vel, linestyle="None", marker="o")
ax2.tick_params(axis="y", labelcolor=color_vel)
ax2.grid(True, linestyle="--", alpha=0.4)
ax2.set_ylim(shared_vel_limits)

# Right axis
ax2b = ax2.twinx()
ax2b.set_ylabel("RMSE", color=color_rmse)
ax2b.plot(chain_p, coh_p, color=color_rmse, linestyle="None", marker="s")
ax2b.tick_params(axis="y", labelcolor=color_rmse)
ax2b.set_ylim(shared_coh_limits)

# Bar plot for point count
ax2c = ax2.twinx()
ax2c.set_ylabel("Point Count", color=color_pts)
ax2c.bar(chain_p, n_p, width=0.8, color=color_pts, alpha=0.25)
ax2c.spines["right"].set_position(("outward", 50))
ax2c.tick_params(axis="y", labelcolor=color_pts)

# Set a common font size for all tick labels
tick_font_size = 14 

ax1.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax1b.tick_params(axis='y', which='major', labelsize=tick_font_size)
ax1c.tick_params(axis='y', which='major', labelsize=tick_font_size)

ax2.tick_params(axis='both', which='major', labelsize=tick_font_size)
ax2b.tick_params(axis='y', which='major', labelsize=tick_font_size)
ax2c.tick_params(axis='y', which='major', labelsize=tick_font_size)

plt.tight_layout()
plt.show()