# =============================================================================
# Fig. 4: cumulative distributions of RMSE (left) and velocity (right),
# EGMS-PS vs PS+DS, all corridors in one figure.                (QGIS console)
#
# Computed directly from the POINT layers (rmse_c3s = EGMS-model RMSE), so the
# QGIS grid-statistics layers do not have to be rebuilt.
#   MODE = "grid"  : median per 100 m x 100 m cell (as in the original Fig. 4)
#   MODE = "points": every measurement point
# =============================================================================

from qgis.core import QgsProject
import matplotlib.pyplot as plt
import numpy as np
import os

# ----------------------------- CONFIG ----------------------------------------
CASES = {
    "Kiruna":     ("egms_kir_coh70_SW99TM",         "PSDS_kir_coh70_SW99TM"),
    "Lilla Edet": ("egms_LE_coh70_SW99TM",          "PSDS_LE_coh70_SW99TM"),
    "Lödöse":     ("egms_LOA_coh70_SW99TM",         "PSDS_LOA_coh70_SW99TM"),
    "Gothenburg": ("egms_gbg_coh60_SW99TM_cor50",   "PSDS_gbg_coh60_SW99TM_cor50"),
    "Stockholm":  ("egms_sthlm_coh60_SW99TM_cor50", "PSDS_sthlm_coh60_SW99TM_cor50"),
}
ORDER = ["Kiruna", "Lilla Edet", "Lödöse", "Gothenburg", "Stockholm"]

MODE  = "points"           # "grid" or "points"
CELL  = 100.0            # m, grid cell size
VEL_FIELDS  = {"egms": ["mean_velocity", "mean_veloc"], "psds": ["velocity", "vel"]}
RMSE_FIELDS = ["rmse_c3s"]

# PS+DS is relative to a scene reference point, EGMS is GNSS-calibrated. With
# ALIGN = True the PS+DS velocities are shifted by the difference of the two
# medians so the SHAPES of the distributions can be compared (the offset itself
# is the datum difference discussed in Sect. 3.1). The shift is printed.
ALIGN = True

RMSE_REF  = 4.0          # mm, vertical reference line + share annotated
VEL_PCT   = (0.5, 99.5)  # velocity x-range per corridor (long tails clipped)
RMSE_PCT  = 99.5         # common RMSE x-range top

OUT_DIR, NAME = "/Users/saam7695/Desktop/fig4", "Fig04"
FORMATS, DPI, TRANSPARENT = ("pdf", "png"), 300, True
FONT, BOLD, LW = 12, True, 2.0
ROW_HEIGHT, WIDTH = 2.8, 12.0
COL = {"egms": "tab:blue", "psds": "tab:orange"}
LABELS = {"egms": "EGMS-PS", "psds": "PS+DS"}
# -----------------------------------------------------------------------------


def layer(name):
    ls = QgsProject.instance().mapLayersByName(name)
    if not ls:
        raise RuntimeError("Layer not found: %s" % name)
    return ls[0]


def field(lyr, candidates):
    names = [f.name() for f in lyr.fields()]
    for c in candidates:
        if c in names:
            return c
    raise RuntimeError("%s: none of %s found" % (lyr.name(), candidates))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def read_points(name, key):
    lyr = layer(name)
    vf, rf = field(lyr, VEL_FIELDS[key]), field(lyr, RMSE_FIELDS)
    xy, vel, rmse = [], [], []
    for f in lyr.getFeatures():
        g = f.geometry()
        if g is None or g.isEmpty():
            continue
        p = g.asPoint() if not g.isMultipart() else g.asMultiPoint()[0]
        xy.append((p.x(), p.y())); vel.append(num(f[vf])); rmse.append(num(f[rf]))
    return np.array(xy), np.array(vel), np.array(rmse)


def grid_medians(xy, values):
    """Median of each value array per CELL x CELL square (cells keyed by index)."""
    ij = np.floor(xy / CELL).astype(np.int64)
    key = ij[:, 0] * 10_000_000 + ij[:, 1]
    order = np.argsort(key)
    key_s = key[order]
    bounds = np.flatnonzero(np.diff(key_s)) + 1
    out = []
    for v in values:
        groups = np.split(v[order], bounds)
        out.append(np.array([np.nanmedian(g) if np.isfinite(g).any() else np.nan
                             for g in groups]))
    return out


def cdf(a):
    a = np.sort(a[np.isfinite(a)])
    return a, np.arange(1, len(a) + 1) / len(a)


# ----------------------------- data ------------------------------------------
data = {}
for case in ORDER:
    d = {}
    for key, name in zip(("egms", "psds"), CASES[case]):
        xy, v, r = read_points(name, key)
        if MODE == "grid":
            v, r = grid_medians(xy, (v, r))
        d[key] = {"vel": v, "rmse": r}
    shift = 0.0
    if ALIGN:
        shift = np.nanmedian(d["egms"]["vel"]) - np.nanmedian(d["psds"]["vel"])
        d["psds"]["vel"] = d["psds"]["vel"] + shift
    d["shift"] = shift
    data[case] = d
    s = ["%s: %s n=%d, median RMSE %.2f mm, <=%.0f mm %.0f%%"
         % (case if k == "egms" else " " * len(case), LABELS[k],
            np.isfinite(d[k]["rmse"]).sum(), np.nanmedian(d[k]["rmse"]), RMSE_REF,
            100 * np.nanmean(d[k]["rmse"][np.isfinite(d[k]["rmse"])] <= RMSE_REF))
         for k in ("egms", "psds")]
    print("\n".join(s) + ("\n   PS+DS velocity shifted by %+.2f mm/yr" % shift if ALIGN else ""))


# ----------------------------- plot ------------------------------------------
if BOLD:
    plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                         "axes.titleweight": "bold"})
n = len(ORDER)
fig, axes = plt.subplots(n, 2, figsize=(WIDTH, ROW_HEIGHT * n), squeeze=False)
if TRANSPARENT:
    fig.patch.set_alpha(0)
allr = np.concatenate([data[c][k]["rmse"] for c in ORDER for k in ("egms", "psds")])
rmax = np.nanpercentile(allr, RMSE_PCT) * 1.05
letters = iter("abcdefghijklmnopqrstuvwxyz")
unit = "cells" if MODE == "grid" else "points"

for row, case in enumerate(ORDER):
    d = data[case]
    axr, axv = axes[row]
    for key in ("egms", "psds"):
        x, y = cdf(d[key]["rmse"])
        axr.plot(x, y, color=COL[key], lw=LW, label=LABELS[key])
        x, y = cdf(d[key]["vel"])
        axv.plot(x, y, color=COL[key], lw=LW, label=LABELS[key])
    # RMSE panel
    axr.axvline(RMSE_REF, color="0.4", ls=":", lw=1.2)
    share = ["%s %.0f%%" % (LABELS[k], 100 * np.nanmean(
        d[k]["rmse"][np.isfinite(d[k]["rmse"])] <= RMSE_REF)) for k in ("egms", "psds")]
    axr.text(0.98, 0.06, "%s ≤ %.0f mm: %s" % (unit, RMSE_REF, ", ".join(share)),
             transform=axr.transAxes, ha="right", va="bottom", fontsize=FONT - 2)
    axr.set_xlim(0, rmax)
    axr.set_title("%s) %s (RMSE)" % (next(letters), case), fontsize=FONT)
    # velocity panel
    allv = np.concatenate([d[k]["vel"] for k in ("egms", "psds")])
    lo, hi = np.nanpercentile(allv, VEL_PCT)
    pad = 0.05 * (hi - lo)
    axv.set_xlim(lo - pad, hi + pad)
    axv.set_title("%s) %s (velocity)" % (next(letters), case), fontsize=FONT)
    for ax in (axr, axv):
        ax.set_facecolor("none")
        ax.set_ylim(0, 1.02)
        ax.grid(True, ls="--", alpha=0.4)
        ax.tick_params(labelsize=FONT)
        ax.set_ylabel("Cumulative prob.", fontsize=FONT)
    if row == n - 1:
        axr.set_xlabel("RMSE (mm)", fontsize=FONT)
        axv.set_xlabel("Velocity (mm/yr)" + (" — PS+DS aligned to EGMS median" if ALIGN else ""),
                       fontsize=FONT)

axes[0, 0].legend(loc="upper left", fontsize=FONT - 1, frameon=False)
fig.tight_layout()
os.makedirs(OUT_DIR, exist_ok=True)
for ext in FORMATS:
    p = os.path.join(OUT_DIR, "%s.%s" % (NAME, ext))
    fig.savefig(p, dpi=DPI, transparent=TRANSPARENT, bbox_inches="tight")
    print("Written: %s" % p)
plt.show()