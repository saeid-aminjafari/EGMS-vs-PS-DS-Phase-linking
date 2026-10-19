# =============================================================================
# Fig. 3 (+ supplement): velocity, RMSE and point count per 100 m chainage
# segment, EGMS-PS (left) vs PS+DS (right), all corridors in one run.
#                                                               (QGIS console)
#
# Statistics are computed here directly from the POINT layers, so the new
# RMSE field (rmse_c3s, EGMS model) is used without redoing the QGIS joins.
# The segment layer only supplies geometry + 'id' + 'order'.
#
# Check printed per corridor: if the segment layer still holds the old summary
# fields (e.g. rmse_count, velocity_median), the recomputed counts/velocities
# are compared against them. They should match; if not, adjust LINE_BUFFER_M.
# =============================================================================

from qgis.core import (QgsProject, QgsSpatialIndex, QgsGeometry, QgsWkbTypes)
import matplotlib.pyplot as plt
import numpy as np
import os

# ----------------------------- CONFIG ----------------------------------------
# case: (segment layer, segment id to use, EGMS point layer, PS+DS point layer)
CASES = {
    "Kiruna":     ("egms_kir_coh70_SW99TM_statSL100",   1, "egms_kir_coh70_SW99TM",         "PSDS_kir_coh70_SW99TM"),
    "Lilla Edet": ("egms_LE_coh70_SW99TM_statsSL100",   1, "egms_LE_coh70_SW99TM",          "PSDS_LE_coh70_SW99TM"),
    "Lödöse":     ("egms_LOA_coh70_SW99TM_statsSL100",  1, "egms_LOA_coh70_SW99TM",         "PSDS_LOA_coh70_SW99TM"),
    "Gothenburg": ("egms_gbg_coh60_SW99TM_cor50_stats100", 1, "egms_gbg_coh60_SW99TM_cor50",   "PSDS_gbg_coh60_SW99TM_cor50"),
    "Stockholm":  ("egms_sthlm_coh60_SW99TM_cor50_stats100", 1, "egms_sthlm_coh60_SW99TM_cor50", "PSDS_sthlm_coh60_SW99TM_cor50"),
}

# Which corridors go into which figure (panel letters follow this order)
FIGURES = {
    "Fig03": ["Kiruna", "Lilla Edet", "Lödöse", "Gothenburg", "Stockholm"],
}

ID_FIELD, ORDER_FIELD = "id", "order"
VEL_FIELDS  = {"egms": ["mean_velocity", "mean_veloc"], "psds": ["velocity", "vel"]}
RMSE_FIELDS = ["rmse_c3s"]                 # EGMS-model RMSE for both datasets
LINE_BUFFER_M = 50.0     # used only if the segment layer is lines (half-width)

# Old summary fields in the segment layer, used only for the consistency check
CHECK_FIELDS = {"egms": ("mean_velocity_median", "rmse_count"),
                "psds": ("velocity_median", "rmse_count")}

OUT_DIR     = "/Users/saam7695/Desktop/fig3"
FORMATS     = ("pdf", "png")
DPI         = 300
TRANSPARENT = True
FONT        = 12
BOLD        = True
MARKER      = 5.0      # marker size of velocity / RMSE points
ROW_HEIGHT  = 3.0
WIDTH       = 13.0

C_VEL, C_RMSE, C_PTS = "#003399", "#CC5500", "#555555"
C_BAR, BAR_ALPHA = "#9E9E9E", 0.45      # light bars so markers stay readable

# Isolated spikes (a segment far from its neighbours) would stretch the axis
# and flatten every other corridor signal. They are kept but drawn as
# triangles on the axis edge with their value. Real deformation that spans
# several segments is NOT a spike and stays on scale.
SPIKE_WINDOW = 5         # segments in the running median
SPIKE_K      = 6.0       # threshold in robust standard deviations
SPIKE_MIN    = 5.0       # mm/yr: never call a deviation below this a spike
RMSE_PCT     = 99.0      # RMSE axis top = this percentile of all segments
LABELS = {"egms": "EGMS-PS", "psds": "PS+DS"}
# -----------------------------------------------------------------------------


# ============================ data (QGIS) ====================================

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


def load_points(lyr, vel_f, rmse_f):
    """Spatial index + attribute arrays, read once per layer."""
    idx = QgsSpatialIndex()
    geom, vel, rmse = {}, {}, {}
    for f in lyr.getFeatures():
        g = f.geometry()
        if g is None or g.isEmpty():
            continue
        idx.addFeature(f)
        geom[f.id()] = g
        vel[f.id()] = num(f[vel_f])
        rmse[f.id()] = num(f[rmse_f])
    return idx, geom, vel, rmse


def segment_stats(seg_geom, pts):
    idx, geom, vel, rmse = pts
    cand = idx.intersects(seg_geom.boundingBox())
    inside = [i for i in cand if seg_geom.intersects(geom[i])]
    v = np.array([vel[i] for i in inside], float)
    r = np.array([rmse[i] for i in inside], float)
    ok = ~np.isnan(r)
    return (np.nanmedian(v) if np.isfinite(v).any() else np.nan,
            np.median(r[ok]) if ok.any() else np.nan,
            int(ok.sum()))


def corridor_profiles(seg_name, seg_id, egms_name, psds_name):
    seg = layer(seg_name)
    is_line = QgsWkbTypes.geometryType(seg.wkbType()) == QgsWkbTypes.LineGeometry
    out = {}
    for key, name in (("egms", egms_name), ("psds", psds_name)):
        pl = layer(name)
        if pl.crs() != seg.crs():
            print("   WARNING: CRS differs (%s vs %s)" % (pl.name(), seg.name()))
        pts = load_points(pl, field(pl, VEL_FIELDS[key]), field(pl, RMSE_FIELDS))
        rows, check = [], []
        for f in seg.getFeatures():
            if num(f[ID_FIELD]) != seg_id:
                continue
            g = f.geometry()
            if is_line:
                g = g.buffer(LINE_BUFFER_M, 8)
            vmed, rmed, n = segment_stats(g, pts)
            rows.append((num(f[ORDER_FIELD]), vmed, rmed, n))
            ov, oc = CHECK_FIELDS[key]
            if seg.fields().indexOf(ov) >= 0 and seg.fields().indexOf(oc) >= 0:
                check.append((vmed, num(f[ov]), n, num(f[oc])))
        if not rows:
            ids = sorted({num(f[ID_FIELD]) for f in seg.getFeatures()})
            raise RuntimeError("%s: no segments with %s = %s. Available ids: %s"
                               % (seg.name(), ID_FIELD, seg_id, ids))
        a = np.array(sorted(rows), float)               # order, vel, rmse, n
        # reference velocity: first segment with data (order 0 if present)
        v0 = a[np.isfinite(a[:, 1]), 1][0] if np.isfinite(a[:, 1]).any() else 0.0
        a[:, 1] -= v0
        out[key] = a
        msg = "   %-7s %d segments, %d points" % (LABELS[key], len(a), int(np.nansum(a[:, 3])))
        if check:
            c = np.array(check, float)
            ok = np.isfinite(c).all(axis=1)
            msg += (" | check vs old summary: velocity |diff| median %.2f mm/yr, "
                    "count diff median %+.0f"
                    % (np.median(np.abs(c[ok, 0] - c[ok, 1])) if ok.any() else np.nan,
                       np.median(c[ok, 2] - c[ok, 3]) if ok.any() else np.nan))
        print(msg)
    return out


# ============================ plotting =======================================

def lims(*arrays, pad=0.05, floor=None):
    v = np.concatenate([np.asarray(x, float).ravel() for x in arrays])
    v = v[np.isfinite(v)]
    lo, hi = v.min(), v.max()
    p = pad * (hi - lo if hi > lo else 1.0)
    lo = lo - p if floor is None else max(floor, lo - p)
    return lo, hi + p


def spikes(v):
    """True where a segment deviates from the running median of its neighbours."""
    v = np.asarray(v, float)
    h = SPIKE_WINDOW // 2
    med = np.array([np.nanmedian(v[max(0, i - h):i + h + 1])
                    if np.isfinite(v[max(0, i - h):i + h + 1]).any() else np.nan
                    for i in range(len(v))])
    r = v - med
    ok = np.isfinite(r)
    if ok.sum() < 5:
        return np.zeros(len(v), bool)
    sig = 1.4826 * np.median(np.abs(r[ok] - np.median(r[ok])))
    thr = max(SPIKE_K * sig, SPIKE_MIN)
    return ok & (np.abs(r) > thr)


def plot_clipped(ax, x, y, lim, color, marker):
    """Plot y; values outside lim become edge triangles labelled with the value."""
    x, y = np.asarray(x), np.asarray(y)
    lo, hi = lim
    inside = (y >= lo) & (y <= hi)
    ax.plot(x[inside], y[inside], ls="None", marker=marker, ms=MARKER, color=color)
    for xi, yi in zip(x[~inside & np.isfinite(y)], y[~inside & np.isfinite(y)]):
        top = yi > hi
        ye = hi if top else lo
        ax.plot(xi, ye, ls="None", marker="^" if top else "v", ms=MARKER + 2.5,
                color=color, clip_on=False, zorder=5)
        ax.annotate("%.0f" % yi, (xi, ye), xytext=(4, -9 if top else 3),
                    textcoords="offset points", fontsize=FONT - 3, color=color,
                    fontweight="bold" if BOLD else "normal")


def plot_figure(name, cases, data):
    if BOLD:
        plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                             "axes.titleweight": "bold"})
    n = len(cases)
    fig, axes = plt.subplots(n, 2, figsize=(WIDTH, ROW_HEIGHT * n), squeeze=False)
    if TRANSPARENT:
        fig.patch.set_alpha(0)
    # one RMSE scale for the whole figure, so corridors are comparable
    allr = np.concatenate([data[c][k][:, 2] for c in cases for k in ("egms", "psds")])
    allr = allr[np.isfinite(allr)]
    rlim = (0, 1.08 * np.percentile(allr, RMSE_PCT))
    letter = iter("abcdefghijklmnopqrstuvwxyz")
    for row, case in enumerate(cases):
        d = data[case]
        inl = [d[k][~spikes(d[k][:, 1]), 1] for k in ("egms", "psds")]
        vlim = lims(*inl, pad=0.08)
        nmax = 1.1 * max(np.nanmax(d["egms"][:, 3]), np.nanmax(d["psds"][:, 3]), 1)
        for col, key in enumerate(("egms", "psds")):
            a = d[key]
            ax = axes[row, col]
            ax.set_facecolor("none")
            axn = ax.twinx()                       # point count (bars, behind)
            axr = ax.twinx()                       # RMSE
            axn.spines["right"].set_position(("outward", 50))
            ax.set_zorder(axr.get_zorder() + 1)    # velocity on top
            ax.patch.set_visible(False)

            axn.bar(a[:, 0], a[:, 3], width=0.8, color=C_BAR, alpha=BAR_ALPHA, lw=0)
            plot_clipped(axr, a[:, 0], a[:, 2], rlim, C_RMSE, "s")
            plot_clipped(ax, a[:, 0], a[:, 1], vlim, C_VEL, "o")

            ax.set_ylim(vlim); axr.set_ylim(rlim); axn.set_ylim(0, nmax)
            ax.grid(True, ls="--", alpha=0.4)
            ax.set_title("%s) %s – %s" % (next(letter), LABELS[key], case), fontsize=FONT)
            ax.set_ylabel("Velocity (mm/yr)", color=C_VEL, fontsize=FONT)
            axr.set_ylabel("RMSE (mm)", color=C_RMSE, fontsize=FONT)
            axn.set_ylabel("Point count", color=C_PTS, fontsize=FONT)
            if row == n - 1:
                ax.set_xlabel("Chainage index (100 m segments)", fontsize=FONT)
            for a_, c_ in ((ax, C_VEL), (axr, C_RMSE), (axn, C_PTS)):
                a_.tick_params(axis="y", labelcolor=c_, labelsize=FONT)
            ax.tick_params(axis="x", labelsize=FONT)
    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    for ext in FORMATS:
        p = os.path.join(OUT_DIR, "%s.%s" % (name, ext))
        fig.savefig(p, dpi=DPI, transparent=TRANSPARENT, bbox_inches="tight")
        print("Written: %s" % p)
    plt.show()


# ============================ run ============================================
needed = [c for cs in FIGURES.values() for c in cs]
data = {}
for case in dict.fromkeys(needed):
    print("CASE: %s" % case)
    data[case] = corridor_profiles(*CASES[case])
for name, cases in FIGURES.items():
    plot_figure(name, cases, data)