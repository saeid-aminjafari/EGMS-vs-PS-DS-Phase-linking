# =============================================================================
# Fig. 5: RMSE per land-cover class (left) and point share + class area
# (right), EGMS-PS vs PS+DS, all corridors in one figure.      (QGIS console)
#
# Land cover is sampled here directly from the NMD raster at every point of the
# BASE point layers (which hold rmse_c3s), so the *_statsLC_NMD layers are not
# needed. The NMD raster of each corridor is found automatically among the
# loaded raster layers whose name contains "NMD" (the one covering the points
# with the smallest extent). Override in RASTERS if the choice printed is wrong.
# =============================================================================

from qgis.core import (QgsProject, QgsRasterLayer, QgsCoordinateTransform,
                       QgsPointXY)
from osgeo import gdal
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

RASTERS = {                     # leave None for automatic choice
    "Kiruna": None, "Lilla Edet": None, "Lödöse": "NMD2023bas_v0_2_LO",
    "Gothenburg": None, "Stockholm": None,
}
RASTER_KEY = "nmd"              # substring used to find land-cover rasters

RMSE_FIELDS = ["rmse_c3s"]
MIN_PTS   = 10                  # classes with fewer points in BOTH datasets are dropped
NODATA    = {0}                 # land-cover codes to ignore (outside corridor)

# x-axis labels: taken from the raster's own legend (symbology), so they
# match the NMD version you use. Shorten any of them here: {code: "label"}.
LABEL_OVERRIDE = {                 # NMD 2023 basic layer, English short names
    3: "Arable land",
    20: "Open wetland",
    23: "Mountain forest, wetland",
    41: "Open land, no veg.",
    42: "Open land, veg.",
    43: "Mountain forest",
    51: "Buildings",
    52: "Artificial, other",
    53: "Roads/railways",
    54: "Peat extraction",
    61: "Inland water",
    62: "Sea",
    110: "Forest",
    118: "Temp. non-forest",
    120: "Forest, wetland",
    128: "Temp. non-forest, wetland",
}
WRAP = 16                       # wrap long legend labels at this many characters

OUT_DIR, NAME = "/Users/saam7695/Desktop/fig5", "Fig05"
FORMATS, DPI, TRANSPARENT = ("pdf", "png"), 300, True
FONT, BOLD = 12, True
ROW_HEIGHT, WIDTH = 3.4, 14.0
COL = {"egms": "tab:blue", "psds": "tab:orange"}
LABELS = {"egms": "EGMS-PS", "psds": "PS+DS"}
# -----------------------------------------------------------------------------


def layer(name):
    ls = QgsProject.instance().mapLayersByName(name)
    if not ls:
        raise RuntimeError("Layer not found: %s" % name)
    return ls[0]


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def read_points(name):
    lyr = layer(name)
    rf = next((c for c in RMSE_FIELDS if lyr.fields().indexOf(c) >= 0), None)
    if rf is None:
        raise RuntimeError("%s: none of %s found" % (name, RMSE_FIELDS))
    pts, rmse = [], []
    for f in lyr.getFeatures():
        g = f.geometry()
        if g is None or g.isEmpty():
            continue
        p = g.asPoint() if not g.isMultipart() else g.asMultiPoint()[0]
        pts.append(QgsPointXY(p)); rmse.append(num(f[rf]))
    return lyr, pts, np.array(rmse)


class Raster:
    """Land-cover raster read once with GDAL, sampled by index (fast)."""
    def __init__(self, rl):
        self.rl = rl
        path = rl.dataProvider().dataSourceUri().split("|")[0]
        ds = gdal.Open(path)
        self.gt = ds.GetGeoTransform()
        self.arr = ds.GetRasterBand(1).ReadAsArray()
        self.px_area = abs(self.gt[1] * self.gt[5])

    def sample(self, pts, src_crs):
        tr = None
        if src_crs != self.rl.crs():
            tr = QgsCoordinateTransform(src_crs, self.rl.crs(), QgsProject.instance())
        xy = np.array([(q.x(), q.y()) for q in
                       ((tr.transform(p) if tr else p) for p in pts)])
        c = np.floor((xy[:, 0] - self.gt[0]) / self.gt[1]).astype(int)
        r = np.floor((xy[:, 1] - self.gt[3]) / self.gt[5]).astype(int)
        ok = (r >= 0) & (r < self.arr.shape[0]) & (c >= 0) & (c < self.arr.shape[1])
        out = np.full(len(xy), -1, int)
        out[ok] = self.arr[r[ok], c[ok]]
        return out

    def area_km2(self):
        u, n = np.unique(self.arr, return_counts=True)
        return {int(k): v * self.px_area / 1e6 for k, v in zip(u, n)}

    def legend(self):
        out = {}
        rend = self.rl.renderer()
        if hasattr(rend, "classes"):
            for c in rend.classes():
                try:
                    out[int(c.value)] = c.label
                except Exception:
                    pass
        return out


def pick_raster(case, lyr, pts):
    if RASTERS.get(case):
        return layer(RASTERS[case])
    cands = [l for l in QgsProject.instance().mapLayers().values()
             if isinstance(l, QgsRasterLayer) and RASTER_KEY in l.name().lower()]
    best = None
    sample = pts[::max(1, len(pts) // 500)]
    for rl in cands:
        tr = QgsCoordinateTransform(lyr.crs(), rl.crs(), QgsProject.instance())
        ext = rl.extent()
        inside = np.mean([ext.contains(tr.transform(p)) for p in sample])
        if inside >= 0.9:
            size = ext.width() * ext.height()
            if best is None or size < best[0]:
                best = (size, rl)
    if best is None:
        raise RuntimeError("%s: no '%s' raster covers the points; set RASTERS"
                           % (case, RASTER_KEY))
    return best[1]


# ----------------------------- data ------------------------------------------
data, legends = {}, {}
for case in ORDER:
    e_lyr, e_pts, e_r = read_points(CASES[case][0])
    p_lyr, p_pts, p_r = read_points(CASES[case][1])
    rl = pick_raster(case, e_lyr, e_pts)
    R = Raster(rl)
    if R.arr.size > 2e8:
        print("   WARNING: %s is large (%.0f Mpx); area per class covers the "
              "whole raster, not only the corridor" % (rl.name(), R.arr.size / 1e6))
    e_lc, p_lc = R.sample(e_pts, e_lyr.crs()), R.sample(p_pts, p_lyr.crs())
    area = R.area_km2()
    legends.update(R.legend())
    d = {"egms": (e_lc, e_r), "psds": (p_lc, p_r)}
    classes = sorted({int(c) for k in d for c in d[k][0]} - NODATA - {-1})
    keep = [c for c in classes
            if max((d[k][0] == c).sum() for k in d) >= MIN_PTS]
    dropped = sorted(set(classes) - set(keep))
    keep.sort(key=lambda c: area.get(c, 0), reverse=True)
    data[case] = {"classes": keep, "area": [area.get(c, 0) for c in keep], "d": d}
    print("%-11s raster %s | classes %s%s" % (case, rl.name(), keep,
          (" | dropped (<%d pts): %s" % (MIN_PTS, dropped)) if dropped else ""))

print("\nClass labels used (code: label) - shorten via LABEL_OVERRIDE if needed:")
for c in sorted({c for v in data.values() for c in v["classes"]}):
    print("  %4d: %s" % (c, LABEL_OVERRIDE.get(c, legends.get(c, "(no label in symbology)"))))


# ----------------------------- plot ------------------------------------------
import textwrap


def tick(c):
    lab = LABEL_OVERRIDE.get(c) or legends.get(c) or str(c)
    return "\n".join(textwrap.wrap(str(lab), WRAP)) or str(c)


if BOLD:
    plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                         "axes.titleweight": "bold"})
n = len(ORDER)
fig, axes = plt.subplots(n, 2, figsize=(WIDTH, ROW_HEIGHT * n), squeeze=False)
if TRANSPARENT:
    fig.patch.set_alpha(0)
allr = np.concatenate([data[c]["d"][k][1] for c in ORDER for k in ("egms", "psds")])
rtop = np.nanpercentile(allr, 99.5) * 1.05
letters = iter("abcdefghijklmnopqrstuvwxyz")

for row, case in enumerate(ORDER):
    v = data[case]
    cls, x = v["classes"], np.arange(len(v["classes"]))
    axb, axs = axes[row]

    # RMSE boxplots
    for off, key in ((-0.2, "egms"), (0.2, "psds")):
        lc, r = v["d"][key]
        groups = [r[(lc == c) & np.isfinite(r)] for c in cls]
        pos = [xi + off for xi, g in zip(x, groups) if len(g)]
        groups = [g for g in groups if len(g)]
        if groups:
            axb.boxplot(groups, positions=pos, widths=0.35, patch_artist=True,
                        boxprops=dict(facecolor=COL[key], edgecolor="black"),
                        medianprops=dict(color="black", lw=1.5),
                        flierprops=dict(marker=".", ms=2, alpha=0.4,
                                        markerfacecolor=COL[key], markeredgecolor="none"))
    axb.set_ylim(0, rtop)
    axb.set_ylabel("RMSE (mm)", fontsize=FONT)
    axb.set_title("%s) %s (RMSE)" % (next(letters), case), fontsize=FONT)

    # point share + area
    bw, smax = 0.38, 0.0
    for off, key in ((-bw / 2, "egms"), (bw / 2, "psds")):
        lc = v["d"][key][0]
        valid = ~np.isin(lc, list(NODATA) + [-1])
        tot = max(valid.sum(), 1)
        cnt = [(lc == c).sum() for c in cls]
        share = [100 * k / tot for k in cnt]
        smax = max([smax] + share)
        bars = axs.bar(x + off, share, bw, color=COL[key])
        for b, k in zip(bars, cnt):
            axs.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.8, str(k),
                     ha="center", va="bottom", fontsize=FONT - 4, rotation=90)
    axs.set_ylim(0, smax * 1.3 + 1)
    axs.set_ylabel("% of points", fontsize=FONT)
    axa = axs.twinx()
    axa.plot(x, v["area"], color="black", marker="o", lw=1.2)
    axa.set_ylim(0, max(v["area"]) * 1.15 if v["area"] else 1)
    axa.set_ylabel("Area (km²)", fontsize=FONT)
    axa.tick_params(labelsize=FONT)
    axs.set_title("%s) %s (point share)" % (next(letters), case), fontsize=FONT)

    for ax in (axb, axs):
        ax.set_facecolor("none")
        ax.set_xticks(x)
        ax.set_xticklabels([tick(c) for c in cls], rotation=35, ha="right",
                           fontsize=FONT - 2)
        ax.tick_params(axis="y", labelsize=FONT)
        ax.grid(axis="y", ls="--", alpha=0.4)

handles = [mpatches.Patch(facecolor=COL[k], edgecolor="black", label=LABELS[k])
           for k in ("egms", "psds")]
handles.append(plt.Line2D([], [], color="black", marker="o", label="Area (km²)"))
fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
           fontsize=FONT, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.985))
os.makedirs(OUT_DIR, exist_ok=True)
for ext in FORMATS:
    p = os.path.join(OUT_DIR, "%s.%s" % (NAME, ext))
    fig.savefig(p, dpi=DPI, transparent=TRANSPARENT, bbox_inches="tight")
    print("Written: %s" % p)
plt.show()