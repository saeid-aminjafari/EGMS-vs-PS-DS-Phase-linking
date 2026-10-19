# =============================================================================
# EGMS-PS  vs  PS+DS   velocity discrepancy diagnostic  (v2)
# Run from the QGIS Python Console (Plugins > Python Console > Show Editor)
#
# NOTE: layers are assumed UNANCHORED - PS+DS is relative to an arbitrary
# reference point, EGMS is GNSS-referenced. So a large median difference is
# EXPECTED and is simply the datum offset. It is reported and then REMOVED.
#
# What actually discriminates (all invariant to a constant datum offset):
#   * spread across LAND-COVER classes -> the two datasets measure different
#                                         physical surfaces (real signal)
#   * spatial trend / smooth gradient  -> residual ramp, orbit or APS related
#   * large MAD, no structure          -> noise
# =============================================================================

from qgis.core import (QgsProject, QgsSpatialIndex, QgsRectangle, QgsPointXY,
                       QgsVectorLayer, QgsField, QgsFeature, QgsGeometry)
from qgis.PyQt.QtCore import QVariant
import statistics as st
import math
import os

# ----------------------------- CONFIG ----------------------------------------
CASES = {
    "Kiruna":      ("egms_kir_coh70_SW99TM",          "PSDS_kir_coh70_SW99TM"),
    "Lilla Edet":  ("egms_LE_coh70_SW99TM",           "PSDS_LE_coh70_SW99TM"),
    "Lodose":      ("egms_LOA_coh70_SW99TM",          "PSDS_LOA_coh70_SW99TM"),
    "Gothenburg":  ("egms_gbg_coh60_SW99TM_cor50",    "PSDS_gbg_coh60_SW99TM_cor50"),
    "Stockholm":   ("egms_sthlm_coh60_SW99TM_cor50",  "PSDS_sthlm_coh60_SW99TM_cor50"),
}

# Velocity field names differ between layers, so give candidates in priority
# order - the first one present in a given layer is used and reported.
EGMS_VEL_FIELDS = ["mean_velocity", "mean_veloc"]
PSDS_VEL_FIELDS = ["velocity", "mean_velocity"]
NMD_LAYER      = "NMD2023bas_v0_2"    # set None to skip land-cover analysis

RADII      = [5, 10, 15, 20.0, 50.0]   # co-location radii, metres (SWEREF99TM)
MIN_NEIGHB = 3              # min PS+DS points near an EGMS point
MIN_PER_LC = 8              # min points to report a land-cover class
MAKE_LAYERS = True
MAKE_PLOTS  = True          # scatter plots (needs matplotlib)
EXPORT_CSV  = True          # write per-corridor paired velocities
OUT_DIR     = "/Users/saam7695/Desktop/insar_regression"   # <-- set your folder
PLOT_RADIUS = 10.0          # which radius to plot / export
# -----------------------------------------------------------------------------


def get_layer(name):
    l = QgsProject.instance().mapLayersByName(name)
    return l[0] if l else None


def landcover_labels(raster):
    """Pull class codes -> labels from a paletted/categorised raster renderer."""
    labels = {}
    if raster is None:
        return labels
    try:
        rend = raster.renderer()
        if hasattr(rend, "classes"):
            for c in rend.classes():
                labels[int(c.value)] = c.label or ("class %d" % int(c.value))
    except Exception as e:
        print("  (could not read land-cover labels: %s)" % e)
    return labels


def resolve_field(layer, candidates):
    """Return the first candidate field present in the layer."""
    names = [f.name() for f in layer.fields()]
    lower = {n.lower(): n for n in names}
    for c in candidates:
        if c in names:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    raise ValueError("none of %s found in '%s'; available: %s"
                     % (candidates, layer.name(), names))


def read_points(layer, vel_field):
    pts = []
    for feat in layer.getFeatures():
        g = feat.geometry()
        if g is None or g.isEmpty():
            continue
        v = feat[vel_field]
        if v is None:
            continue
        try:
            v = float(v)
        except (TypeError, ValueError):
            continue
        p = g.asPoint()
        pts.append((p.x(), p.y(), v))
    return pts


def detect_unit_scale(ev, pv):
    """SARvey may write m/yr while EGMS is mm/yr."""
    if not ev or not pv:
        return 1.0, ""
    e = st.median([abs(v) for v in ev]) or 1e-9
    p = st.median([abs(v) for v in pv]) or 1e-9
    if e / p > 100:
        return 1000.0, "   [PS+DS appears m/yr -> x1000]"
    if p / e > 100:
        return 0.001, "   [PS+DS appears um/yr -> x0.001]"
    return 1.0, ""


def mad(vals, med):
    return 1.4826 * st.median([abs(v - med) for v in vals]) if len(vals) > 1 else 0.0


def pearson(a, b):
    n = len(a)
    if n < 3:
        return float("nan")
    ma, mb = st.mean(a), st.mean(b)
    num = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    da = sum((v - ma) ** 2 for v in a) ** 0.5
    db = sum((v - mb) ** 2 for v in b) ** 0.5
    return num / (da * db) if da > 0 and db > 0 else float("nan")


def ols_fit(x, y):
    """Ordinary least squares: y = slope*x + intercept."""
    n = len(x); mx, my = st.mean(x), st.mean(y)
    sxy = sum((x[i]-mx)*(y[i]-my) for i in range(n))
    sxx = sum((v-mx)**2 for v in x)
    if sxx == 0:
        return float("nan"), float("nan")
    b = sxy/sxx
    return b, my - b*mx


def orthogonal_fit(x, y):
    """Total least squares (Deming with equal error variances).
    Corrects the slope attenuation that OLS suffers when BOTH axes carry
    comparable measurement error - which is the case here."""
    n = len(x); mx, my = st.mean(x), st.mean(y)
    sxx = sum((v-mx)**2 for v in x)/(n-1)
    syy = sum((v-my)**2 for v in y)/(n-1)
    sxy = sum((x[i]-mx)*(y[i]-my) for i in range(n))/(n-1)
    if sxy == 0:
        return float("nan"), float("nan")
    b = (syy - sxx + math.sqrt((syy-sxx)**2 + 4*sxy**2)) / (2*sxy)
    return b, my - b*mx


def rankdata(a):
    idx = sorted(range(len(a)), key=lambda i: a[i])
    r = [0.0]*len(a); i = 0
    while i < len(idx):
        j = i
        while j+1 < len(idx) and a[idx[j+1]] == a[idx[i]]:
            j += 1
        avg = (i+j)/2.0 + 1
        for k in range(i, j+1):
            r[idx[k]] = avg
        i = j+1
    return r


def spearman(x, y):
    return pearson(rankdata(x), rankdata(y))


def percentile(sorted_vals, q):
    """Linear-interpolated percentile of an already-sorted list."""
    if not sorted_vals:
        return float("nan")
    k = (len(sorted_vals) - 1) * q
    f = int(math.floor(k)); c = min(f + 1, len(sorted_vals) - 1)
    return sorted_vals[f] + (sorted_vals[c] - sorted_vals[f]) * (k - f)


def theil_sen(x, y, max_pairs=200000, seed=0):
    """Robust (outlier-resistant) slope: median of all pairwise slopes.
    For large n a random subsample of pairs is used."""
    import random as _random
    n = len(x)
    slopes = []
    if n * (n - 1) // 2 <= max_pairs:
        for i in range(n):
            xi, yi = x[i], y[i]
            for j in range(i + 1, n):
                dx = x[j] - xi
                if abs(dx) > 1e-9:
                    slopes.append((y[j] - yi) / dx)
    else:
        rs = _random.Random(seed)
        for _ in range(max_pairs):
            i = rs.randrange(n); j = rs.randrange(n)
            if i == j:
                continue
            dx = x[j] - x[i]
            if abs(dx) > 1e-9:
                slopes.append((y[j] - y[i]) / dx)
    if not slopes:
        return float("nan"), float("nan")
    b = st.median(slopes)
    return b, st.median([y[i] - b * x[i] for i in range(n)])


def make_outlier_layer(name, recs, crs):
    """recs: (x, y, egms_v, psds_v, diff, n_mad)"""
    lyr = QgsVectorLayer("Point?crs=%s" % crs, name, "memory")
    dp = lyr.dataProvider()
    dp.addAttributes([QgsField("egms_v", QVariant.Double),
                      QgsField("psds_v", QVariant.Double),
                      QgsField("diff", QVariant.Double),
                      QgsField("n_mad", QVariant.Double)])
    lyr.updateFields()
    fs = []
    for x, y, ev, pv, d, nm in recs:
        f = QgsFeature(lyr.fields())
        f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
        f["egms_v"] = float(ev); f["psds_v"] = float(pv)
        f["diff"] = float(d); f["n_mad"] = float(nm)
        fs.append(f)
    dp.addFeatures(fs)
    lyr.updateExtents()
    QgsProject.instance().addMapLayer(lyr)
    return len(fs)


def sample_lc(raster, x, y):
    if raster is None:
        return None
    val, ok = raster.dataProvider().sample(QgsPointXY(x, y), 1)
    if not ok or val is None:
        return None
    try:
        return int(val)
    except (TypeError, ValueError):
        return None


def make_diff_layer(name, recs, crs):
    """recs: (x, y, raw_diff, residual_diff, lc)"""
    lyr = QgsVectorLayer("Point?crs=%s" % crs, name, "memory")
    dp = lyr.dataProvider()
    dp.addAttributes([QgsField("diff_raw", QVariant.Double),
                      QgsField("diff_resid", QVariant.Double),
                      QgsField("landcover", QVariant.Int)])
    lyr.updateFields()
    fs = []
    for x, y, draw, dres, lc in recs:
        f = QgsFeature(lyr.fields())
        f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
        f["diff_raw"] = float(draw)
        f["diff_resid"] = float(dres)
        f["landcover"] = int(lc) if lc is not None else -1
        fs.append(f)
    dp.addFeatures(fs)
    lyr.updateExtents()
    QgsProject.instance().addMapLayer(lyr)


# ================================ MAIN =======================================
nmd = get_layer(NMD_LAYER) if NMD_LAYER else None
if NMD_LAYER and nmd is None:
    print("WARNING: land-cover layer '%s' not found - skipping that section.\n" % NMD_LAYER)
LC = landcover_labels(nmd)
if LC:
    print("Read %d land-cover class labels from raster.\n" % len(LC))

print("=" * 78)
print("PART 1 - RAW DATUM OFFSET  (layers are unanchored: offset is EXPECTED)")
print("=" * 78)
print("%-13s %10s %10s %10s %8s %8s" %
      ("case", "EGMS med", "PSDS med", "offset", "nE", "nP"))
print("-" * 78)

store = {}
for case, (en, pn) in CASES.items():
    el, pl = get_layer(en), get_layer(pn)
    if el is None or pl is None:
        print("%-13s MISSING: %s" % (case, en if el is None else pn))
        continue
    try:
        ef = resolve_field(el, EGMS_VEL_FIELDS)
        pf = resolve_field(pl, PSDS_VEL_FIELDS)
        ep = read_points(el, ef)
        pp = read_points(pl, pf)
    except ValueError as err:
        print("%-13s %s" % (case, err))
        continue
    fieldnote = "" if (ef, pf) == (EGMS_VEL_FIELDS[0], PSDS_VEL_FIELDS[0]) \
                else "   [fields: %s / %s]" % (ef, pf)
    if not ep or not pp:
        print("%-13s no valid points" % case)
        continue
    sc, note = detect_unit_scale([v for _, _, v in ep], [v for _, _, v in pp])
    pp = [(x, y, v * sc) for x, y, v in pp]
    em, pm = st.median([v for _, _, v in ep]), st.median([v for _, _, v in pp])
    print("%-13s %10.2f %10.2f %10.2f %8d %8d%s%s"
          % (case, em, pm, pm - em, len(ep), len(pp), note, fieldnote))
    store[case] = (ep, pp, el.crs().authid())

print("\nThis offset mixes the arbitrary PS+DS reference with the GNSS-tied EGMS")
print("datum. It is NOT evidence of a processing difference. Part 2 removes it.")

print("\n" + "=" * 78)
print("PART 2 - STRUCTURE OF THE RESIDUAL  (after removing the datum offset)")
print("=" * 78)
print("datum-only  -> residual MAD small, land-cover spread small, no xy trend")
print("surfaces    -> LARGE spread across land-cover classes")
print("ramp/APS    -> |corr(x)| or |corr(y)| high")
print("noise       -> large MAD, no structure\n")

PAIRS = {}
PAIRS_XY = {}
CRS_OF = {}
for case, (ep, pp, crs) in store.items():
    print("-" * 78)
    print("CASE: %s" % case)

    idx = QgsSpatialIndex()
    look = {}
    for i, (x, y, v) in enumerate(pp):
        f = QgsFeature(i)
        f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
        idx.addFeature(f)
        look[i] = (x, y, v)

    for r in RADII:
        recs = []
        r2 = r * r
        for (ex, ey, ev) in ep:
            near = []
            for cid in idx.intersects(QgsRectangle(ex - r, ey - r, ex + r, ey + r)):
                px, py, pv = look[cid]
                if (px - ex) ** 2 + (py - ey) ** 2 <= r2:
                    near.append(pv)
            if len(near) < MIN_NEIGHB:
                continue
            # store: x, y, EGMS velocity, median PS+DS velocity, land cover
            recs.append([ex, ey, ev, st.median(near), sample_lc(nmd, ex, ey)])

        if len(recs) < 10:
            print("  r=%4.0f m : only %d pairs - too few" % (r, len(recs)))
            continue

        raw = [pv - ev for _, _, ev, pv, _ in recs]
        offset = st.median(raw)                      # datum offset estimate
        resid = [d - offset for d in raw]            # structure of interest
        spread = mad(resid, 0.0)
        xs = [x for x, _, _, _, _ in recs]
        ys = [y for _, y, _, _, _ in recs]
        rx, ry = pearson(xs, resid), pearson(ys, resid)

        print("  r=%4.0f m : n=%5d | datum offset = %+6.2f mm/yr"
              % (r, len(recs), offset))
        print("             residual MAD = %5.2f mm/yr | corr(x)=%+.2f corr(y)=%+.2f"
              % (spread, rx, ry))

        if nmd is not None:
            by = {}
            for _, _, ev, pv, lc in recs:
                if lc is not None:
                    by.setdefault(lc, []).append((pv - ev) - offset)
            sel = [(lc, st.median(v), len(v)) for lc, v in by.items()
                   if len(v) >= MIN_PER_LC]
            if sel:
                sel.sort(key=lambda t: -t[2])
                print("             residual by land cover (>=%d pts):" % MIN_PER_LC)
                for lc, m, n in sel[:9]:
                    print("               %-26s %+6.2f  (n=%d)"
                          % (LC.get(lc, "class %d" % lc)[:26], m, n))
                sp = [m for _, m, _ in sel]
                print("             >>> LAND-COVER SPREAD = %.2f mm/yr  %s"
                      % (max(sp) - min(sp),
                         "<-- STRONG, likely different surfaces"
                         if (max(sp) - min(sp)) > 1.5 else ""))

        if abs(r - PLOT_RADIUS) < 1e-6:
            PAIRS[case] = [(ev, pv) for _, _, ev, pv, _ in recs]
            PAIRS_XY[case] = [(x, y, ev, pv) for x, y, ev, pv, _ in recs]
            CRS_OF[case] = crs

        if MAKE_LAYERS and r == RADII[-1]:
            out = [(x, y, pv - ev, (pv - ev) - offset, lc)
                   for x, y, ev, pv, lc in recs]
            make_diff_layer("diff_%s_r%d" % (case.replace(" ", "_"), int(r)),
                            out, crs)


# ============================================================================
# PART 3 - REGRESSION OF PS+DS AGAINST EGMS  (reviewer comment 3)
# ============================================================================
# Reported per corridor:
#   OLS slope        - conventional, but ATTENUATED when both axes carry error
#   Orthogonal slope - total least squares, corrects that attenuation
#   intercept        - should reproduce the datum offset from Part 2
#   Pearson / Spearman r, R^2, n, and the EGMS dynamic range
#
# NOTE ON INTERPRETATION: r depends on how much real signal there is to
# correlate. A corridor with a small velocity range and ~0.6 mm/yr noise will
# show a lower r than one with an 8 mm/yr subsidence bowl, even if the two
# products agree equally well. ALWAYS report r together with the range.
# ============================================================================

print("\n" + "=" * 78)
print("PART 3 - REGRESSION:  PS+DS  =  slope x EGMS  +  intercept   (r=%.0f m)"
      % PLOT_RADIUS)
print("=" * 78)
print("%-12s %6s %8s %8s %9s %7s %8s %14s" %
      ("case", "n", "OLS b", "orth b", "intcpt", "r", "rho", "EGMS range"))
print("-" * 78)

summary = {}
for case, pr in PAIRS.items():
    if len(pr) < 10:
        print("%-12s too few pairs" % case)
        continue
    x = [a for a, _ in pr]      # EGMS velocity
    y = [b for _, b in pr]      # PS+DS velocity
    b_ols, a_ols = ols_fit(x, y)
    b_orth, a_orth = orthogonal_fit(x, y)
    r = pearson(x, y)
    rho = spearman(x, y)
    rng = max(x) - min(x)
    summary[case] = dict(n=len(pr), b_ols=b_ols, a_ols=a_ols, b_orth=b_orth,
                         a_orth=a_orth, r=r, rho=rho, rng=rng,
                         xmin=min(x), xmax=max(x))
    print("%-12s %6d %8.3f %8.3f %9.2f %7.3f %8.3f %6.1f to %5.1f" %
          (case, len(pr), b_ols, b_orth, a_ols, r, rho, min(x), max(x)))

print("\nR^2 (= r^2):  " + "   ".join("%s %.3f" % (c, v["r"]**2)
                                       for c, v in summary.items()))
print("\nA slope near 1 means the two products report the same deformation")
print("MAGNITUDE; the intercept is the datum offset. Compare OLS vs orthogonal:")
print("a much lower OLS slope indicates regression dilution, not disagreement.")

# ---------------------------- export + plots --------------------------------
if EXPORT_CSV or MAKE_PLOTS:
    try:
        os.makedirs(OUT_DIR, exist_ok=True)
    except Exception as e:
        print("\nCould not create OUT_DIR (%s): %s" % (OUT_DIR, e))

if EXPORT_CSV:
    for case, pr in PAIRS.items():
        fp = os.path.join(OUT_DIR, "pairs_%s.csv" % case.replace(" ", "_"))
        try:
            with open(fp, "w") as fh:
                fh.write("egms_vel_mmyr,psds_vel_mmyr\n")
                for a, b in pr:
                    fh.write("%.4f,%.4f\n" % (a, b))
        except Exception as e:
            print("CSV write failed for %s: %s" % (case, e))
    print("\nCSVs written to %s" % OUT_DIR)

if MAKE_PLOTS:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print("matplotlib unavailable (%s) - skipping plots" % e)
    else:
        n_case = len(PAIRS)
        ncol = 3
        nrow = int(math.ceil(n_case / float(ncol)))
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.0*ncol, 3.8*nrow))
        axes = axes.ravel() if n_case > 1 else [axes]

        for ax, (case, pr) in zip(axes, PAIRS.items()):
            if case not in summary:
                ax.axis("off"); continue
            s_ = summary[case]
            x = [a for a, _ in pr]; y = [b for _, b in pr]
            ax.scatter(x, y, s=5, alpha=0.25, edgecolors="none", color="#1f77b4")

            lo = min(min(x), min(y)); hi = max(max(x), max(y))
            pad = 0.05*(hi-lo) if hi > lo else 1.0
            lo -= pad; hi += pad

            ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="1:1")
            xx = [lo, hi]
            ax.plot(xx, [s_["b_ols"]*v + s_["a_ols"] for v in xx],
                    "-", color="#d62728", lw=1.4, label="OLS")
            ax.plot(xx, [s_["b_orth"]*v + s_["a_orth"] for v in xx],
                    "-", color="#2ca02c", lw=1.4, label="orthogonal")

            ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
            ax.set_title(case, fontsize=11)
            ax.set_xlabel("EGMS-PS velocity (mm/yr)", fontsize=9)
            ax.set_ylabel("PS+DS velocity (mm/yr)", fontsize=9)
            ax.tick_params(labelsize=8)
            ax.text(0.03, 0.97,
                    "n=%d\nb(orth)=%.2f\nintercept=%.2f\nr=%.2f" %
                    (s_["n"], s_["b_orth"], s_["a_orth"], s_["r"]),
                    transform=ax.transAxes, va="top", ha="left", fontsize=8,
                    bbox=dict(fc="white", ec="0.7", alpha=0.85, pad=3))
            ax.legend(fontsize=7, loc="lower right", framealpha=0.85)
            ax.grid(alpha=0.25, lw=0.5)

        for ax in axes[len(PAIRS):]:
            ax.axis("off")

        fig.tight_layout()
        out_png = os.path.join(OUT_DIR, "regression_egms_vs_psds.png")
        out_pdf = os.path.join(OUT_DIR, "regression_egms_vs_psds.pdf")
        try:
            fig.savefig(out_png, dpi=300)
            fig.savefig(out_pdf)
            print("Plots written:\n  %s\n  %s" % (out_png, out_pdf))
        except Exception as e:
            print("Plot save failed: %s" % e)
        plt.close(fig)


# ============================================================================
# PART 4 - ROBUST REGRESSION AND OUTLIER DIAGNOSTIC
# ============================================================================
# Part 3 uses moment statistics (covariance, Pearson r) which are dominated by
# extreme values, whereas Part 2 uses robust statistics describing the bulk.
# When they disagree, it is because (a) a small tail of outliers inflates the
# variance of one axis without contributing covariance, and/or (b) the corridor
# contains too little real deformation for a correlation to be meaningful.
#
# This section separates those effects:
#   SNR  = IQR of EGMS velocities / residual MAD.  This compares the real
#          deformation signal against the noise floor. SNR below ~2 means the
#          corridor is essentially stable and NO correlation can be expected,
#          however well the two products agree.
#   Theil-Sen and 3xMAD-trimmed fits are outlier-resistant and recover the
#          true slope and the datum offset where the raw fit fails.
# ============================================================================

TRIM_NMAD  = 3.0      # trim pairs deviating more than this many MAD
EXPORT_OUTLIERS = True

print("\n" + "=" * 78)
print("PART 4 - ROBUST REGRESSION AND OUTLIERS  (r=%.0f m)" % PLOT_RADIUS)
print("=" * 78)

robust = {}
for case, pr in PAIRS.items():
    if len(pr) < 20:
        print("%-12s too few pairs" % case); continue
    x = [a for a, _ in pr]; y = [b for _, b in pr]
    d = [y[i] - x[i] for i in range(len(x))]
    med = st.median(d)
    spread = mad(d, med)
    ds = sorted(d)
    xs_sorted = sorted(x)
    iqr_x = percentile(xs_sorted, 0.75) - percentile(xs_sorted, 0.25)
    snr = iqr_x / spread if spread > 0 else float("nan")

    keep = [i for i in range(len(x)) if abs(d[i] - med) <= TRIM_NMAD * spread]
    xt = [x[i] for i in keep]; yt = [y[i] for i in keep]
    n_out = len(x) - len(keep)

    b_ts, a_ts = theil_sen(x, y)
    if len(xt) >= 10:
        b_ot, a_ot = ols_fit(xt, yt)
        b_tt, a_tt = orthogonal_fit(xt, yt)
        r_t = pearson(xt, yt)
    else:
        b_ot = a_ot = b_tt = a_tt = r_t = float("nan")

    robust[case] = dict(med=med, mad=spread, iqr=iqr_x, snr=snr,
                        n=len(x), n_out=n_out, b_ts=b_ts, a_ts=a_ts,
                        b_ot=b_ot, a_ot=a_ot, b_tt=b_tt, a_tt=a_tt, r_t=r_t,
                        p01=percentile(ds, 0.01), p99=percentile(ds, 0.99))

    print("-" * 78)
    print("CASE: %s   (n=%d)" % (case, len(x)))
    print("  difference percentiles (mm/yr):  1%%=%+.2f  50%%=%+.2f  99%%=%+.2f"
          % (percentile(ds, 0.01), med, percentile(ds, 0.99)))
    print("  residual MAD = %.2f | EGMS IQR = %.2f | SNR = IQR/MAD = %.1f  %s"
          % (spread, iqr_x, snr,
             "<-- TOO LOW: correlation not meaningful" if snr < 2.0 else ""))
    print("  outliers beyond %.0f x MAD: %d of %d (%.1f%%)"
          % (TRIM_NMAD, n_out, len(x), 100.0 * n_out / len(x)))
    print("  Theil-Sen (robust)      slope %6.3f   intercept %+7.2f" % (b_ts, a_ts))
    print("  trimmed OLS             slope %6.3f   intercept %+7.2f" % (b_ot, a_ot))
    print("  trimmed orthogonal      slope %6.3f   intercept %+7.2f   r = %.3f"
          % (b_tt, a_tt, r_t))
    print("  (for reference, Part 2 datum offset was %+.2f)" % med)

    if EXPORT_OUTLIERS and n_out > 0 and case in PAIRS_XY:
        orecs = []
        for (px, py, ev, pv) in PAIRS_XY[case]:
            dd = pv - ev
            if spread > 0 and abs(dd - med) > TRIM_NMAD * spread:
                orecs.append((px, py, ev, pv, dd, (dd - med) / spread))
        if orecs:
            make_outlier_layer("outliers_%s" % case.replace(" ", "_"),
                               orecs, CRS_OF.get(case, "EPSG:3006"))

print("\n" + "=" * 78)
print("SUMMARY - which corridors support a meaningful regression")
print("=" * 78)
print("%-12s %6s %6s %7s %9s %10s %9s" %
      ("case", "SNR", "out%", "TheilSen", "trim orth", "trim r", "offset"))
print("-" * 78)
for case, v in robust.items():
    verdict = "yes" if v["snr"] >= 2.0 else "NO (flat)"
    print("%-12s %6.1f %5.1f%% %8.3f %9.3f %9.3f %9.2f   %s" %
          (case, v["snr"], 100.0*v["n_out"]/v["n"], v["b_ts"], v["b_tt"],
           v["r_t"], v["med"], verdict))
print("\nReport the TRIMMED / Theil-Sen slopes, not the raw Part 3 values.")
print("For corridors flagged NO, state that the deformation range is comparable")
print("to the noise floor so correlation is uninformative; agreement there is")
print("demonstrated by the offset and residual statistics instead.")
print("Outlier layers added to the project - inspect them for unwrapping errors")
print("or genuinely fast-moving targets (relevant to reviewer comment 1).")

print("\n" + "=" * 78)
print("Layers added: style by 'diff_resid' (diverging ramp, centred on 0).")
print("  uniform colour        -> datum only, nothing further to explain")
print("  follows land use      -> different physical surfaces (real signal)")
print("  smooth gradient       -> residual ramp / atmosphere")
print("  salt-and-pepper       -> noise")
print("=" * 78)