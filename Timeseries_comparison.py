# =============================================================================
# EGMS-PS  vs  PS+DS   displacement TIME-SERIES comparison       (QGIS console)
# Reviewer 2: "provide a comparison based on the displacement time series and
#              show a few representative time series from both"
#
# Method (same co-location as the velocity analysis):
#   1. Each EGMS point is paired with all PS+DS points within RADIUS; the PS+DS
#      series of the pair is the per-date MEDIAN of those neighbours.
#   2. Only acquisition dates present in BOTH datasets are used, and both series
#      are set to zero at the first common date.
#   3. DATUM: PS+DS is relative to an arbitrary scene reference, EGMS is
#      GNSS-tied. The difference d_i(t) = PSDS_i(t) - EGMS_i(t) therefore holds
#      a common-mode series c(t) (reference-point motion incl. GIA trend).
#      c(t) = median over all pairs of d_i(t) is estimated and REMOVED, the
#      time-series counterpart of the constant velocity offset.
#      Its linear trend should reproduce the velocity offset (consistency check).
#   4. Per pair: RMS of the residual difference, and Pearson r between the
#      LOCAL series (scene-wide median EGMS motion removed from both). Without
#      that removal, the regional uplift shared by every series inflates r.
#   5. Representative pairs are chosen by RULE, not by eye:
#        typical      - RMS closest to the corridor median
#        deforming    - largest |EGMS velocity| among pairs with RMS <= P75
#        disagreement - RMS at the 95th percentile
#
# Field formats:  EGMS dates "20200608", PS+DS dates "D20200608"
# =============================================================================

from qgis.core import (QgsProject, QgsSpatialIndex, QgsRectangle, QgsPointXY,
                       QgsVectorLayer, QgsField, QgsFeature, QgsGeometry)
from qgis.PyQt.QtCore import QVariant
import statistics as st
import datetime as dt
import math
import os
import re

# ----------------------------- CONFIG ----------------------------------------
CASES = {
    "Kiruna":      ("egms_kir_coh70_SW99TM",          "PSDS_kir_coh70_SW99TM"),
    "Lilla Edet":  ("egms_LE_coh70_SW99TM",           "PSDS_LE_coh70_SW99TM"),
    "Lodose":      ("egms_LOA_coh70_SW99TM",          "PSDS_LOA_coh70_SW99TM"),
    "Gothenburg":  ("egms_gbg_coh60_SW99TM_cor50",    "PSDS_gbg_coh60_SW99TM_cor50"),
    "Stockholm":   ("egms_sthlm_coh60_SW99TM_cor50",  "PSDS_sthlm_coh60_SW99TM_cor50"),
}

EGMS_DATE_RE = re.compile(r"^(\d{8})$")      # 20200608
PSDS_DATE_RE = re.compile(r"^[Dd](\d{8})$")  # D20200608

MODE          = "grid"   # "grid":  median EGMS vs median PS+DS in the same
                         #          non-overlapping square cell (symmetric,
                         #          every cell an independent comparison)
                         # "point": single EGMS point vs median PS+DS within RADIUS
CELL          = 30.0     # grid cell size, metres (MODE = "grid")
MIN_EGMS      = 1        # min EGMS points in a cell (MODE = "grid")
RADIUS        = 15.0     # co-location radius, metres (MODE = "point")
MIN_NEIGHB    = 3        # min PS+DS points per cell / near an EGMS point
MIN_DATES     = 10       # min common dates for a pair to be used
V_SIG         = 5.0      # mm/yr: local velocity above which a cell has real local signal
DATE_TOL_DAYS = 1        # allow +-1 day when matching dates (UTC midnight shifts)
MAKE_LAYERS   = True     # memory layer of per-pair RMS / r, to map
MAKE_PLOTS    = True
EXPORT_CSV    = True
OUT_DIR       = "/Users/saam7695/Desktop/insar_timeseries"   # <-- set your folder
# -----------------------------------------------------------------------------


# ============================ pure helpers ===================================
# (no QGIS in this block - kept separate so it can be tested on its own)

def parse_date(s):
    return dt.datetime.strptime(s, "%Y%m%d").date()


def decimal_year(d):
    start = dt.date(d.year, 1, 1)
    n = (dt.date(d.year + 1, 1, 1) - start).days
    return d.year + (d - start).days / float(n)


def match_dates(egms_dates, psds_dates, tol_days=1):
    """Return list of (egms_key, psds_key, date) for dates present in both.
    Exact matches first; otherwise nearest PS+DS date within tol_days."""
    ed = {k: parse_date(k) for k in egms_dates}
    pd = {k: parse_date(k) for k in psds_dates}
    by_date = {v: k for k, v in pd.items()}
    used, out = set(), []
    for ek, edate in sorted(ed.items(), key=lambda t: t[1]):
        pk = by_date.get(edate)
        if pk is None and tol_days > 0:
            best = None
            for k, pdate in pd.items():
                if k in used:
                    continue
                dd = abs((pdate - edate).days)
                if dd <= tol_days and (best is None or dd < best[0]):
                    best = (dd, k)
            pk = best[1] if best else None
        if pk is not None and pk not in used:
            used.add(pk)
            out.append((ek, pk, edate))
    return out


def unit_scale(egms_series, psds_series):
    """SARvey may store metres while EGMS stores mm. Compare typical ranges."""
    def typical_range(series_list):
        rs = []
        for s in series_list:
            v = [x for x in s if x is not None]
            if len(v) > 2:
                rs.append(max(v) - min(v))
        return st.median(rs) if rs else 0.0
    re_ = typical_range(egms_series)
    rp = typical_range(psds_series)
    if rp > 0 and re_ / rp > 100:
        return 1000.0, "[PS+DS appears to be in metres -> x1000]"
    if re_ > 0 and rp / re_ > 100:
        return 0.001, "[PS+DS appears to be in micrometres -> x0.001]"
    return 1.0, ""


def median_series(series_list):
    """Per-date median across several series (None ignored)."""
    n = len(series_list[0])
    out = []
    for t in range(n):
        v = [s[t] for s in series_list if s[t] is not None]
        out.append(st.median(v) if v else None)
    return out


def rereference(series):
    """Set series to zero at its first date; None if first value missing."""
    if series[0] is None:
        return None
    s0 = series[0]
    return [None if v is None else v - s0 for v in series]


def ols_slope(t, y):
    pts = [(a, b) for a, b in zip(t, y) if b is not None]
    if len(pts) < 3:
        return float("nan")
    mt = st.mean([a for a, _ in pts]); my = st.mean([b for _, b in pts])
    sxx = sum((a - mt) ** 2 for a, _ in pts)
    if sxx == 0:
        return float("nan")
    return sum((a - mt) * (b - my) for a, b in pts) / sxx


def pearson(a, b):
    pts = [(x, y) for x, y in zip(a, b) if x is not None and y is not None]
    if len(pts) < 3:
        return float("nan")
    ma = st.mean([x for x, _ in pts]); mb = st.mean([y for _, y in pts])
    num = sum((x - ma) * (y - mb) for x, y in pts)
    da = math.sqrt(sum((x - ma) ** 2 for x, _ in pts))
    db = math.sqrt(sum((y - mb) ** 2 for _, y in pts))
    return num / (da * db) if da > 0 and db > 0 else float("nan")


def percentile(vals, q):
    s = sorted(vals)
    if not s:
        return float("nan")
    k = (len(s) - 1) * q
    f = int(math.floor(k)); c = min(f + 1, len(s) - 1)
    return s[f] + (s[c] - s[f]) * (k - f)


def compare_pairs(pairs, tyears):
    """pairs: list of dict(egms=[...], psds=[...]) re-referenced, aligned.
    Returns common-mode series c(t) and per-pair stats (in place)."""
    n_t = len(tyears)
    diffs = []
    for p in pairs:
        diffs.append([None if (p["egms"][t] is None or p["psds"][t] is None)
                      else p["psds"][t] - p["egms"][t] for t in range(n_t)])
    c = []
    for t in range(n_t):
        v = [d[t] for d in diffs if d[t] is not None]
        c.append(st.median(v) if v else None)
    # scene-wide EGMS motion m(t) (regional uplift, seasonal, ...). It is shared
    # by every pair, so it inflates correlations; the LOCAL series below have it
    # removed and show only what differs from place to place.
    m = []
    for t in range(n_t):
        v = [p["egms"][t] for p in pairs if p["egms"][t] is not None]
        m.append(st.median(v) if v else None)
    for p, d in zip(pairs, diffs):
        resid = [None if (d[t] is None or c[t] is None) else d[t] - c[t]
                 for t in range(n_t)]
        # RMS about the pair's own mean, so the noise of the single reference
        # date (t0) does not inflate the result
        rv = [r for r in resid if r is not None]
        if len(rv) > 2:
            mr = st.mean(rv)
            p["rms"] = math.sqrt(sum((r - mr) ** 2 for r in rv) / len(rv))
        else:
            p["rms"] = float("nan")
        p["psds_corr"] = [None if (p["psds"][t] is None or c[t] is None)
                          else p["psds"][t] - c[t] for t in range(n_t)]
        p["r_total"] = pearson(p["egms"], p["psds_corr"])
        loc_e = [None if (p["egms"][t] is None or m[t] is None)
                 else p["egms"][t] - m[t] for t in range(n_t)]
        loc_p = [None if (p["psds_corr"][t] is None or m[t] is None)
                 else p["psds_corr"][t] - m[t] for t in range(n_t)]
        p["r"] = pearson(loc_e, loc_p)                 # LOCAL correlation
        p["v_egms"] = ols_slope(tyears, p["egms"])
        p["v_local"] = ols_slope(tyears, loc_e)        # velocity relative to scene
        e = [x for x in loc_e if x is not None]
        p["egms_span"] = (max(e) - min(e)) if e else float("nan")
    return c


def pick_representatives(pairs):
    """Return dict label -> pair, chosen by explicit rules."""
    ok = [p for p in pairs if not math.isnan(p["rms"])]
    if not ok:
        return {}
    rms = [p["rms"] for p in ok]
    med, p75, p95 = percentile(rms, 0.5), percentile(rms, 0.75), percentile(rms, 0.95)
    typical = min(ok, key=lambda p: abs(p["rms"] - med))
    good = [p for p in ok if p["rms"] <= p75 and not math.isnan(p["v_local"])]
    deforming = max(good, key=lambda p: abs(p["v_local"])) if good else typical
    disagree = min(ok, key=lambda p: abs(p["rms"] - p95))
    return {"typical (median RMS)": typical,
            "deforming (max local |v|, RMS <= P75)": deforming,
            "disagreement (RMS at P95)": disagree}


# ============================ QGIS reading ===================================

def get_layer(name):
    l = QgsProject.instance().mapLayersByName(name)
    return l[0] if l else None


def date_fields(layer, regex):
    """Map date string YYYYMMDD -> actual field name."""
    out = {}
    for f in layer.fields():
        m = regex.match(f.name())
        if m:
            out[m.group(1)] = f.name()
    return out


def read_series(layer, field_names):
    """Return list of (x, y, [values in field_names order])."""
    pts = []
    for feat in layer.getFeatures():
        g = feat.geometry()
        if g is None or g.isEmpty():
            continue
        vals = []
        for fn in field_names:
            v = feat[fn]
            try:
                vals.append(None if v is None else float(v))
            except (TypeError, ValueError):
                vals.append(None)
        p = g.asPoint()
        pts.append((p.x(), p.y(), vals))
    return pts


def make_example_layers(case, reps, crs):
    """Points (and, in grid mode, cell squares) for the plotted examples,
    with the same labels as the figure panels."""
    pt = QgsVectorLayer("Point?crs=%s" % crs, "ts_examples_%s" % case.replace(" ", "_"),
                        "memory")
    dp = pt.dataProvider()
    dp.addAttributes([QgsField("label", QVariant.String),
                      QgsField("kind", QVariant.String),
                      QgsField("rms_mm", QVariant.Double),
                      QgsField("r", QVariant.Double),
                      QgsField("n_egms", QVariant.Int),
                      QgsField("n_psds", QVariant.Int)])
    pt.updateFields()
    poly = None
    if MODE == "grid":
        poly = QgsVectorLayer("Polygon?crs=%s" % crs,
                              "ts_example_cells_%s" % case.replace(" ", "_"), "memory")
        pp_ = poly.dataProvider()
        pp_.addAttributes([QgsField("label", QVariant.String)])
        poly.updateFields()
    feats, cells = [], []
    for kind, p in reps.items():
        f = QgsFeature(pt.fields())
        f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(p["x"], p["y"])))
        f["label"] = p["label"]; f["kind"] = kind
        f["rms_mm"] = float(p["rms"]); f["r"] = None if math.isnan(p["r"]) else float(p["r"])
        f["n_egms"] = int(p["ne"]); f["n_psds"] = int(p["n"])
        feats.append(f)
        if poly is not None and p["cell"] is not None:
            i, j = p["cell"]
            g = QgsGeometry.fromRect(QgsRectangle(i * CELL, j * CELL,
                                                  (i + 1) * CELL, (j + 1) * CELL))
            c = QgsFeature(poly.fields()); c.setGeometry(g); c["label"] = p["label"]
            cells.append(c)
    dp.addFeatures(feats); pt.updateExtents()
    QgsProject.instance().addMapLayer(pt)
    if poly is not None:
        poly.dataProvider().addFeatures(cells); poly.updateExtents()
        QgsProject.instance().addMapLayer(poly)


def make_pair_layer(name, pairs, crs):
    lyr = QgsVectorLayer("Point?crs=%s" % crs, name, "memory")
    dp = lyr.dataProvider()
    dp.addAttributes([QgsField("rms_mm", QVariant.Double),
                      QgsField("r", QVariant.Double),
                      QgsField("v_egms", QVariant.Double),
                      QgsField("n_neigh", QVariant.Int)])
    lyr.updateFields()
    fs = []
    for p in pairs:
        f = QgsFeature(lyr.fields())
        f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(p["x"], p["y"])))
        f["rms_mm"] = float(p["rms"]) if not math.isnan(p["rms"]) else None
        f["r"] = float(p["r"]) if not math.isnan(p["r"]) else None
        f["v_egms"] = float(p["v_egms"]) if not math.isnan(p["v_egms"]) else None
        f["n_neigh"] = int(p["n"])
        fs.append(f)
    dp.addFeatures(fs)
    lyr.updateExtents()
    QgsProject.instance().addMapLayer(lyr)


# ================================ MAIN =======================================
if EXPORT_CSV or MAKE_PLOTS:
    try:
        os.makedirs(OUT_DIR, exist_ok=True)
    except Exception as e:
        print("Could not create OUT_DIR (%s): %s" % (OUT_DIR, e))

print("=" * 86)
if MODE == "grid":
    print("TIME-SERIES COMPARISON  EGMS-PS vs PS+DS   (grid cells %.0f m: median vs median)"
          % CELL)
else:
    print("TIME-SERIES COMPARISON  EGMS-PS vs PS+DS   (EGMS point vs PS+DS median, r=%.0f m)"
          % RADIUS)
print("=" * 86)

results = {}
for case, (en, pn) in CASES.items():
    el, pl = get_layer(en), get_layer(pn)
    if el is None or pl is None:
        print("%-12s MISSING LAYER: %s" % (case, en if el is None else pn)); continue

    efields = date_fields(el, EGMS_DATE_RE)
    pfields = date_fields(pl, PSDS_DATE_RE)
    if not efields or not pfields:
        print("%-12s no date fields found (EGMS %d, PS+DS %d)"
              % (case, len(efields), len(pfields))); continue

    matched = match_dates(list(efields), list(pfields), DATE_TOL_DAYS)
    if len(matched) < MIN_DATES:
        print("%-12s only %d common dates (EGMS %d, PS+DS %d) - skipped"
              % (case, len(matched), len(efields), len(pfields))); continue

    dates = [m[2] for m in matched]
    tyears = [decimal_year(d) for d in dates]
    e_names = [efields[m[0]] for m in matched]
    p_names = [pfields[m[1]] for m in matched]

    ep = read_series(el, e_names)
    pp = read_series(pl, p_names)
    sc, note = unit_scale([s for _, _, s in ep[:500]], [s for _, _, s in pp[:500]])
    if sc != 1.0:
        pp = [(x, y, [None if v is None else v * sc for v in s]) for x, y, s in pp]

    # ---- build comparison units -------------------------------------------
    # each unit: (x, y, list of EGMS series, list of PS+DS series, cell or None)
    units = []
    if MODE == "grid":
        cells = {}
        for (x, y, s) in ep:
            cells.setdefault((int(math.floor(x / CELL)), int(math.floor(y / CELL))),
                             ([], []))[0].append((x, y, s))
        for (x, y, s) in pp:
            k = (int(math.floor(x / CELL)), int(math.floor(y / CELL)))
            if k in cells:
                cells[k][1].append(s)
        for (i, j), (eg, ps) in cells.items():
            if len(eg) >= MIN_EGMS and len(ps) >= MIN_NEIGHB:
                cx, cy = (i + 0.5) * CELL, (j + 0.5) * CELL
                units.append((cx, cy, [s for _, _, s in eg], ps, (i, j)))
    else:
        idx = QgsSpatialIndex()
        for i, (x, y, _) in enumerate(pp):
            f = QgsFeature(i)
            f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
            idx.addFeature(f)
        r2 = RADIUS * RADIUS
        for (ex, ey, es) in ep:
            near = [pp[cid][2] for cid in idx.intersects(
                        QgsRectangle(ex - RADIUS, ey - RADIUS, ex + RADIUS, ey + RADIUS))
                    if (pp[cid][0] - ex) ** 2 + (pp[cid][1] - ey) ** 2 <= r2]
            if len(near) >= MIN_NEIGHB:
                units.append((ex, ey, [es], near, None))

    def band(series_list, q):
        out = []
        for t in range(len(dates)):
            v = [s[t] - s[0] for s in series_list if s[t] is not None and s[0] is not None]
            out.append(percentile(v, q) if v else None)
        return out

    pairs = []
    for (ux, uy, eg_list, ps_list, cell) in units:
        e_ref = rereference(median_series(eg_list))
        p_ref = rereference(median_series(ps_list))
        if e_ref is None or p_ref is None:
            continue
        if sum(1 for v in e_ref if v is not None) < MIN_DATES:
            continue
        pairs.append(dict(x=ux, y=uy, egms=e_ref, psds=p_ref, n=len(ps_list),
                          ne=len(eg_list), cell=cell,
                          psds_lo=band(ps_list, 0.25), psds_hi=band(ps_list, 0.75)))
    if len(pairs) < 10:
        print("%-12s only %d co-located pairs - skipped" % (case, len(pairs))); continue

    c = compare_pairs(pairs, tyears)
    c_trend = ols_slope(tyears, c)
    rms = [p["rms"] for p in pairs if not math.isnan(p["rms"])]
    rr = [p["r"] for p in pairs if not math.isnan(p["r"])]
    reps = pick_representatives(pairs)
    abbrev = "".join(w[0] for w in case.split()).upper()[:2] if " " in case else case[:3].upper()
    for k, (kind, p) in enumerate(reps.items()):
        p["label"] = "%s-%s" % (abbrev, "abc"[k])
    results[case] = dict(dates=dates, c=c, pairs=pairs, reps=reps,
                         crs=el.crs().authid())

    print("-" * 86)
    print("CASE: %s   %s" % (case, note))
    print("  common dates %d (%s to %s) | pairs %d"
          % (len(dates), dates[0], dates[-1], len(pairs)))
    print("  common-mode datum series: linear trend %+.2f mm/yr"
          "  (should match the velocity offset)" % c_trend)
    print("  RMS of residual difference (mm): median %.2f | P25 %.2f | P75 %.2f | P95 %.2f"
          % (st.median(rms), percentile(rms, .25), percentile(rms, .75), percentile(rms, .95)))
    rt = [p["r_total"] for p in pairs if not math.isnan(p["r_total"])]
    print("  Pearson r, LOCAL (scene motion removed): median %.2f | share r > 0.7: %.0f%%"
          % (st.median(rr), 100.0 * sum(1 for v in rr if v > 0.7) / len(rr)))
    print("  Pearson r, total (incl. shared regional motion): median %.2f  [inflated]"
          % st.median(rt))
    sig = [p for p in pairs if not math.isnan(p["r"]) and not math.isnan(p["v_local"])
           and abs(p["v_local"]) >= V_SIG]
    if sig:
        rs = [p["r"] for p in sig]
        print("  r, LOCAL, cells with |local v| >= %.0f mm/yr: median %.2f | share r > 0.7:"
              " %.0f%% | n = %d (%.0f%% of cells)"
              % (V_SIG, st.median(rs), 100.0 * sum(1 for v in rs if v > 0.7) / len(rs),
                 len(sig), 100.0 * len(sig) / len(pairs)))
    else:
        print("  r, LOCAL: no cells with |local v| >= %.0f mm/yr" % V_SIG)
    print("  examples (x, y in %s):" % el.crs().authid())
    for lab, p in reps.items():
        print("    %-6s %-38s RMS %.2f mm  r_loc %.2f  v_EGMS %+.1f (local %+.1f) mm/yr"
              "  x=%.0f y=%.0f  (EGMS %d, PS+DS %d pts)"
              % (p["label"], lab, p["rms"], p["r"], p["v_egms"], p["v_local"],
                 p["x"], p["y"], p["ne"], p["n"]))

    if MAKE_LAYERS:
        make_pair_layer("ts_pairs_%s" % case.replace(" ", "_"), pairs, el.crs().authid())
        make_example_layers(case, reps, el.crs().authid())

    if EXPORT_CSV:
        fp = os.path.join(OUT_DIR, "ts_pairs_%s.csv" % case.replace(" ", "_"))
        try:
            with open(fp, "w") as fh:
                fh.write("label,x,y,n_egms,n_psds,rms_mm,r_local,r_total,"
                         "v_egms_mmyr,v_local_mmyr,local_span_mm\n")
                for p in pairs:
                    fh.write("%s,%.2f,%.2f,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n"
                             % (p.get("label", ""), p["x"], p["y"], p["ne"], p["n"],
                                p["rms"], p["r"], p["r_total"], p["v_egms"],
                                p["v_local"], p["egms_span"]))
        except Exception as e:
            print("  CSV write failed: %s" % e)

# ------------------------------- summary --------------------------------------
print("\n" + "=" * 86)
print("SUMMARY (for the manuscript table)")
print("%-12s %6s %6s %11s %11s %9s %13s" %
      ("case", "pairs", "dates", "med RMS mm", "P95 RMS mm", "r_loc sig.", "c(t) trend"))
print("-" * 86)
for case, R in results.items():
    rms = [p["rms"] for p in R["pairs"] if not math.isnan(p["rms"])]
    rr = [p["r"] for p in R["pairs"] if not math.isnan(p["r"])
          and not math.isnan(p["v_local"]) and abs(p["v_local"]) >= V_SIG] or [float("nan")]
    print("%-12s %6d %6d %11.2f %11.2f %9.2f %+10.2f mm/yr" %
          (case, len(R["pairs"]), len(R["dates"]), st.median(rms),
           percentile(rms, .95), st.median(rr),
           ols_slope([decimal_year(d) for d in R["dates"]], R["c"])))

# ------------------------------- plots ----------------------------------------
if MAKE_PLOTS and results:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print("matplotlib unavailable (%s) - skipping plots" % e)
    else:
        nrow = len(results)
        fig, axes = plt.subplots(nrow, 3, figsize=(13, 2.6 * nrow), squeeze=False,
                                 sharex=False)
        for i, (case, R) in enumerate(results.items()):
            for j, (lab, p) in enumerate(R["reps"].items()):
                ax = axes[i][j]
                d = R["dates"]
                ev = [(t, v) for t, v in zip(d, p["egms"]) if v is not None]
                pv = [(t, v) for t, v in zip(d, p["psds_corr"]) if v is not None]
                band = [(t, lo - (cv or 0), hi - (cv or 0))
                        for t, lo, hi, cv in zip(d, p["psds_lo"], p["psds_hi"], R["c"])
                        if lo is not None and hi is not None and cv is not None]
                if band:
                    ax.fill_between([b[0] for b in band], [b[1] for b in band],
                                    [b[2] for b in band], color="#d95f02", alpha=0.15,
                                    lw=0, label="PS+DS neighbours (IQR)")
                def gapped(seq, max_gap=60):
                    """Insert a break so lines are not drawn across winter gaps."""
                    xs, ys = [], []
                    for k, (t, v) in enumerate(seq):
                        if k and (t - seq[k - 1][0]).days > max_gap:
                            xs.append(t); ys.append(float("nan"))
                        xs.append(t); ys.append(v)
                    return xs, ys
                ex_, ey_ = gapped(ev); px_, py_ = gapped(pv)
                ax.plot(ex_, ey_, "o-", ms=2.5, lw=0.8,
                        color="#1b4f8a", label="EGMS-PS")
                ax.plot(px_, py_, "s-", ms=2.5, lw=0.8,
                        color="#d95f02", label="PS+DS (median, datum-corrected)")
                ax.set_title("%s  %s - %s\nRMS %.1f mm, local r = %.2f  (x %.0f, y %.0f)"
                             % (p["label"], case, lab, p["rms"], p["r"], p["x"], p["y"]),
                             fontsize=8)
                ax.set_ylabel("LOS displacement (mm)", fontsize=8)
                ax.tick_params(labelsize=7)
                ax.grid(alpha=0.25, lw=0.5)
                if i == 0 and j == 0:
                    ax.legend(fontsize=6.5, loc="best", framealpha=0.85)
        fig.autofmt_xdate()
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fp = os.path.join(OUT_DIR, "timeseries_representative.%s" % ext)
            try:
                fig.savefig(fp, dpi=300 if ext == "png" else None)
                print("Plot written: %s" % fp)
            except Exception as e:
                print("Plot save failed: %s" % e)
        plt.close(fig)

print("\nLayers 'ts_pairs_<case>' added: style by rms_mm to map where the series agree.")