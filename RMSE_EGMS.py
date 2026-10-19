# =============================================================================
# RMSE of displacement time series with the EGMS definition       (QGIS console)
#
# EGMS Product Description: rmse is "evaluated on the time series residuals
# after applying a regression model of a third order polynomial plus a
# seasonal (sinusoidal) component";  RMSE = sqrt(mean((TS - model)^2)), mm.
#
# Applies exactly that model to BOTH datasets:
#   - PS+DS : new field  rmse_c3s  (comparable with EGMS rmse)
#   - EGMS  : new field  rmse_c3s  + check against the official 'rmse' field
#             (if they agree, the definition is reproduced correctly)
#
# Date fields:  EGMS "20200608",  PS+DS "D20200608"
# =============================================================================

from qgis.core import QgsProject, QgsField
from qgis.PyQt.QtCore import QVariant
import numpy as np
import datetime as dt
import re

# ----------------------------- CONFIG ----------------------------------------
CASES = {
    "Kiruna":      ("egms_kir_coh70_SW99TM",          "PSDS_kir_coh70_SW99TM"),
    "Lilla Edet":  ("egms_LE_coh70_SW99TM",           "PSDS_LE_coh70_SW99TM"),
    "Lodose":      ("egms_LOA_coh70_SW99TM",          "PSDS_LOA_coh70_SW99TM"),
    "Gothenburg":  ("egms_gbg_coh60_SW99TM_cor50",    "PSDS_gbg_coh60_SW99TM_cor50"),
    "Stockholm":   ("egms_sthlm_coh60_SW99TM_cor50",  "PSDS_sthlm_coh60_SW99TM_cor50"),
}
EGMS_DATE_RE = re.compile(r"^(\d{8})$")
PSDS_DATE_RE = re.compile(r"^[Dd](\d{8})$")
EGMS_RMSE_FIELDS = ["rmse", "RMSE"]
NEW_FIELD = "rmse_c3s"
WRITE_FIELD = True          # False: only print statistics
THRESHOLDS = (3.0, 4.0)     # report share of points at or below these (mm)
# -----------------------------------------------------------------------------


# ============================ pure helpers ===================================

def decimal_year(d):
    start = dt.date(d.year, 1, 1)
    n = (dt.date(d.year + 1, 1, 1) - start).days
    return d.year + (d - start).days / float(n)


def design_matrix(t):
    """Cubic polynomial + annual sinusoid, t in years."""
    t = np.asarray(t, float) - np.mean(t)
    w = 2 * np.pi * t
    return np.column_stack([np.ones_like(t), t, t**2, t**3, np.sin(w), np.cos(w)])


def rmse_egms(Y, t):
    """Y: (n_points, n_dates) displacement in mm (NaN allowed).
    Returns RMSE (mm) of residuals about the EGMS model, per point."""
    A = design_matrix(t)
    out = np.full(Y.shape[0], np.nan)
    full = ~np.isnan(Y).any(axis=1)
    if full.any():                                  # fast path: same dates
        coef, *_ = np.linalg.lstsq(A, Y[full].T, rcond=None)
        res = Y[full].T - A @ coef
        out[full] = np.sqrt(np.mean(res**2, axis=0))
    for i in np.where(~full)[0]:                    # rows with gaps
        m = ~np.isnan(Y[i])
        if m.sum() > A.shape[1] + 2:
            c, *_ = np.linalg.lstsq(A[m], Y[i, m], rcond=None)
            out[i] = np.sqrt(np.mean((Y[i, m] - A[m] @ c) ** 2))
    return out


def to_mm(Y):
    """SARvey may export metres. Detect from the typical range of the series."""
    rng = np.nanmedian(np.nanmax(Y, axis=1) - np.nanmin(Y, axis=1))
    if rng < 0.5:
        return Y * 1000.0, "metres -> mm"
    return Y, "mm"


# ============================ QGIS part ======================================

def num(v):
    """QGIS NULL comes back as a QVariant: treat it, None and text as NaN."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def get_layer(name):
    ls = QgsProject.instance().mapLayersByName(name)
    if not ls:
        raise RuntimeError("Layer not found: %s" % name)
    return ls[0]


def read_series(layer, date_re):
    fields = []
    for f in layer.fields():
        m = date_re.match(f.name())
        if m:
            fields.append((f.name(), dt.datetime.strptime(m.group(1), "%Y%m%d").date()))
    fields.sort(key=lambda x: x[1])
    names = [n for n, _ in fields]
    t = np.array([decimal_year(d) for _, d in fields])
    ids, rows = [], []
    for feat in layer.getFeatures():
        ids.append(feat.id())
        rows.append([num(feat[n]) for n in names])
    return ids, np.array(rows, float), t


def write_field(layer, ids, values):
    pr = layer.dataProvider()
    if layer.fields().indexOf(NEW_FIELD) < 0:
        pr.addAttributes([QgsField(NEW_FIELD, QVariant.Double)])
        layer.updateFields()
    idx = layer.fields().indexOf(NEW_FIELD)
    changes = {fid: {idx: (None if np.isnan(v) else float(v))}
               for fid, v in zip(ids, values)}
    ok = pr.changeAttributeValues(changes)
    layer.triggerRepaint()
    return ok


def summary(v):
    v = v[~np.isnan(v)]
    s = "median %.2f | P80 %.2f | max %.2f" % (np.median(v), np.percentile(v, 80), v.max())
    for th in THRESHOLDS:
        s += " | <=%.0f mm: %.0f%%" % (th, 100 * np.mean(v <= th))
    return s


print("=" * 86)
print("RMSE with the EGMS model (cubic polynomial + annual sinusoid)")
print("=" * 86)
ONLY = None               # e.g. ["Stockholm"] to rerun one corridor
for case, (egms_name, psds_name) in CASES.items():
    if ONLY and case not in ONLY:
        continue
    print("-" * 86)
    print("CASE: %s" % case)
    for label, name, rx in (("EGMS ", egms_name, EGMS_DATE_RE),
                            ("PS+DS", psds_name, PSDS_DATE_RE)):
        lyr = get_layer(name)
        ids, Y, t = read_series(lyr, rx)
        Y, unit = to_mm(Y)
        r = rmse_egms(Y, t)
        print("  %s n=%d dates=%d (%s)  %s" % (label, len(ids), len(t), unit, summary(r)))
        if label == "EGMS ":
            fn = next((f for f in EGMS_RMSE_FIELDS if lyr.fields().indexOf(f) >= 0), None)
            if fn:
                off = np.array([num(f[fn]) for f in lyr.getFeatures()], float)
                ok = ~np.isnan(off) & ~np.isnan(r)
                d = r[ok] - off[ok]
                print("         check vs official '%s': median diff %+.2f mm, "
                      "median |diff| %.2f mm, r = %.3f"
                      % (fn, np.median(d), np.median(np.abs(d)),
                         np.corrcoef(r[ok], off[ok])[0, 1]))
                print("         official %s" % summary(off))
        if WRITE_FIELD:
            ok = write_field(lyr, ids, r)
            print("         field '%s' %s" % (NEW_FIELD, "written" if ok else "NOT written (layer read-only?)"))
print("\nDone.")