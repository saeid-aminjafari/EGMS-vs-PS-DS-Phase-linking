# =============================================================================
# NKG 3-D velocity validation of the EGMS-PS datum offset        (QGIS console)
#
# Projects the FULL NKG velocity vector (East, North, Up) into the satellite
# line of sight, instead of assuming purely vertical motion.
#
#   v_LOS = vE*sin(inc)*sin(head-90) + vN*sin(inc)*cos(head-90) + vU*cos(inc)
#
# where head is the satellite heading (flight azimuth, clockwise from north)
# and the sensor is right-looking. The script reports BOTH the vertical-only
# and the full 3-D projection so the contribution of the horizontal components
# is explicit.
#
# !! REFERENCE FRAME WARNING !!
# The horizontal components are only meaningful here if NKG_RF17vel is
# expressed in a Eurasia-fixed frame (ETRF/NKG_ETRF). If the grids are in
# ITRF, the horizontals contain ~25 mm/yr of Eurasian plate motion, which does
# NOT appear in EGMS (referenced to ETRS89) and would corrupt the comparison.
# The script prints the median horizontal rate: a few mm/yr means Eurasia-fixed
# and is fine; tens of mm/yr means ITRF - stop and convert the frame first.
# =============================================================================

from qgis.core import (QgsProject, QgsCoordinateTransform,
                       QgsCoordinateReferenceSystem, QgsPointXY)
import statistics as st
import math

# ----------------------------- CONFIG ----------------------------------------
GRI_U = "/Users/saam7695/Library/CloudStorage/OneDrive-Chalmers/FinalAll/nkg_rf17vel/NKG_RF17vel_U.GRI"   # up       component
GRI_E = "/Users/saam7695/Library/CloudStorage/OneDrive-Chalmers/FinalAll/nkg_rf17vel/NKG_RF17vel_E.GRI"   # east     component
GRI_N = "/Users/saam7695/Library/CloudStorage/OneDrive-Chalmers/FinalAll/nkg_rf17vel/NKG_RF17vel_N.GRI"   # north    component

# orbit per corridor: used only if no heading field is found in the layer
NOMINAL_HEADING = {"ascending": 347.0, "descending": 193.0}

CASES = {   # case: (layer name, orbit)
    "Kiruna":      ("egms_kir_coh70_SW99TM",         "ascending"),
    "Lilla Edet":  ("egms_LE_coh70_SW99TM",          "ascending"),
    "Lodose":      ("egms_LOA_coh70_SW99TM",         "ascending"),
    "Gothenburg":  ("egms_gbg_coh60_SW99TM_cor50",   "descending"),
    "Stockholm":   ("egms_sthlm_coh60_SW99TM_cor50", "descending"),
}

EGMS_VEL_FIELDS = ["mean_velocity", "mean_veloc"]
INC_FIELDS      = ["incidence_angle", "incidence_", "incidence"]
HEAD_FIELDS     = ["track_angle", "track_angl", "heading_angle", "heading",
                   "head_angle", "azimuth_angle", "satellite_heading"]
# -----------------------------------------------------------------------------


class GravsoftGrid:
    def __init__(self, path):
        tok = open(path).read().split()
        (self.lat_s, self.lat_n, self.lon_w,
         self.lon_e, self.dlat, self.dlon) = [float(v) for v in tok[:6]]
        vals = [float(v) for v in tok[6:]]
        self.nlat = int(round((self.lat_n - self.lat_s) / self.dlat)) + 1
        self.nlon = int(round((self.lon_e - self.lon_w) / self.dlon)) + 1
        if self.nlat * self.nlon != len(vals):
            raise ValueError("grid mismatch in %s: header implies %d, file has %d"
                             % (path, self.nlat * self.nlon, len(vals)))
        self.g = [vals[i * self.nlon:(i + 1) * self.nlon]
                  for i in range(self.nlat)]
        self.vmin, self.vmax = min(vals), max(vals)

    def sample(self, lat, lon):
        if not (self.lat_s <= lat <= self.lat_n and
                self.lon_w <= lon <= self.lon_e):
            return None
        fr = (self.lat_n - lat) / self.dlat
        fc = (lon - self.lon_w) / self.dlon
        r0, c0 = int(math.floor(fr)), int(math.floor(fc))
        r1, c1 = min(r0 + 1, self.nlat - 1), min(c0 + 1, self.nlon - 1)
        wr, wc = fr - r0, fc - c0
        return ((1 - wr) * ((1 - wc) * self.g[r0][c0] + wc * self.g[r0][c1]) +
                wr * ((1 - wc) * self.g[r1][c0] + wc * self.g[r1][c1]))


def los_unit(inc_deg, head_deg):
    """Unit vector (E,N,U) from ground target to satellite, right-looking SAR."""
    th = math.radians(inc_deg)
    phi = math.radians(head_deg - 90.0)
    return (math.sin(th) * math.sin(phi),
            math.sin(th) * math.cos(phi),
            math.cos(th))


def get_layer(name):
    l = QgsProject.instance().mapLayersByName(name)
    return l[0] if l else None


def resolve_field(layer, candidates):
    names = [f.name() for f in layer.fields()]
    lower = {n.lower(): n for n in names}
    for c in candidates:
        if c in names:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def mad(vals, med):
    return 1.4826 * st.median([abs(v - med) for v in vals]) if len(vals) > 1 else 0.0


def as_degrees(v):
    """For incidence angles: values below 2 must be radians."""
    return math.degrees(v) if abs(v) < 2.0 else v


def norm360(a):
    a = a % 360.0
    return a + 360.0 if a < 0 else a


def ang_diff(a, b):
    """Smallest absolute difference between two azimuths, in degrees."""
    d = abs(norm360(a) - norm360(b)) % 360.0
    return min(d, 360.0 - d)


def resolve_heading(raw, nominal):
    """A stored heading may be in degrees or radians, signed or 0-360.
    Both readings are tested and the one consistent with the known orbit
    direction is returned, together with the angular mismatch."""
    cands = [norm360(raw), norm360(math.degrees(raw))]
    best = min(cands, key=lambda c: ang_diff(c, nominal))
    return best, ang_diff(best, nominal)


# ================================ MAIN =======================================
gU = GravsoftGrid(GRI_U)
try:
    gE = GravsoftGrid(GRI_E)
    gN = GravsoftGrid(GRI_N)
    have3d = True
except Exception as e:
    print("Could not load E/N grids (%s) - running vertical-only.\n" % e)
    gE = gN = None
    have3d = False

print("NKG U grid: %dx%d  range %.2f..%.2f mm/yr" % (gU.nlat, gU.nlon, gU.vmin, gU.vmax))
if have3d:
    print("NKG E grid: range %.2f..%.2f mm/yr" % (gE.vmin, gE.vmax))
    print("NKG N grid: range %.2f..%.2f mm/yr" % (gN.vmin, gN.vmax))
    hmax = max(abs(gE.vmin), abs(gE.vmax), abs(gN.vmin), abs(gN.vmax))
    if hmax > 10.0:
        print("\n*** WARNING: horizontal rates reach %.1f mm/yr. This looks like"
              " ITRF (plate\n    motion included), NOT a Eurasia-fixed frame."
              " EGMS is ETRS89-referenced,\n    so the horizontal components must"
              " be converted before use. ***\n" % hmax)
    else:
        print("horizontal magnitudes look intraplate (max %.1f mm/yr) - OK\n" % hmax)

wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")
ctx = QgsProject.instance().transformContext()

geom_note = {}
print("%-13s %6s %7s %8s %8s %8s %9s %9s %8s %7s %6s" %
      ("case", "inc", "head", "vU", "vE", "vN", "pred3D", "predU", "obs",
       "res3D", "resU"))
print("-" * 104)

for case, (lname, orbit) in CASES.items():
    lyr = get_layer(lname)
    if lyr is None:
        print("%-13s MISSING LAYER" % case); continue
    vf = resolve_field(lyr, EGMS_VEL_FIELDS)
    incf = resolve_field(lyr, INC_FIELDS)
    headf = resolve_field(lyr, HEAD_FIELDS)
    if vf is None or incf is None:
        print("%-13s missing velocity/incidence field" % case); continue

    tr = QgsCoordinateTransform(lyr.crs(), wgs84, ctx)
    nominal = NOMINAL_HEADING[orbit]

    incs, heads, vUs, vEs, vNs, head_miss = [], [], [], [], [], []
    p3, pU, obs, r3, rU = [], [], [], [], []

    for feat in lyr.getFeatures():
        geom = feat.geometry()
        if geom is None or geom.isEmpty():
            continue
        v, inc = feat[vf], feat[incf]
        if v is None or inc is None:
            continue
        try:
            v = float(v); inc = as_degrees(float(inc))
        except (TypeError, ValueError):
            continue
        if not (10.0 < inc < 80.0):
            continue

        head, hmiss = nominal, 0.0
        if headf is not None and feat[headf] is not None:
            try:
                head, hmiss = resolve_heading(float(feat[headf]), nominal)
            except (TypeError, ValueError):
                head, hmiss = nominal, 0.0
        head_miss.append(hmiss)

        p = geom.asPoint()
        ll = tr.transform(QgsPointXY(p.x(), p.y()))
        vu = gU.sample(ll.y(), ll.x())
        if vu is None:
            continue
        eU, nU, uU = los_unit(inc, head)

        pred_u = vu * math.cos(math.radians(inc))
        if have3d:
            ve = gE.sample(ll.y(), ll.x()); vn = gN.sample(ll.y(), ll.x())
            if ve is None or vn is None:
                continue
            pred_3 = ve * eU + vn * nU + vu * uU
            vEs.append(ve); vNs.append(vn)
        else:
            pred_3 = pred_u

        incs.append(inc); heads.append(head); vUs.append(vu)
        p3.append(pred_3); pU.append(pred_u); obs.append(v)
        r3.append(v - pred_3); rU.append(v - pred_u)

    if len(obs) < 10:
        print("%-13s only %d usable points" % (case, len(obs))); continue

    m3, mu = st.median(r3), st.median(rU)
    mh = st.median(heads); mi = st.median(incs)
    eU_, nU_, uU_ = los_unit(mi, mh)
    exp_e = -0.6 if orbit == "ascending" else 0.6
    warn = ""
    if (eU_ < 0) != (exp_e < 0):
        warn = ("   <-- LOS east component %+0.2f has the WRONG SIGN for %s; "
                "the heading field may be a LOOK azimuth, not a flight heading"
                % (eU_, orbit))
    if head_miss and st.median(head_miss) > 30.0:
        warn += ("   <-- heading differs from nominal by %.0f deg"
                 % st.median(head_miss))
    geom_note[case] = ("   LOS unit vector (E,N,U) = (%+.3f, %+.3f, %+.3f) "
                       "at inc %.1f, head %.1f%s" % (eU_, nU_, uU_, mi, mh, warn))
    print("%-13s %6.1f %7.1f %8.2f %8s %8s %9.2f %9.2f %8.2f %7.2f %6.2f" %
          (case, st.median(incs), st.median(heads), st.median(vUs),
           ("%.2f" % st.median(vEs)) if vEs else "-",
           ("%.2f" % st.median(vNs)) if vNs else "-",
           st.median(p3), st.median(pU), st.median(obs), m3, mu))

print()
for c, note in geom_note.items():
    print("%-13s%s" % (c, note))
print("\nExpected LOS east component: about -0.60 ascending, +0.60 descending.")

print("\npred3D / res3D : full (E,N,U) projection")
print("predU  / resU  : vertical-only, as reported previously")
print("If res3D and resU differ noticeably, the horizontal GIA components matter")
print("and the 3-D projection is the one to report. Note that a horizontal rate")
print("projects with OPPOSITE sign on ascending and descending tracks, so it can")
print("shift the descending corridors relative to the ascending ones.")