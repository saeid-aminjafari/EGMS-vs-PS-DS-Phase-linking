# =============================================================================
# Switch the EGMS-equivalent RMSE filter on or off for all PS+DS layers.
#                                                               (QGIS console)
# All analysis scripts read features through QGIS, so they respect the layer
# filter. Run this with RMSE_MAX = 4.7 before re-running an analysis for the
# "RMSE-filtered" version, and with RMSE_MAX = None to go back to all points.
# Requires the field rmse_c3s (run rmse_egms_model.py first).
# =============================================================================
from qgis.core import QgsProject

RMSE_MAX = 4.7      # mm (largest EGMS-PS RMSE in our corridors); None = no filter
FIELD    = "rmse_c3s"
PSDS = ["PSDS_kir_coh70_SW99", "PSDS_le_coh70_SW99", "PSDS_lo_coh70_SW99",
        "PSDS_gbg_coh70_sw99_cor50", "PSDS_sthlm_coh70_sw99_cor50"]

for name in PSDS:
    lyr = QgsProject.instance().mapLayersByName(name)[0]
    lyr.setSubsetString("")
    n_all = lyr.featureCount()
    if RMSE_MAX is not None:
        lyr.setSubsetString('"%s" <= %s' % (FIELD, RMSE_MAX))
    n = lyr.featureCount()
    print("%-30s %s  kept %6d of %6d (%.1f%%)"
          % (name, ("RMSE <= %.1f mm" % RMSE_MAX) if RMSE_MAX else "no filter",
             n, n_all, 100.0 * n / max(n_all, 1)))