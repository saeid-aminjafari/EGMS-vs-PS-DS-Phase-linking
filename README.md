# Python codes to compute and visualize point count, RMSE, and ground motion velocity along linear infrastructure.

The codes should be run in QGIS IDE. The name of the layers should be set first. All the layers have been unified with same geographic extent and coordinate systems. The layers could be inported from csv, or shp, or gpkg files with ground motion data such as those from the European Ground Motion Service (EGMS) or with export commands in MintPy.

The output of the codes are figures. If you want to save the data, you need to add a command for that.

## Stage 1 — Preparation
Project each product to SW99TM: CRS EPSG:3006 (or merge egms data: Merge vector layers)
Create chainage segments (Split line by maximum length 100)
For Lilla Edet, Lödöse, and kiruna, additionally, use this too, because I processed PSDS for an area not for a corridor:
Vector → Research Tools → Create Grid
Grid type: Rectangle (polygon)
X Spacing: 100
Y Spacing: 100
Set CRS to the same projected CRS as your layer (e.g., SWEREF99 TM).
Extent: Choose the bounding box of your point layer.
First buffer (Buffer your corridor line without 100m chainages e.g., 50 m on each side, 100m for Kiruna because the line is between E45 and Malmbanan)
Second buffer (Buffer your corridor line with 100m chainage e.g., 50 m on each side).

Tool: Vector geometry → Buffer

Only for gbg and sthlm, Use Extract by location to select InSAR points within the first buffer for each dataset.

Set Temp_coherence of egms to >0.7 and and Amplitude Dispersion to <0.42:
In the attribute table → select by expression:
 "temporal_coherence" >= 0.7 AND "amplitude_dispersion" <= 0.42
Export → Save Selected Features As
 Keep the CRS and attribute fields.
 
## Stage 2 — Visualize coverage
Visual density maps
Use Heatmap (Kernel density). (radius: 50, raster size: 0.5)

## Stage 3 — Group-based comparison (not 1-to-1 points)
A. Aggregate InSAR statistics per segment
Use Join attributes by location (summary) for each product (For Lilla Edet, Lödöse, and kiruna do the same also with 100m grid):
Input: buffer_split.gpkg
Join layer: e.g., egms_corridor.gpkg
Fields to summarize: EGMS = height_wgs84,rmse, temporal_coherence, mean velocity
PSDS = velocity, coherence, dem_error, dem, rmse 

Summaries (Count, range, mean, median, stddev, q1, q3):
Do the same for PS+DS.
Each line segment now has aggregated InSAR metrics derived from all points within its buffer.

Figure 3: Plot chainage vs. mean_vel → trend plots (also plot rmse on the right axis to see if change of velocity is related to change of rmse). Python code is: Chanage_Velo_Coh.py

B. Compare distributions
Compute coverage-weighted statistics (explanation to be included in the manuscript comes at the end):

Figure 5: 4 Columns (CDF of Coh and V_std for EGMS and PSDS) - 5 rows (five cases)
Code: Hist_rmse_veloSTD
Figure 6: The same as Fig 4 but coverage weighted
Code: Hist_rmse_veloSTD_coverageWeight.py

## Stage 4 — Landcover-based density and geolocation realism
A. Clip land cover with the first buffer:
NMD swedish land cover raster: clip raster by mask extent
Input: corridor_firstBuffer(without 100m chainage)
Overlay: landcover polygons
 Output: lc_corridor
You can load the default style in the landcover folder to know what each landcover code means.
B. for raster Swedish landcover NMD2023: “Sample raster values”
Join to features in: lc_corridor

Compare to: InSAR points (PSDS and EGMS)

Fields to summarize:
EGMS = height_wgs84,rmse, temporal_coherence, mean velocity
PSDS = velocity, coherence, dem_error, dem, rmse
Summaries (Count, range, mean, median, stddev, q1, q3)

Repeat for EGMS + PSDS
The code below calculates and plots median of rmse in each land cover and Density per land cover (urban, agri, forest, wetland, bare).
Code for raster landcover: bar_stats_LC_raster.py
Figure 7: Plot as grouped bars (or Boxplots by landcover/segment type.) Code for shape file: bar_stats_LC.py: You can tweak the code to add rmse as well



## Stage 5 — “Closeness to the corridor” analysis
Compute the distance-to-infrastructure cleanly for rail (jl_riks) and road (vl_riks) separately.
1. Use “Shortest line between features”
This finds the perpendicular distance to the polyline.
Menu:
 Vector → Analysis Tools → Shortest line between features
Parameters:
Source points: PSDS (or EGMS) InSAR points
Destination hubs: jl and vl (rail and road)

Measurement: meter


Save as: egms60_gbg_rail.gpkg (or similar)
Sometimes it fails when running as batch, in this case you need to save as shp file or run one by one and then save the layer with the option of removing fid.
Repeat for PSDS.
Repeat the process for roads
2. Run this code:  Distance_to_line.py
Figure 8. Two panels for each case study: Distance to infrastructure distribution, bar chart (for points closer than 20m to the centerline) and box plot
Figure 9. Two panels for each case study: The relationship between the distance to the centerline and the coherence and the DEM.
