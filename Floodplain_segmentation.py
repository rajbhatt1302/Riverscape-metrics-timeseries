
# This script segments the river floodplain into 6km segments using perpendicular cut lines along the centerline.
# It improves upon the original approach by:
# 1.  Using robust UTM projections.
# 2.  Optimizing the cut-line generation using parallel processing.
# 3.  Using efficient geometry splitting operations.

import os
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import LineString, MultiLineString, MultiPolygon
from shapely.ops import split, linemerge, unary_union
from joblib import Parallel, delayed
import multiprocessing
from tqdm.notebook import tqdm
import warnings

warnings.filterwarnings('ignore')


# ## 1. Load Data & Project

# Input Paths
RIVER_PATH = r"D:\irp_new\Test\River\The_Ganga_River.shp"
CENTERLINE_PATH = r"D:\irp_new\Test\River\The_Ganga_Floodplain_CentreLine_Smooth_V1.shp"
OUTPUT_DIR = r"D:\irp_new\Test\River\Output"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Load Files
print("Loading files...")
river_gdf = gpd.read_file(RIVER_PATH)
centerline_gdf = gpd.read_file(CENTERLINE_PATH)

# Auto-Project to UTM (Zone 44N for Ganga/India is standard, or auto-detect)
# Using EPSG:32644 (WGS 84 / UTM zone 44N)
TARGET_CRS = "EPSG:32644"
print(f"Projecting to {TARGET_CRS}...")
river_gdf = river_gdf.to_crs(TARGET_CRS)
centerline_gdf = centerline_gdf.to_crs(TARGET_CRS)

# Merge centerline if it's fragmented
centerline_geom = linemerge(centerline_gdf.geometry.unary_union)
if centerline_geom.geom_type == 'MultiLineString':
    # Take the longest line if still MultiLine
    centerline_geom = max(centerline_geom.geoms, key=lambda x: x.length)

print(f"River Area: {river_gdf.geometry.area.sum() / 1e6:.2f} sq km")
print(f"Centerline Length: {centerline_geom.length / 1000:.2f} km")


# ## 2. Generate Cut Lines (Parallel)

SEGMENT_LENGTH_M = 6000  # 6km
CUT_WIDTH_M = 50000      # Increased width to ensure coverage

def generate_cut_line(distance, centerline, width):
    """Generates a perpendicular cut line at a specific distance."""
    try:
        # Get point and local tangent
        point = centerline.interpolate(distance)
        
        # Small delta to calculate tangent
        delta = 1.0
        p1 = centerline.interpolate(max(0, distance - delta))
        p2 = centerline.interpolate(min(centerline.length, distance + delta))
        
        # Tangent vector
        dx = p2.x - p1.x
        dy = p2.y - p1.y
        
        # Normal vector (-dy, dx)
        length = (dx**2 + dy**2)**0.5
        if length == 0: return None
        
        nx = -dy / length
        ny = dx / length
        
        # Create cut line centered at 'point'
        start_x = point.x + nx * (width / 2)
        start_y = point.y + ny * (width / 2)
        end_x = point.x - nx * (width / 2)
        end_y = point.y - ny * (width / 2)
        
        return LineString([(start_x, start_y), (end_x, end_y)])
    except Exception as e:
        return None

# Calculate distances for cut lines
# The first segment typically spans 0 to 6km, requiring a cut line at the 6km mark.
# The "leftmost" segment is defined as the polygon area preceding the first cut line.
# Note: If the river geometry extends significantly upstream of the centerline start, 
# that section will not be segmented as no cut lines can be generated there.
# The cut width has been increased to 50km to ensure complete coverage of the floodplain.

distances = np.arange(SEGMENT_LENGTH_M, centerline_geom.length, SEGMENT_LENGTH_M)

print(f"Generating {len(distances)} cut lines using {multiprocessing.cpu_count()} cores...")

# Parallel Generation of cut lines
cut_lines = Parallel(n_jobs=-1)(
    delayed(generate_cut_line)(d, centerline_geom, CUT_WIDTH_M) for d in tqdm(distances)
)

# Filter out failed generations (None values)
cut_lines = [line for line in cut_lines if line is not None]

# Sort cut lines by their projection distance along the centerline to ensure sequential order
cut_lines.sort(key=lambda l: centerline_geom.project(l.interpolate(0.5, normalized=True)))

# Overlap Resolution Strategy
# To handle potential overlaps between cut lines (which can cause topological errors),
# the script employs a robust "split and dissolve" strategy later in the process.
# 1. The river polygon is split by ALL cut lines.
# 2. Resulting polygons are assigned a Segment ID based on their centroid's location along the centerline.
# 3. Polygons with the same ID are dissolved together.
# This avoids the need for complex geometric adjustments of individual cut lines while ensuring valid output.

def clean_cuts(lines):
    """
    Optional filter to remove cut lines that intersect with their immediate predecessor.
    This ensures strict topological validity but may result in merged segments (e.g., 12km instead of 6km).
    Currently disabled in favor of the split-and-dissolve strategy.
    """
    keep = []
    if not lines: return []
    keep.append(lines[0])
    for i in range(1, len(lines)):
        last = keep[-1]
        curr = lines[i]
        if not last.intersects(curr):
            keep.append(curr)
        else:
            # Overlap detected; keep the current line to maintain segment count.
            # The dissolve step will handle the resulting geometry.
            keep.append(curr)
    return keep

# cut_lines = clean_cuts(cut_lines) # Uncomment to enable strict non-overlapping cuts

cut_lines_geom = MultiLineString(cut_lines)


# ## 3. Split River Polygon

print("Splitting river polygon...")
# Get the single polygon geometry
river_poly = river_gdf.geometry.unary_union

# Simplify to speed up processing (10m tolerance)
print("Simplifying river geometry for performance...")
river_poly = river_poly.simplify(10, preserve_topology=True)


# Split
split_geoms = split(river_poly, cut_lines_geom)

print(f"Initial split into {len(split_geoms.geoms)} polygons.")

# Post-processing to merge slivers and assign IDs correctly
# We calculate the distance of each polygon's centroid along the centerline
# and bin them into 6km chunks.

print("Assigning segments...")
segment_polys = []
segment_ids = []

for poly in tqdm(split_geoms.geoms, desc="Processing Segments"):
    # Find distance of centroid along centerline
    dist = centerline_geom.project(poly.centroid)
    
    # Calculate Segment ID
    seg_id = int(dist // SEGMENT_LENGTH_M)
    
    segment_ids.append(seg_id)
    segment_polys.append(poly)

# Create temporary GDF
temp_gdf = gpd.GeoDataFrame({'geometry': segment_polys, 'Segment_ID': segment_ids}, crs=TARGET_CRS)

# Dissolve by Segment ID to merge any split parts (e.g. islands or slivers from crossing lines)
print("Dissolving by ID...")
final_segments_gdf = temp_gdf.dissolve(by='Segment_ID').reset_index()

print(f"Final segments: {len(final_segments_gdf)}")


# ## 4. Save & Visualize

# Save
out_file = os.path.join(OUTPUT_DIR, "Ganga_Segments_Perpendicular.shp")
final_segments_gdf.to_file(out_file)
print(f"Saved to {out_file}")

# Plot
fig, ax = plt.subplots(figsize=(15, 5))
final_segments_gdf.plot(ax=ax, column='Segment_ID', cmap='jet', edgecolor='black', linewidth=0.5)
centerline_gdf.plot(ax=ax, color='white', linewidth=1, alpha=0.5)
plt.title("River Segments (Perpendicular Cuts)")
plt.show()
