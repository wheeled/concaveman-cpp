import csv
import os
import numpy as np
import pandas as pd
import geopandas as gpd
from concaveman import concaveman2d
from scipy.spatial import ConvexHull
from shapely import geometry


CRS = "epsg:4283"
INFILE = '/Users/wheeled/GEOSPATIAL/SDS boundaries/Prototype.gpkg'
# INFILE = '/Users/wheeled/PycharmProjects/concaveman-cpp/src/main/python/data/boundary_prototype.csv'
OUTFILE = '/Users/wheeled/GEOSPATIAL/SDS boundaries/Prototype.gpkg'
CSVPATH = '/Users/wheeled/PycharmProjects/concaveman-cpp/src/main/python/data'
CONCAVITY = 1.2
CORRECTIONS = os.path.join(CSVPATH, 'corrections.csv')


def unique(points):
    # eliminate duplicate geometries
    unique_points = []
    for point in points:
        if point not in unique_points:
            unique_points.append(point)
    return unique_points


def save_as_csv(filename, points):
    # save list of points as an indexed CSV file (that could be imported in QGIS)
    with open(os.path.join(CSVPATH, filename), 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        csv_writer.writerow(['id', 'X', 'Y'])
        for point in points:
            csv_writer.writerow(point)


class GeoPackageFile(object):
    # this object will contain the intermediate geometries as a link to QGIS
    def __init__(self, filepath):
        self.filepath = filepath

    def add_point_layer(self, layername, indexed_list, crs=CRS):
        df = pd.DataFrame({
            'id': [row[0] for row in indexed_list],
            'X': [row[1] for row in indexed_list],
            'Y': [row[2] for row in indexed_list],
        })
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.X, df.Y), crs=crs)
        gdf.to_file(self.filepath, layer=layername, driver="GPKG")

    def add_boundary_polygon(self, layername, indexed_list, crs=CRS):
        polygon = geometry.Polygon([[p[1], p[2]] for p in indexed_list])
        gdf = gpd.GeoDataFrame(pd.DataFrame({'geometry': [polygon]}), crs=crs)
        gdf.to_file(self.filepath, layer=layername, driver="GPKG")


def apply_corrections(points, missing):
    if not os.path.exists(CORRECTIONS):
        print("No corrections possible: 'corrections.csv' file missing")
        return points

    df = pd.read_csv(CORRECTIONS)
    array = np.array(df)

    if array.shape[0] < 3:
        print("No corrections yet or insufficient info in 'corrections.csv'")
        return points

    if np.max([_ for _ in array[:, 0] if not np.isnan(_)]) > len(points) - 1:
        print("An index value in column 0 of 'corrections.csv' exceeds the number of reordered hull points")
        return points

    if np.max([_ for _ in array[:, 1] if not np.isnan(_)]) > len(missing) - 1:
        print("An index value in column 1 of 'corrections.csv' exceeds the number of missing points")
        return points

    print('Applying corrections and inserting %s missing points' % len(missing))
    corrected_points = []
    pointer = 0
    for idx, row in enumerate(array):
        active_column = np.logical_not(np.isnan(row))
        if active_column.all() or not active_column.any():
            print("Invalid format: each row of 'corrections.csv' must have one of "
                  "hull index or missing index (%s: %s)" % (idx, row))
            return points

        if active_column[0] and row[0] >= pointer and pointer not in array[idx + 1:, 0]:
            corrected_points.extend(points[pointer:int(row[0] + 1)])
            pointer = int(row[0] + 1)

        elif active_column[0] and pointer in array[:, 0]:
            corrected_points.append(points[int(row[0])])
            pointer += 1

        elif active_column[1]:
            corrected_points.append(missing[int(row[1])])

        else:
            print("Corrections out of order at row %s: pointer (%s) has already "
                  "passed this value %s" % (idx, pointer, row))
            return points

    if pointer < len(points) - 1:
        corrected_points.extend(points[pointer:])

    # # printout before renumbering to confirm that the points were remapped correctly
    # print("Corrected points before reindexing (%s)" % len(corrected_points))
    # for point in corrected_points:
    #     print(point)

    corrected_points = [[_, p[1], p[2]] for _, p in enumerate(corrected_points)]

    return corrected_points


def main():
    infile_format = INFILE.split('.')[-1].lower()

    if infile_format == 'gpkg':
        gdf = gpd.read_file(INFILE, driver="GPKG", layer="boundary_points")
        raw_points = [[point.x, point.y] for point in gdf.geometry]
        points = unique(raw_points)
        print("Opening boundary_points with CRS = %s (%s points)" % (gdf.crs['init'], len(points)))
        if gdf.crs['init'] != CRS:
            print("GeoPackage CRS (%s) is not the same as default CRS (%s)" % (gdf.crs['init'], CRS))

    elif infile_format == 'csv':
        with open(INFILE, 'r') as csvfile:
            csv_reader = csv.reader(csvfile)
            points = [row for row in csv_reader]
            points.pop(0)
            raw_points = [[float(row[0]), float(row[1])] for row in points]
            points = unique(raw_points)

    else:
        print("INFILE format not supported")
        return

    prototype = GeoPackageFile(OUTFILE)

    # create a convex hull to seed the concave hull
    hull = ConvexHull(points)

    hull_points = [[idx, hull.points[idx][0], hull.points[idx][1]] for idx in hull.vertices]
    prototype.add_point_layer("hull_points", hull_points)
    prototype.add_boundary_polygon("hull_polygon", hull_points)

    result = concaveman2d(points, hull.vertices, concavity=CONCAVITY)
    indexed_result = [[_, *point] for _, point in enumerate(result)]
    save_as_csv('concave_boundary.csv', indexed_result)

    missing = [[_, p[0], p[1]] for _, p in enumerate([point for point in points if point not in result])]
    prototype.add_point_layer("missing_points", missing)
    save_as_csv('missing_points.csv', missing)

    print(len(raw_points), '-->', len(points), '-->', len(indexed_result))
    indexed_result = apply_corrections(indexed_result, missing)

    prototype.add_point_layer("reordered_points", indexed_result)
    prototype.add_boundary_polygon("constructed_polygon", indexed_result)

    print("Updated '%s' with 'constructed_polygon' layer (%s points)" % (
        os.path.split(OUTFILE)[1], len(indexed_result)
    ))


main()

