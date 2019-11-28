import csv
import os
import pandas as pd
import geopandas as gpd
from concaveman import concaveman2d
from scipy.spatial import ConvexHull
from shapely import geometry


CRS = "EPSG:4283"
INFILE = '/Users/wheeled/GEOSPATIAL/SDS boundaries/Prototype.gpkg'
# INFILE = '/Users/wheeled/PycharmProjects/concaveman-cpp/src/main/python/data/boundary_prototype.csv'
OUTFILE = '/Users/wheeled/GEOSPATIAL/SDS boundaries/Prototype.gpkg'
CSVPATH = '/Users/wheeled/PycharmProjects/concaveman-cpp/src/main/python/data'


def unique(points):
    # eliminate duplicate geometries
    unique_points = []
    for point in points:
        if point not in unique_points:
            unique_points.append(point)
    return unique_points


class GeoPackageFile(object):
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


def main():
    infile_format = INFILE.split('.')[-1]
    if infile_format == 'gpkg':
        gdf = gpd.read_file(INFILE, driver="GPKG", layer="boundary_points")
        points = unique([[point.x, point.y] for point in gdf.geometry])
        print("Opening boundary_points with CRS = %s (%s points)" % (gdf.crs, len(points)))
        if gdf.crs != CRS:
            print("GeoPackage CRS (%s) is not the same as default CRS (%s)" % (gdf.crs, CRS))


    elif infile_format == 'csv':
        with open(INFILE, 'r') as csvfile:
            csv_reader = csv.reader(csvfile)
            points = [row for row in csv_reader]
            points.pop(0)
            points = unique([[float(row[0]), float(row[1])] for row in points])

    else:
        print("INFILE format not supported")
        return

    prototype = GeoPackageFile(OUTFILE)

    # create a convex hull to seed the concave hull
    hull = ConvexHull(points)

    hull_points = [[idx, hull.points[idx][0], hull.points[idx][1]] for idx in hull.vertices]
    prototype.add_point_layer("hull_points", hull_points)
    prototype.add_boundary_polygon("hull_polygon", hull_points)

    result = concaveman2d(points, hull.vertices, concavity=1.2)
    indexed_result = [[_, *point] for _, point in enumerate(result)]
    # print(indexed_result[:10])
    print(len(points), '-->', len(points), '-->', len(indexed_result))

    with open(os.path.join(CSVPATH, 'concave_boundary.csv'), 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        csv_writer.writerow(['id', 'X', 'Y'])
        for point in indexed_result:
            csv_writer.writerow(point)

    missing = [point for point in points if point not in result]
    print("%s points missing:\n%s" % (
        len(missing),
        "".join(["  %s\n" % point for point in missing])
    ))
    with open(os.path.join(CSVPATH, 'missing_points.csv'), 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        csv_writer.writerow(['id', 'X', 'Y'])
        for _, point in enumerate(missing):
            csv_writer.writerow([_, point[0], point[1]])

    # TODO: rather than manually perform the manipulation in QGIS:
    #   insert the "missing" vertices in logical places, then

    prototype.add_point_layer("reordered_points", indexed_result)
    prototype.add_boundary_polygon("constructed_polygon", indexed_result)

    # TODO: provide a tool to do the reordering based on a csv list of old and new ids


main()

