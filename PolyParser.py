# Stolen from the transnet repo.
from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon
from shapely.geometry import Point


class PolyParser:
    def __init__(self):
        pass

    @staticmethod
    def poly_to_polygon(file_name):
        with open(file_name) as f:
            content = f.readlines()

        return PolyParser.parse_poly(content)

    @staticmethod
    def parse_poly(lines):
        """Parse an Osmosis polygon filter file.

        Accept a sequence of lines from a polygon file, return a shapely.geometry.MultiPolygon object.

        http://wiki.openstreetmap.org/wiki/Osmosis/Polygon_Filter_File_Format
        """
        ring = None
        in_ring = False
        coords = []

        for index, line in enumerate(lines):
            if index == 0:
                # first line is junk.
                continue

            elif index == 1:
                # second line is the first polygon ring.
                coords.append([[], []])
                ring = coords[-1][0]
                in_ring = True

            elif in_ring and line.strip() == "END":
                # we are at the end of a ring, perhaps with more to come.
                in_ring = False

            elif in_ring:
                # we are in a ring and picking up new coordinates.
                ring.append(map(float, line.split()))

            elif not in_ring and line.strip() == "END":
                # we are at the end of the whole polygon.
                break

            elif not in_ring and line.startswith("!"):
                # we are at the start of a polygon part hole.
                coords[-1][1].append([])
                ring = coords[-1][1][-1]
                in_ring = True

            elif not in_ring:
                # we are at the start of a polygon part.
                coords.append([[], []])
                ring = coords[-1][0]
                in_ring = True
        # print('coords')
        # print(coords)
        # np.save('coords',coords)

        # fixes an issue with creating multipolygon from map objects.
        assert len(coords[0]) == 2
        assert len(coords[0][1]) == 0

        newcoords = []
        # print(len(coords[0][0]))
        for shape in coords[0][0]:
            point = Point([x for x in shape])

            if point == Point():
                print("THIS IS VERY BAD")
                quit()
                return MultiPolygon([])
            newcoords.append(point)
        # np.save('newcoords',newcoords)
        # print('success')
        # print(newcoords)
        # print('newcoords')
        # print(newcoords)
        newcoords.append(
            newcoords[0]
        )  # repeat the first point to create a 'closed loop'

        poly = Polygon([[p.x, p.y] for p in newcoords])
        # plt.plot(*poly.exterior.xy)
        # plt.show()
        return MultiPolygon([poly])
