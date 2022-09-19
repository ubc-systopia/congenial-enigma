import subprocess
from shapely.ops import unary_union
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from shapely.geometry import Polygon
import geopandas as gpd
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from shapely.geometry import Polygon
from PIL import Image


class PatchPolygon:

    def __init__(self, polygon, **kwargs):
        polygon_path = self.pathify(polygon)
        self._patch = PathPatch(polygon_path, **kwargs)

    @property
    def patch(self):
        return self._patch

    @staticmethod
    def pathify(polygon):
        ''' Convert coordinates to path vertices. Objects produced by Shapely's
            analytic methods have the proper coordinate order, no need to sort.

            The codes will be all "LINETO" commands, except for "MOVETO"s at the
            beginning of each subpath
        '''
        vertices = list(polygon.exterior.coords)
        codes = [Path.MOVETO if i == 0 else Path.LINETO
                 for i in range(len(polygon.exterior.coords))]

        for interior in polygon.interiors:
            vertices += list(interior.coords)
            codes += [Path.MOVETO if i == 0 else Path.LINETO
                      for i in range(len(interior.coords))]

        return Path(vertices, codes)


# Plots a Polygon to pyplot `ax`
def plot_polygon(ax, poly, **kwargs):
    path = Path.make_compound_path(
        Path(np.asarray(poly.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in poly.interiors])

    patch = PathPatch(path, **kwargs)
    collection = PatchCollection([patch], **kwargs)

    ax.add_collection(collection, autolim=True)
    ax.autoscale_view()
    return collection


def hilbert_debug(n):
    hilbert_debug_exec = "/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/sb_furhilbert"

    args = [
        hilbert_debug_exec, "-n", str(n)
    ]
    subprocess.check_output(args)


def rename_files(dir):
    import os
    path = '/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/'
    for filename in os.listdir(path):
        print(filename)
        num, name = filename[:-4].split('_')
        num = num.zfill(4)
        new_filename = num + "_" + name + ".png"
        os.rename(os.path.join(path, filename), os.path.join(path, new_filename))


def create_shapely_line(start_x, start_y, end_x, end_y, width):
    is_vertical = start_x == end_x
    is_left = end_x < start_x
    is_right = start_x < end_x
    is_down = end_y < start_y

    if is_left:
        poly = Polygon([
            (end_x - width, end_y - width),
            (start_x - width, start_y - width),
            (start_x - width, start_y + width),
            (end_x - width, end_y + width),
        ])
    elif is_right:
        poly = Polygon([
            (end_x + width, end_y - width),
            (start_x - width, start_y - width),
            (start_x - width, start_y + width),
            (end_x + width, end_y + width),
        ])
    elif is_down:
        poly = Polygon([
            (start_x + width, start_y - width),
            (end_x + width, end_y - width),
            (end_x - width, end_y - width),
            (start_x - width, start_y - width),
        ])
    elif is_vertical:
        poly = Polygon([
            (start_x + width, start_y + width),
            (end_x + width, end_y + width),
            (end_x - width, end_y + width),
            (start_x - width, start_y + width),
        ])
    else:
        poly = Polygon([
            (end_x + width, end_y - width),
            (start_x + width, start_y - width),
            (start_x + width, start_y + width),
            (end_x + width, end_y + width),
        ])
    return poly


def plot_shapely_lines(n, lines, path, width=.1):
    """
    Given the lines of a furhilbert curve traversal, plot thin rectangles (i.e. "lines") that follow
    that curve
    """

    systopia_logo_path = "/home/atrostan/Downloads/image1.png"
    img = Image.open(systopia_logo_path).resize((n, n))
    img = img.rotate(180)
    im = np.asarray(img)
    fig, ax = plt.subplots()
    n = n - 1
    square_scale = 20
    square = Polygon([
        (n + width * square_scale, 0 - width * square_scale),
        (n + width * square_scale, n + width * square_scale),
        (0 - width * square_scale, n + width * square_scale),
        (0 - width * square_scale, 0 - width * square_scale),
    ])

    polygons = []
    diff = square
    for (start_x, start_y), (end_x, end_y) in lines:
        poly = create_shapely_line(start_x, start_y, end_x, end_y, width)
        polygons.append(poly)
        # p = gpd.GeoSeries(poly)
        # plot_polygon(ax, poly, facecolor='lightblue')
        # ax.plot(p)
        # ax.plot(*curve.exterior.xy)
        # diff = diff.difference(poly)
    # ax.set_aspect('equal', 'datalim')
    curve = unary_union(polygons)
    diff = square.intersection(curve)
    #
    # for geom in curve.geoms:
    #     xs, ys = geom.exterior.xy
    #     ax.plot(xs, ys, alpha=0.5, color='b')
    #
    # xs, ys = square.exterior.xy
    # ax.fill(xs, ys, alpha=0.5, color='b')
    ax.set_xlim([0 - width * square_scale, n + width * square_scale])
    ax.set_ylim([0 - width * square_scale, n + width * square_scale])
    patch = PatchPolygon(diff, facecolor='blue', edgecolor='red', alpha=0.1, transform=ax.transData).patch
    ax.add_patch(patch)

    im = ax.imshow(im)
    im.set_clip_path(patch)
    # for geom in diff.geoms:
    # xs, ys = diff.exterior.xy
    # ax.fill(xs, ys, alpha=0.5, color='b')
    ax.axis('off')
    fig.savefig(path, dpi=200, bbox_inches='tight', pad_inches=0)
    plt.close(fig)

    return


def main():
    hilbert_debug_path = "/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/hilbert"

    # rename_files("/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/")
    width_incr = 0.001
    width = .1
    for n in range(4, 400):
        hilbert_shapely_path = f"/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/{str(n).zfill(4)}_shapely.png"
        hilbert_debug(n)

        order = np.loadtxt(hilbert_debug_path).astype(np.uint32)
        mat = np.zeros((n, n))

        fig, ax = plt.subplots()

        ax.matshow(mat)

        lines = []

        for previous, current in zip(order, order[1:]):
            lines.append([
                [previous[1], previous[0]],
                # [previous[0], previous[1]],
                [current[1], current[0]],
                # [current[0], current[1]],
            ])

        plot_shapely_lines(n - 1, lines, hilbert_shapely_path, width)
        width += width_incr
        print(width)

        # width = 2
        # lc = LineCollection(
        #     lines,
        #     linewidths=width
        # )
        # ax.add_collection(lc)
        # ax.axis('off')
        # plt.tight_layout()
        # fig.savefig(
        #     f'/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/{n}_hilbert.png',
        #     bbox_inches='tight', pad_inches=0,
        #     dpi=200)
        plt.close(fig)

    return


if __name__ == "__main__":
    main()
