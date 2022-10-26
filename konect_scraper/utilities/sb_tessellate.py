import os

from matplotlib import pyplot as plt, patches

from konect_scraper import config
from konect_scraper.io import get_adj_mat_from_edge_list, read_iso_map
from konect_scraper.sql import single_val_get
from konect_scraper.util import get_directed, translate_adj_mat

import numpy as np
import pandas as pd


def __init__(self):
    config.init()


def main():
    graph_name = "maayan-faa"
    graphs_dir = "/media/atrostan/patterson_backup/data/graphs/"
    graph_dir = f"{graphs_dir}/{graph_name}"
    graph_path = f"{graph_dir}/comp.net"
    iso_path = f"{graph_dir}/parsb"
    directed = get_directed(graph_name)
    adj_mat = get_adj_mat_from_edge_list(graph_path, directed)
    iso_path = os.path.join(graph_dir, iso_path)
    iso_map = read_iso_map(iso_path)

    map_mat = translate_adj_mat(adj_mat, iso_map)
    sb_k = single_val_get('par_sb_k', 'statistics', graph_name)
    sb_n_iters = single_val_get('par_sb_n_iters', 'statistics', graph_name)
    sqr_width = sb_k * sb_n_iters
    rect = patches.Rectangle((0, sqr_width), sqr_width, -sqr_width,
                             linewidth=1, edgecolor='r', alpha=0.3, zorder=2, facecolor='r')
    sb_rows = map_mat[:sqr_width, sqr_width:]
    sb_cols = map_mat[sqr_width:, :sqr_width]
    max_sb_cols = np.argmax(sb_cols, axis=1)
    n = map_mat.shape[0]
    max_i = np.zeros(sqr_width).astype(np.uint32)
    max_j = np.zeros(n - sqr_width).astype(np.uint32)
    for i in range(sqr_width - 1, -1, -1):
        for j in range(sqr_width, n):
            if map_mat[i, j]:
                max_i[i] = j

    for i in range(sqr_width, n):
        for j in range(0, sqr_width):
            if map_mat[i, j]:
                max_j[i - sqr_width] = j

    print(max_j)
    series_i = pd.Series(reversed(max_i))
    series_j = pd.Series(reversed(max_j))
    fig, ax = plt.subplots()
    for i, j in zip(range(sqr_width - 1, -1, -1), series_i.cummax()):
        map_mat[i, int(j)] = 2

    for i, j in zip(range(n - 1, -1, -1), series_j.cummax()):
        map_mat[i, int(j)] = 2



    mn = np.full(n - sqr_width, np.inf)
    mx = np.zeros(n - sqr_width)
    for i in range(sqr_width, n):
        for j in range(sqr_width, n):
            if map_mat[i, j]:
                if j < mn[i - sqr_width]:
                    mn[i - sqr_width] = j
                if j > mx[i - sqr_width]:
                    mx[i - sqr_width] = j


    for ii, (i, j) in enumerate(zip(mn, mx)):
        print(ii, i, j)
        if i == np.inf:
            continue
        map_mat[ii + sqr_width, int(i)] = 2
        map_mat[ii + sqr_width, int(j)] = 2

    ax.spy(map_mat == 1, markersize=1, color='b')
    ax.spy(map_mat == 2, markersize=1, color='r')
    ax.add_patch(rect)
    plt.tight_layout()
    plt.savefig('./tmp', dpi=200)
    # read graph

    # read wing width ratio
    # read parsb map

    # highlight bounds on sb map
    return


if __name__ == '__main__':
    config.init()
    settings = config.settings

    main()
