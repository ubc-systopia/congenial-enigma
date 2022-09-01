import sqlite3
import os
import json
from pathlib import Path
import logging
import igraph as ig
import numpy as np
import matplotlib.pyplot as plt
from konect_scraper import config


def __init__(self):
    config.init()


def read_iso_map(path):
    with open(path) as f:
        content = f.readlines()

    n = int(content[0])
    m = int(content[1])

    iso_map = {}
    for line in content[2:]:
        x, y = list(map(int, line.split()))
        iso_map[x] = y
    return iso_map


def get_adj_mat_from_edge_list(path, directed):
    g = ig.Graph.Read_Edgelist(path, directed=directed)

    return np.array(g.get_adjacency().data).astype(np.uint32)

def save_spy_plots(fig, axes, graph_names, order_names, path):
    settings = config.settings
    dpi = settings['plot']['dpi']

    # label the rows and columns

    for ax, col in zip(axes[0], order_names):
        ax.set_title(col)

    for ax, row in zip(axes[:, 0], graph_names):
        ax.set_ylabel(row, rotation=45, size='large')

    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    return