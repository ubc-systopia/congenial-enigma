import sqlite3
import os
import json
from pathlib import Path
import logging
import igraph as ig
import igraph._igraph
import numpy as np
import matplotlib.pyplot as plt
from konect_scraper import config
from konect_scraper.download_and_extract import get_n_comments_at_top_of_file


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
    settings = config.settings
    comments = settings['comment_strings']

    n_comments = get_n_comments_at_top_of_file(path)
    with open(path) as f:
        for i in range(n_comments):
            line = f.readline()
    n_values = len(line.split())
    dtype_cols = []
    for i in range(n_values):
        dtype_cols += [(f'c{i}', np.uint32)]
    # arr = np.fromfile(path, sep=' ').astype(np.uint32)
    arr = np.loadtxt(path, comments=comments).astype(np.uint32)
    # print(arr)
    # if arr.shape[0] % n_values != 0:  # values missing fill with zeros
    #     filled_len = int(np.ceil(arr.shape[0] / n_values) * n_values)
    #
    #     filled = np.zeros(filled_len).astype(np.uint32)
    #     filled[:arr.shape[0]] = arr
    #     arr = filled
    #
    # print(arr)
    # reshaped = arr.reshape((int(arr.shape[0] / n_values), n_values))[:, :2]
    max_vid = np.max(arr.flatten())
    adj_mat = np.zeros((max_vid + 1, max_vid + 1))

    for e in arr:
        adj_mat[e[0], e[1]] = 1
        if not directed:
            adj_mat[e[1], e[0]] = 1

    return adj_mat

def read_image(path, fmt):
    match fmt:
        case 'png':
            return plt.imread(path)
        case 'pdf':
            raise Exception("PDF format - UNIMPLEMENTED!")

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