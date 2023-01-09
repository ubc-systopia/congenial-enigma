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
import struct
import scipy.sparse as ss


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


def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

def file_exists_and_is_nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0

def load_quads(filename, is_rect=False):
    quads = []
    print(filename)
    with open(filename, 'rb') as f:
        n_quads = struct.unpack('I', f.read(4))[0]
        q_side_len = struct.unpack('I', f.read(4))[0]
        wing_width = struct.unpack('I', f.read(4))[0]
        n = struct.unpack('I', f.read(4))[0]
        m = struct.unpack('L', f.read(8))[0]
        for i in range(n_quads):
            qx = struct.unpack('I', f.read(4))[0]
            qy = struct.unpack('I', f.read(4))[0]
            if is_rect:
                q_idx = struct.unpack('I', f.read(4))[0]
            nnz = struct.unpack('I', f.read(4))[0]
            edges = np.fromfile(
                f,
                dtype=np.dtype('u4'),
                count=nnz*2
            ).reshape(nnz, 2)
            if is_rect:
                quads.append((qx, qy, q_idx, nnz, edges))
            else:
                quads.append((qx, qy, nnz, edges))
            if qx == 5:
                for src, dest in edges:
                    if qx * q_side_len + src == 1282:
                        # print(src, dest)
                        if 'rw' in filename:
                            print(qy * q_side_len +dest+wing_width)
                            if qy * q_side_len +dest+wing_width == 2940:
                                print(src, dest)
                            print(qx, qy)
                        else:
                            print(qy * q_side_len +dest)

    return quads, q_side_len, wing_width, n, m


def load_mat(filename):
    with open(filename, 'rb') as f:
        rows = struct.unpack('q', f.read(8))[0]
        cols = struct.unpack('q', f.read(8))[0]
        nnzs = struct.unpack('q', f.read(8))[0]
        outS = struct.unpack('q', f.read(8))[0]
        innS = struct.unpack('q', f.read(8))[0]
        indptr = np.fromfile(
            f, dtype=np.dtype('u4'), count=outS).reshape(outS)
        indices = np.fromfile(f, dtype=np.dtype('u4'), count=nnzs).reshape(nnzs)
        logging.info(f"{filename=}: {rows=} {cols=} {nnzs=} {outS=} {innS=}")
        # print(f'{indices=}')
        # print(f'{indptr=}')
        assert(rows == cols)

        M = ss.csr_matrix((
            np.ones(nnzs), # 1s as data seems redundant
            indices,
            # need to append the end of the indptr array to denote the last
            # possible elt to point to
            np.concatenate((indptr, [nnzs])), 
        ), 
        shape=(rows, cols)
        )
        assert M.shape[0] == M.shape[1]
        return M
        # data = np.fromfile(f, dtype=np.float64).reshape((rows, cols))
    # return data
