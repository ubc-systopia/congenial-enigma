import numpy as np
import scipy
from scipy.special import comb
import scipy.sparse as ss
import pandas as pd
from konect_scraper.sql import is_bipartite, single_val_get
from konect_scraper.util import get_n_m, get_directed
from scipy.sparse.linalg import eigs
from scipy.sparse.csgraph import laplacian, connected_components
import igraph as ig
import networkx as nx

def algebraic_connectivity(csr_mat):
    """The Algebraic Connectivity equals the second smallest nonzero eigenvalue
    of the Laplacian matrix of a directed graph's Giant Connected Component

    Args:
        csr_mat (scipy.sparse.csr_matrix): A CSR matrix of a directed graph
    """

    ccs = connected_components(csr_mat, connection='weak', return_labels=True)
    print(f'{ccs=}')
    L = laplacian(csr_mat.astype(np.int64), dtype=np.int64)

    return


def get_laplacian(csr_mat):
    """Return the laplacian of a directed graph

    Args:
        csr_mat (scipy.sparse.csr_matrix): A CSR matrix of a directed graph
    """
    return laplacian(csr_mat, dtype=np.int64)


def spectral_separation(l1, l2):
    """The spectral separation (|λ1[A] / λ2[A]|) equals the largest absolute 
    eigenvalue of the adjacency matrix divided by the second largest absolute 
    eigenvalue.

    Args:
        l1 (np.complex128): largest eigenvalue of adj mat
        l2 (np.complex128): second largest eigenvalue of adj mat
    """

    assert l1.imag == 0.0 and l2.imag == 0.0

    return np.abs(l1) / np.abs(l2)

def get_eigenvalues(mat, k=10):
    return eigs(mat, k, which='LM', maxiter=1_000_000_000)


def read_data_file_as_coo_matrix(filename):
    "Read data file and return sparse matrix in coordinate format."
    data = pd.read_csv(filename, sep=' ', header=None, dtype=np.uint32)
    rows = data[0]  # Not a copy, just a reference.
    cols = data[1]
    print(f'{data=}')
    print(rows, cols)
    ones = np.ones(len(rows), np.uint32)
    matrix = ss.coo_matrix((ones, (rows, cols)))
    return matrix

def read_nx_edgelist_to_csr(path):
    g = nx.read_edgelist(path, create_using=nx.DiGraph)
    
    return nx.to_scipy_sparse_matrix(g)

def save_csr_matrix(filename, matrix):
    """Save compressed sparse row (csr) matrix to file.

    Based on http://stackoverflow.com/a/8980156/232571

    """
    assert filename.endswith('.npz')

    ss.save_npz(filename, matrix)


def load_csr_matrix(filename):
    """Load compressed sparse row (csr) matrix from file.

    Based on http://stackoverflow.com/a/8980156/232571

    """
    assert filename.endswith('.npz')
    return ss.load_npz(filename)


def verify_stat(col, graph_name):
    """
    If the statistic "col" needs to be recomputed:
      - either Null, 0 (where inappropriate), or is different from a manually computed value
    Recompute it.
    Otherwise, return the existing Konect stat
    """
    match col:
        case 'fill':
            return verify_float(col, graph_name)


def verify_float(col, graph_name):
    existing_val = single_val_get(col, 'statistics', graph_name)
    computed_val = float(calc_stat(col, graph_name))
    print(f"{existing_val=}, {computed_val=}, {graph_name=}, {is_bipartite(graph_name)}")
    if existing_val == 0:
        return computed_val
    # elif not np.allclose(computed_val, existing_val): trust konect
    #     return computed_val
    elif is_bipartite(graph_name):
        return existing_val
    else:
        return existing_val


def calc_stat(col, graph_name):
    match col:
        case 'fill':
            return fill(graph_name)


def loops_allowed(graph_name):
    loops = single_val_get('loops', 'metadata', graph_name)
    if loops:
        return loops == 'Contains loops'
    else:
        return loops


def fill(graph_name):
    n, m = get_n_m(graph_name)
    directed = get_directed(graph_name)
    loops = loops_allowed(graph_name)
    if directed:
        n_possible_edges = comb(n, 2, exact=True)
    else:
        n_possible_edges = comb(n, 2, exact=True)
    if loops:
        n_possible_edges += n

    return m / n_possible_edges
