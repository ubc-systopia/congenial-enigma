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
from numba import njit, prange

def algebraic_connectivity(csr_mat):
    """The Algebraic Connectivity equals the second smallest nonzero eigenvalue
    of the Laplacian matrix of a directed graph's Giant Connected Component

    Args:
        csr_mat (scipy.sparse.csr_matrix): A CSR matrix of a directed graph
    """

    n_ccs, comp = connected_components(csr_mat, connection='weak', return_labels=True)
    component_counts = np.unique(comp, return_counts=True)
    gcc_id = component_counts[0][np.argmax(component_counts[1])]

    L = extract_gcc_from_csr(gcc_id, comp, csr_mat)
    lap = get_laplacian(L.astype(np.int64))
    
    return


def drop_rows_from_csr(mat, rows_to_keep):
    return ss.lil_matrix(mat[rows_to_keep, :]).tocsr()

def drop_cols_from_csr(mat, cols_to_keep):
    return ss.lil_matrix(mat[:, cols_to_keep]).tocsr()

def extract_gcc_from_csr(gcc_id, comp, csr_mat):
    vids_in_gcc_ids = np.where(comp == gcc_id)[0]
    mat = drop_rows_from_csr(csr_mat, vids_in_gcc_ids)
    mat = drop_cols_from_csr(mat, vids_in_gcc_ids)
    return mat

def get_laplacian(csr_mat):
    """Return the laplacian of a directed graph

    Args:
        csr_mat (scipy.sparse.csr_matrix): A CSR matrix of a directed graph
    """
    return laplacian(csr_mat)


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

def get_eigenvalues(mat, k=10, which='LM'):
    return eigs(mat, k, which=which)


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


def controllability(graph_name):

    return 

def size_of_lscc(csr_mat):
    n_ccs, comp = connected_components(csr_mat, connection='strong', return_labels=True)
    component_counts = np.unique(comp, return_counts=True)
    return np.max(component_counts[1])

def diameter(nx_graph):
    try:
        diam = nx.diameter(nx_graph)
        return diam
    except:
        largest = max(nx.strongly_connected_components(nx_graph), key=len)
        return diameter(nx_graph.subgraph(largest))


def get_degrees(csr_mat, mode='out'):
    if mode == 'out':
        return ss.csgraph.laplacian(csr_mat, return_diag=True, use_out_degree=True)
    else:
        return ss.csgraph.laplacian(csr_mat, return_diag=True, use_out_degree=False)

@njit(parallel=True)
def power_law_estimate(ds, xmin=None):
    ds = ds[ds > 0]

    mn = np.min(ds)
    n = ds.shape[0]
    return 1 + n * (np.reciprocal(np.sum(np.log(ds / mn))))


@njit(parallel=True)
def tail_power_law_estimate(ds, xmin):
    ds = ds[ds > 0]
    # filter out all values less than xmin
    # this achieves us only taking into account the tail of the distribution
    ds = ds[ds < xmin]
    print(ds)
    mn = np.min(ds)
    n = ds.shape[0]
    return 1 + n * (np.reciprocal(np.sum(np.log(ds / mn))))

def plfit_stats(csr_mat):
    out_degrees = get_degrees(csr_mat, mode='out')[1]
    in_degrees = get_degrees(csr_mat, mode='in')[1]
    np.savetxt('./in_degs.tmp', in_degrees.astype(np.uint32), fmt='%d')
    # print(f'{out_degrees=}')
    # plfit.plfit_discrete(out_degrees)
    # plfit.plfit_discrete(out_degrees)
    # print(plfit.alpha)
    return 

def compute_scipy_stats(graph_path):
    """Compute graph's summary stats using scipy. Namely:
    - 

    Args:
        graph_path (_type_): _description_
    """
    
    return 
