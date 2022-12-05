import ast
import itertools
import logging
import os
import subprocess

import numpy as np
import scipy
from scipy.special import comb
import scipy.sparse as ss
import pandas as pd
from scipy import stats
from konect_scraper import config
from konect_scraper.sql import is_bipartite, single_val_get
from konect_scraper.util import get_n_m, get_directed, get_n, get_m, get_webgraph_jars
from scipy.sparse.linalg import eigs
from scipy.sparse.linalg import eigsh
from scipy.sparse.csgraph import laplacian, connected_components
import igraph as ig
import networkx as nx
from numba import njit, prange
from scipy.special import comb


def percentile_effective_diameter(p, cumulative_distances):
    """Given a percentile and an array of cumulative distances, compute
    p-percentile effective diameter: the minimal distance such that (p*100)% of 
    node pairs are at most at that distance from each other.

    Args:
        p (float): (0, 1) - a fraction corresponding to the percentile 
            (e.g 0.5 == 50%)
        cumulative_distances (list[Int]): a list of cumulative distances:
        the value at index i is the estimate of the number of pairs of nodes 
        (x, y) s.t. that the distance from x to y is _at most_ i


    Returns:
        float: the p-percentile effective diameter
    """
    # the pth percentile lies between perc_idx and perc_idx + 1
    perc_val = int((cumulative_distances[-1] + 1) * p)
    perc_idx = np.searchsorted(cumulative_distances, perc_val) - 1
    start_x = cumulative_distances[perc_idx]
    end_x = cumulative_distances[perc_idx + 1]
    interp = (perc_val - start_x) / (end_x - start_x)
    return perc_idx + interp


def algebraic_connectivity(lcc_mat):
    """The Algebraic Connectivity equals the second smallest nonzero eigenvalue
    of the Laplacian matrix of a directed graph's Giant Connected Component

    Args:
        lcc_mat (scipy.sparse.csr_matrix): A CSR matrix of the largest 
        connected component of a directed graph
    """
    logging.info("\tLaplacian..")
    lap = get_laplacian(lcc_mat.astype(np.int64))
    logging.info("\tComputed Laplacian. Computing Eigenvalues..")

    egvals, eg_vecs = get_eigenvalues(lap.astype('float'), which='SM')

    return egvals[1].real


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


def get_eigenvalues(mat, k=2, which='LM'):
    if which == 'SM':
        try:
            return eigs(mat, k, which='LM', sigma=0.0)
        except RuntimeError:
            return eigs(mat, k, which='SM')
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
    G = ig.Graph.Read_Edgelist(path, directed=True)
    g = G.to_networkx()
    # g = nx.read_edgelist(path, create_using=nx.DiGraph)

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
    n_ccs, comp = connected_components(
        csr_mat, connection='strong', return_labels=True)
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


def save_as_scipy_csr(graph_dir, edgelist_path, mat_name):
    settings = config.settings
    scipy_csr_suffix = settings['scipy_csr_suffix']
    mat_path = os.path.join(
        graph_dir, f'{mat_name}.{scipy_csr_suffix}')
    csr_mat = read_nx_edgelist_to_csr(edgelist_path)
    save_csr_matrix(mat_path, csr_mat)
    return


def parse_plfit(arr_path, plfit_stats, deg_type_str):
    settings = config.settings
    args = ['plfit', '-M', '-p', 'exact', arr_path]
    res = subprocess.check_output(args)
    logging.info(f"Running: {' '.join(args)}")
    stat_str = res.decode('utf-8')
    moment_stat_strs, mle_stat_strs = stat_str.split('\tDiscrete MLE')
    moment_stat_strs = [s.replace('\t', '')
                        for s in moment_stat_strs.split('\n')]
    mle_stat_strs = [s.replace('\t', '') for s in mle_stat_strs.split('\n')]

    moment_stat_kvs = settings['plfit']['stats']['moments']
    mle_stat_kvs = settings['plfit']['stats']['mle']

    for kvs, strs in zip([moment_stat_kvs, mle_stat_kvs], [moment_stat_strs, mle_stat_strs]):
        for stat_str in strs:
            for k, v in kvs.items():
                if k in stat_str:
                    stat_str = stat_str \
                        .replace(k, '') \
                        .strip() \
                        .replace('=', '') \
                        .strip()
                    if 'nan' in stat_str:
                        plfit_stats[f'{deg_type_str}_{v}'] = np.nan
                    else:
                        plfit_stats[f'{deg_type_str}_{v}'] = ast.literal_eval(stat_str)
    return plfit_stats


def compute_plfit_stats(graph_name):
    settings = config.settings
    graphs_dir = settings['graphs_dir']
    graph_dir = os.path.join(graphs_dir, graph_name)
    in_degs_path = os.path.join(graph_dir, f'in_degs')
    out_degs_path = os.path.join(graph_dir, f'out_degs')
    degs_path = os.path.join(graph_dir, f'degs')
    plfit_stats = {}
    paths = [in_degs_path, out_degs_path, degs_path]
    deg_type_strs = ['in', 'out', 'undir']
    for path, deg_type_str in zip(paths, deg_type_strs):
        plfit_stats = parse_plfit(path, plfit_stats, deg_type_str)
    return plfit_stats


# @njit(parallel=True)
# todo -- the gini coefficient of a graph depends on the id assignment
# todo this value will vary between diffferent isomorphism
def gini_coef(deg_arr, deg_str):
    ids = np.arange(0, deg_arr.shape[0])
    deg_arr = np.sort(deg_arr)


@njit(parallel=True)
def edge_dist_entropy(deg_arr, m):
    n = deg_arr.shape[0]
    lnn_recip = 1 / np.log(n)
    degs_2m = deg_arr / (2 * m)
    H_er = lnn_recip * (-1 * degs_2m * np.log(degs_2m)).sum()
    return H_er


def reciprocity(csr_mat):
    cx = csr_mat.tocoo()
    total = 0
    recip = 0
    for i, j, v in zip(cx.row, cx.col, cx.data):
        # print(i, j, v)
        total += 1
        if csr_mat[j, i]:
            recip += 1

    return recip / total


def balanced_inequality_ratio(deg_arr):
    total = deg_arr.sum()
    n = deg_arr.shape[0]
    vs = np.arange(1, n + 1)
    vs = vs / n
    deg_arr[::-1].sort()
    cumsum_deg = np.cumsum(deg_arr)
    cumsum_prop = 1 - (cumsum_deg / total)
    # for x, y in zip(vs, cumsum_prop):
    #     print(f"{x, y=}")


def in_out_deg_corr(out_degs, in_degs):
    inds = np.log(1 + in_degs)
    outds = np.log(1 + out_degs)
    return stats.pearsonr(inds, outds)


def compute_deg_stats(graph_name):
    settings = config.settings
    graphs_dir = settings['graphs_dir']
    compressed_fname = settings['compressed_el_file_name']
    edgelist_file_suffix = settings['edgelist_file_suffix']
    m = get_m(graph_name)
    graph_dir = os.path.join(graphs_dir, graph_name)
    in_degs_path = os.path.join(graph_dir, f'in_degs')
    out_degs_path = os.path.join(graph_dir, f'out_degs')
    degs_path = os.path.join(graph_dir, f'degs')
    out_degs = np.loadtxt(out_degs_path, dtype=np.uint32)
    in_degs = np.loadtxt(in_degs_path, dtype=np.uint32)
    degs = np.loadtxt(degs_path, dtype=np.uint32)
    deg_stats = {}
    for deg_arr, deg_str in zip(
            [out_degs, in_degs, degs],
            ['out', 'in', 'undir']):
        deg_stats[f'{deg_str}_max'] = np.max(deg_arr)
        deg_stats[f'{deg_str}_median'] = np.median(deg_arr)
        # gini_coef(deg_arr, deg_str)
        if deg_str == 'undir':
            deg_stats['relative_edge_distribution_entropy'] = edge_dist_entropy(deg_arr, m)

    coef, pval = in_out_deg_corr(out_degs, in_degs)
    deg_stats['in_out_degree_corr_coef'] = coef
    deg_stats['in_out_degree_corr_pval'] = pval
    return deg_stats


def compute_distances(adj_mat):
    print(adj_mat)
    dist_matrix = ss.csgraph.dijkstra(adj_mat, directed=True)
    dist_matrix = dist_matrix[dist_matrix != np.inf]
    print(dist_matrix)
    n = dist_matrix.shape[0]
    print(f'{np.max(dist_matrix)=}')
    print(f'{np.mean(dist_matrix)=}')
    print(f'{np.median(dist_matrix)=}')

    return


def compute_radii(graph_name, order_str='dbg'):
    settings = config.settings
    order_str_dict = settings['dbg']['order_str_dict']
    order_idx_dict = settings['dbg']['order_idx_dict']
    degree_used_for_reordering = settings['dbg']['degree_used_for_reordering']
    dbg_order_idx = order_idx_dict[order_str_dict[order_str]]
    dbg_apps_dir = settings['dbg_apps_dir']
    dbg_home = settings['dbg_home']
    dbg_datasets_dir = settings['dbg_datasets_dir']
    make_executable = settings['make_executable']
    max_iters = settings['dbg']['max_iters']
    sqlite3_db_path = settings['sqlite3']['sqlite3_db_path']
    order_file = os.path.join(
        dbg_datasets_dir, f"{graph_name}.{dbg_order_idx}.map")

    graphs_dir = settings['graphs_dir']
    graph_dir = os.path.join(graphs_dir, graph_name)
    order_file = os.path.join(graph_dir, f"{dbg_order_idx}.map")
    # order_file = os.path.join(dbg_datasets_dir, f"{graph_name}.{dbg_order_idx}.map")
    n = get_n(graph_name)
    m = get_m(graph_name)
    args = [
        make_executable,
        f"REORDERING_ALGO={dbg_order_idx}",
        f"DEGREE_USED_FOR_REORDERING={degree_used_for_reordering}",
        f"DATASET={graph_name}",
        f"MAXITERS={max_iters}",
        f"ORDER_FILE={order_file}",
        f"NUM_VERTICES={n}",
        f"NUM_EDGES={m}",
        f"GRAPH_NAME={graph_name}",
        f"SQLITE_DB_PATH={sqlite3_db_path}",
        f"GPATH={graph_dir}",
        f"DBG_ROOT={dbg_home}",
        f"COMPRESSED_EDGE_LIST_FNAME={settings['compressed_el_file_name']}",
        "run-Radii",
    ]
    env = {
        **os.environ,
        # "DBG_ROOT": dbg_home,
    }
    print(f'{dbg_home=}')
    print(f'{dbg_apps_dir=}')
    logging.info(f"Executing: " + ' '.join(args))
    print(f"Executing: " + ' '.join(args))
    res = subprocess.check_output(args, cwd=dbg_apps_dir, env=env)
    print(res.decode('utf-8'))


def get_count_of_most_common_elt(a):
    values, counts = np.unique(a, return_counts=True)
    return np.max(counts)

def symmetrize(csr_mat):
    lil_mat = csr_mat.tolil()
    rs, cs = lil_mat.nonzero()
    symm_lil_mat = lil_mat.copy()
    symm_lil_mat[cs, rs] = csr_mat[rs, cs]
    
    return symm_lil_mat.tocsr()

def compute_scipy_stats(graph_name):
    """Compute graph's summary stats using scipy. Namely:
    - directed degree statistics moments (mean, variance, std. dev, etc)
    - distances summary stats (mean, median distance)
    - algebraic stats (cyclic eigenvalue, op 2 norm)

    Args:
        graph_name (string): name of graph to compute
    """
    settings = config.settings
    graphs_dir = settings['graphs_dir']
    compressed_fname = settings['compressed_el_file_name']
    edgelist_file_suffix = settings['edgelist_file_suffix']
    scipy_csr_suffix = settings['scipy_csr_suffix']
    graph_dir = os.path.join(graphs_dir, graph_name)
    graph_path = os.path.join(
        graph_dir, f'{compressed_fname}.{edgelist_file_suffix}')

    lcc_path = os.path.join(graph_dir, f'lcc.{edgelist_file_suffix}')
    lscc_path = os.path.join(graph_dir, f'lscc.{edgelist_file_suffix}')

    nx_graph = nx.read_edgelist(graph_path, create_using=nx.DiGraph)
    mat_path = os.path.join(
        graph_dir, f'{compressed_fname}.{scipy_csr_suffix}')
    lcc_mat_path = os.path.join(graph_dir, f'lcc.{scipy_csr_suffix}')
    lscc_mat_path = os.path.join(graph_dir, f'lscc.{scipy_csr_suffix}')

    mat_names = ['comp', 'lcc', 'lscc']
    edge_list_paths = [graph_path, lcc_path, lscc_path]
    for mat_name, path in zip(mat_names, edge_list_paths):
        # if not os.path.isfile(mat_path):
        save_as_scipy_csr(graph_dir, path, mat_name)

    csr_mat = load_csr_matrix(mat_path)
    logging.info("\tComputing Reciprocity..")
    recip = reciprocity(csr_mat)
    print(f"{recip=}")
    logging.info("\tComputing SCCs..")
    n_sccs, scc_comp = ss.csgraph.connected_components(csr_mat, directed=True,
                                                       connection='strong', return_labels=True)

    logging.info("\tComputing CCs..")
    n_ccs, cc_comp = ss.csgraph.connected_components(csr_mat, directed=True,
                                                     connection='weak', return_labels=True)

    lscc_size = get_count_of_most_common_elt(scc_comp)
    lcc_size = get_count_of_most_common_elt(cc_comp)

    lcc_mat = load_csr_matrix(lcc_mat_path)
    lscc_mat = load_csr_matrix(lscc_mat_path)
    # compute_distances(csr_mat)
    n = csr_mat.shape[0]
    m = csr_mat.count_nonzero()
    # symmetrize directed csr
    logging.info("\tSymmetrizing..")
    symm_csr_mat = symmetrize(csr_mat)

    # symm_csr_mat = csr_mat.copy()
    # rs, cs = csr_mat.nonzero()
    # symm_csr_mat[cs, rs] = csr_mat[rs, cs]
    # test_eq = (symm_csr_mat != tmp).nnz == 0  
    # print(f'{test_eq =}')

    logging.info("\tSymmetric Eigenvalues..")
    # symm_egvals, symm_eg_vecs = get_eigenvalues(symm_csr_mat.astype('float'))
    symm_egvals, symm_eg_vecs = eigsh(
        A=symm_csr_mat.astype('float'),
        k=2,
        which='LM'
    )

    logging.info("\tEigenvalues..")
    egvals, eg_vecs = get_eigenvalues(csr_mat.astype('float'))

    u, s, vh = ss.linalg.svds(csr_mat.astype('float'))
    op_2_norm = s[-1]
    cyclic_eval = np.abs(egvals[0])
    spectral_norm = symm_egvals[0].real
    spectral_sep = symm_egvals[0].real / symm_egvals[1].real
    n = csr_mat.shape[0]
    m = csr_mat.count_nonzero()
    f = m / ((n * (n - 1)) / 2)
    # al_conn = algebraic_connectivity(lcc_mat)
    stats = {
        'n_vertices': n,
        'n_edges': m,
        'op_2_norm': op_2_norm,
        'cyclic_eval': cyclic_eval,
        # 'al_conn': al_conn,
        'spectral_norm': spectral_norm,
        'spectral_separation': spectral_sep,
        'fill': f,
        'n_sccs': n_sccs,
        'n_ccs': n_ccs,
        'lscc_size': lscc_size,
        'lcc_size': lcc_size,
    }

    return stats


def hyperball(graph_name):
    settings = config.settings
    graph_dir = os.path.join(settings['graphs_dir'], graph_name)
    webgraph_dir = config.settings['webgraph_dir']
    webgraph_jars = get_webgraph_jars(webgraph_dir)

    heap_size = ''
    if 'slurm_params' in settings:
        mem = int(settings['slurm_params']['mem'].replace('G', ''))
        if mem == 187:
            heap_size = 160
        else:  
            heap_size = int(settings['slurm_params']['mem'].replace('G', '')) - 4
    else:
        heap_size = 8
    print(f"{heap_size=}")
    command = f"""
        java \
        -Xss256K \
        -Xms{heap_size}G \
        -XX:PretenureSizeThreshold=512M \
        -XX:MaxNewSize=4G \
        -XX:+UseConcMarkSweepGC \
        -XX:CMSInitiatingOccupancyFraction=99 \
        -XX:+UseCMSInitiatingOccupancyOnly \
        -XX:+UseNUMA \
        -XX:+UseTLAB \
        -XX:+ResizeTLAB \
        it.unimi.dsi.webgraph.algo.HyperBall -l 10 {graph_dir}/webgraph -n {graph_dir}/neigh_func 
    """
    logging.info(f"Running: {command}..")
    env = os.environ.copy()
    if 'CLASSPATH' in env.keys():
        env['CLASSPATH'] = ':'.join(webgraph_jars) + ':' + env['CLASSPATH']
    else:
        env['CLASSPATH'] = ':'.join(webgraph_jars) 
    ret = subprocess.run(command, capture_output=True, shell=True, cwd=webgraph_dir, env=env)
    if ret.returncode == 1:  # failed
        print(f"{ret.stderr.decode()=}")
        raise Exception(f"{command} did not complete!")
    res = ret.stdout.decode()
    logging.info(res)


def compute_distance_stats(graph_name):
    settings = config.settings
    graph_dir = os.path.join(settings['graphs_dir'], graph_name)
    hyperball(graph_name)

    # the neighbourhood function of a graph is the function returning for
    # each t the number of pairs of nodes at distance at most t,
    neighbourhood_fn = np \
        .loadtxt(f"{graph_dir}/neigh_func") \
        .astype(int)
    diff = np.diff(neighbourhood_fn)
    diff = diff[diff != 0]
    n_dists = neighbourhood_fn[0]
    total = 0
    med_idx = int((neighbourhood_fn[-1] + 1) / 2)
    med_dist = np.searchsorted(neighbourhood_fn, med_idx)
    for i, v in enumerate(diff):
        total += (i + 1) * v
        n_dists += v
    mean_dist = total / n_dists
    diam = len(diff)
    perc_50 = percentile_effective_diameter(0.5, neighbourhood_fn)
    perc_90 = percentile_effective_diameter(0.9, neighbourhood_fn)
    stats = {
        'diameter': diam,
        'mean_distance': mean_dist,
        'median_distance': med_dist,
        'percentile_effective_diameter_50': perc_50,
        'percentile_effective_diameter_90': perc_90,
    }
    return stats


def n_chords_in_k_cycle(k):
    return (k * (k - 3)) / 2


def max_n_cycles(n, k):
    """
    Args:
        n
        k

    Returns:
        The maximum number of k cycles on a graph with n vertices
    """
    # todo

    return


def normalization_constant(n, k, name):
    # hack - only compute for 4motifs - will need to generalize if more
    # motifs computed
    if name.endswith('wedge'):
        return (n * (n - 1) * (n - 2)) / 2

    elif name == 'tailed_triangle':
        return comb(n, 3) * 3 * (n - 3)

    elif name == 'triangle':
        return comb(n, k)

    elif 'cycle' in name:
        n_cycles = (n * (n - 1) * (n - 2) * (n - 3)) / 4
        if 'chordal' in name:
            return n_cycles * n_chords_in_k_cycle(k)

        else:
            return n_cycles

    elif 'star' in name:
        return n * comb(n - 1, k - 1)

    elif 'path' in name:
        return (n * (n - 1) * (n - 2) * (n - 3)) / 2

    elif 'clique' in name:
        return comb(n, k)

    else:
        raise "invalid motif"
        return


def parse_motifs(k, s, n):
    d = config.settings['peregrine']['motifs'][str(k)]
    ls = s.split('\n')[3:-1]

    stats = {}
    for line in ls:
        motif_name = d[line.split(':')[0]]
        n_motifs = int(line.split(':')[1].strip())
        stats[motif_name] = n_motifs
        # norm_const = normalization_constant(n, k, motif_name)
        # print(f"{norm_const=}")
    return stats


def compute_motif_stats(graph_name):
    """
    Using peregrine, compute counts of all 3, 4-motifs in the graph
    Args:
        graph_name:

    Returns:
    """
    settings = config.settings
    graph_dir = os.path.join(settings['graphs_dir'], graph_name)
    p_data = os.path.join(graph_dir, 'peregrine')
    prgrn_count_executable = settings['peregrine']['prgrn_count_executable']
    n_threads = settings['n_threads']
    n = get_n(graph_name)
    stats = {}
    for k in [3, 4]:
        args = [prgrn_count_executable, p_data, f'{k}-motifs', str(n_threads)]

        print(f"Executing: " + ' '.join(args))
        res = subprocess.check_output(args)
        stats.update(parse_motifs(k, res.decode('utf-8'), n))
        print(res.decode('utf-8'))

    # compute edge_induced wedge
    tmp_pattern = './pattern'
    with open(tmp_pattern, 'w') as f:
        f.write('1 3\n')
        f.write('2 3\n')

    args = [prgrn_count_executable, p_data, tmp_pattern, str(n_threads)]
    res = subprocess.check_output(args)
    stats.update(parse_motifs(3, res.decode('utf-8'), n))
    print(res.decode('utf-8'))
    os.remove(tmp_pattern)
    stats['global_clustering_coefficient'] = (3 * stats['triangle']) / stats['edge_induced_wedge']
    return stats


def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)
