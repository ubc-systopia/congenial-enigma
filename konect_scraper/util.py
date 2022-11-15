import shutil
import sqlite3
import os
import json
import subprocess
from pathlib import Path
import logging

import igraph
import networkx as nx
import numpy as np
import psutil

from konect_scraper import config, column_names
import pandas as pd
import math
from konect_scraper.sql import connect, single_val_get


def __init__(self):
    config.init()


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def column_exists(column, table):
    """
    true, if column in table, false otherwise
    """
    conn = connect()
    sql = f"SELECT * FROM pragma_table_info('{table}') WHERE name=?"
    cursor = conn.cursor()

    cursor.execute(sql, (column,))
    res = cursor.fetchall()
    return len(res) == 1


def cast_np_dtypes(df, dtypes):
    # TODO imputing with zero - is this appropriate?
    for col in df.columns:

        # columns contain large values so cast them using larger dtypes
        if col in ['wedge_count', 'claw_count', 'cross_count', 'triangle_count']:
            # arr = np \
            #     .nan_to_num(df[col].values) \
            #     .astype(np.ulonglong)
            # df[col] = pd.to_numeric(arr, errors='coerce').astype(np.ulonglong)
            # if count is null, replace with zero
            df[col] = df[col].fillna(0)
            df[col] = df[col].astype(object)

            continue
        #     # return

        dtype = dtypes[col]
        if dtype == 'Int64':
            arr = np \
                .nan_to_num(df[col].values) \
                .astype(int)
            df[col] = pd.to_numeric(arr, errors='coerce').astype(int)
        if dtype == 'Float64':
            arr = np \
                .nan_to_num(df[col].values) \
                .astype(np.float64)
            df[col] = pd.to_numeric(arr, errors='coerce').astype(np.float64)

    return


def update_table_schema(table, col_names):
    for col in col_names:
        dtype = col_names[col]
        add_column_if_not_exists(col, table, dtype)

    return


def add_column_if_not_exists(column, table, dtype):
    # check if column exists in table
    if column_exists(column, table):
        return
    conn = connect()
    sql = f"alter table {table} add column {column} {dtype} default null"
    cursor = conn.cursor()

    cursor.execute(sql)
    conn.commit()
    conn.close()

    return


def create_pr_expt_table():
    conn = connect()

    # graph_name, datetime, expt_num should be unique
    column_dict = column_names.pr_expts_col_names
    unique_cols = ["graph_name", "datetime",
                   "expt_num", "vertex_order", "edge_order"]
    sql_create_table = f"""CREATE TABLE IF NOT EXISTS pr_expts ("""
    for k in list(column_dict.keys()):
        v = column_dict[k]
        sql_create_table += f"{k} {v.lower()}, "

    sql_create_table += f"UNIQUE({', '.join(unique_cols)}));"
    cursor = conn.cursor()
    cursor.execute(sql_create_table)
    conn.commit()


def create_sql_table(conn, table_name, column_dict):
    sql_create_table = f"""CREATE TABLE IF NOT EXISTS {table_name} ("""

    for k in list(column_dict.keys())[:-1]:
        v = column_dict[k]
        if k == 'graph_name':
            sql_create_table += f"{k} {v.lower()} PRIMARY KEY, "
            continue
        sql_create_table += f"{k} {v.lower()}, "
    k = list(column_dict.keys())[-1]
    v = column_dict[k]
    sql_create_table += f"{k} {v.lower()}"
    sql_create_table += ");"
    cursor = conn.cursor()
    cursor.execute(sql_create_table)
    return


def get_stat_dtype(col_name):
    return column_names.sql_to_np_dtypes[column_names.stat_col_names[col_name]]


def get_additional_stats(cols, graph_name):
    d = {}
    for col in cols:
        val = single_val_get(col, 'statistics', graph_name)
        d[col] = val
    return d


def delete_graphs_db():
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    logging.info(f"Deleting {db_path}..")
    try:
        os.remove(db_path)
    except OSError:
        pass


def get_graph_dir(graph_name):
    settings = config.settings
    graphs_dir = settings['graphs_dir']
    return os.path.join(graphs_dir, graph_name)


def get_plot_dir(graph_name, plot_type):
    settings = config.settings
    plots_dir = config.settings['plots_dir']
    plot_dir = os.path.join(plots_dir, graph_name)
    return os.path.join(plot_dir, plot_type)


def get_plotting_dirs(graph_name):
    graph_dir = get_graph_dir(graph_name)


def create_plot_dirs_if_not_exists(graph_name):
    plots_dir = config.settings['plots_dir']
    plot_dir = os.path.join(plots_dir, graph_name)
    adj_mat_dir = os.path.join(plot_dir, "adj_mat")
    spy_dir = os.path.join(plot_dir, "spy")

    Path(plot_dir).mkdir(parents=True, exist_ok=True)
    Path(adj_mat_dir).mkdir(parents=True, exist_ok=True)
    Path(spy_dir).mkdir(parents=True, exist_ok=True)


def remove_all_files_in_directory(folder):
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


def create_log_dir_if_not_exists():
    slurm_log_dir = config.settings['logging']['slurm_log_dir']
    data_dirs = [
        slurm_log_dir,
        config.settings['logging']['log_dir'],
    ] + [
        os.path.join(slurm_log_dir, m) for m in
        config.settings['compute_canada']['execution_modes']
    ]

    [Path(data_dir).mkdir(parents=True, exist_ok=True)
     for data_dir in data_dirs]


def create_dir_if_not_exists(path_to_dir):
    Path(path_to_dir).mkdir(parents=True, exist_ok=True)


def create_data_dirs_if_not_exists():
    data_dir_names = [
        "data_dir",
        "dataframes_dir",
        "graphs_dir",
        "plots_dir",
    ]
    data_dirs = [
        config.settings[n] for n in data_dir_names
    ]

    for data_dir in data_dirs:
        Path(data_dir).mkdir(parents=True, exist_ok=True)


def translate_adj_mat(arr, iso_map):
    n = arr.shape[0]
    map_arr = np.zeros((n, n))
    nnz = np.nonzero(arr)
    for src, dest in zip(nnz[0], nnz[1]):
        map_arr[
            iso_map[src],
            iso_map[dest]
        ] = 1

    return map_arr


def verify_graphs_in_json(konect_names):
    settings = config.settings

    datasets_json_path = settings['datasets_json_path']
    with open(datasets_json_path) as f:
        datasets = json.load(f)
    # verify that all graph names exist in datasets.json

    graph_names = set([d['name'] for d in datasets['datasets']])
    return all([name in graph_names for name in konect_names])


def valid_orderings(orders):
    settings = config.settings

    orderings = settings['orderings']
    all_orders = set([k for k in orderings.keys()])
    return all([o in all_orders for o in orders])


def valid_edge_orderings(orders):
    settings = config.settings
    orderings = settings['edge_orderings']
    all_orders = set([k for k in orderings.keys()])
    return all([o in all_orders for o in orders])


def get_all_konect_names():
    settings = config.settings

    datasets_json_path = settings['datasets_json_path']
    with open(datasets_json_path) as f:
        datasets = json.load(f)
    # verify that all graph names exist in datasets.json

    return [d['name'] for d in datasets['datasets']]


def delete_all_rows(table):
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path, timeout=config.settings['sqlite3']['timeout'])
    c = conn.cursor()
    c.execute(f'DELETE FROM {table};', )
    logging.info(f'{c.rowcount} rows deleted from the {table}')

    # commit the changes to db
    conn.commit()
    # close the connection
    conn.close()
    return


def get_all_rows(table):
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path, timeout=config.settings['sqlite3']['timeout'])
    conn.row_factory = sqlite3.Row
    sql = f"select * from {table}"

    cursor = conn.execute(sql)
    rows = cursor.fetchall()
    return rows


def get_datasets(graph_names=[]):
    settings = config.settings

    dataframes_dir = settings['dataframes_dir']
    datasets_json_path = settings['datasets_json_path']
    with open(datasets_json_path) as f:
        datasets = json.load(f)
    # print(json.dumps(data_json, indent=4, sort_keys=True))

    if graph_names:
        return [e for e in datasets['datasets'] if e['name'] in graph_names]
    else:
        return datasets['datasets']


def init_logger(log_file_name):
    fmt = config.settings['logging']['log_format']

    formatter = logging.Formatter(fmt)

    logging.basicConfig(filename=log_file_name, encoding='utf-8',
                        level=logging.DEBUG, format=fmt, filemode='a')


def get_query_vals_str(n_vals):
    return '(' + \
           ','.join(['?'] * n_vals) + \
           ')'


def set_n(graph_name, val):
    return single_val_numeric_set('n', 'n_m', graph_name, val)


def set_m(graph_name, val):
    return single_val_numeric_set('m', 'n_m', graph_name, val)


def set_n_m(graph_name, n, m):
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path, timeout=config.settings['sqlite3']['timeout'])
    cursor = conn.cursor()
    query = f"insert or replace into n_m (graph_name, n, m) values (?,?,?)"
    cursor.execute(query, [graph_name, n, m])
    conn.commit()
    conn.close()


def single_val_numeric_set(col_name, table_name, graph_name, val):
    """
    Find the graph row in the corresponding table and update its column value
    :param col_name:
    :param table_name:
    :param graph_name:
    :param val:
    :return:
    """
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path, timeout=config.settings['sqlite3']['timeout'])
    cursor = conn.cursor()
    query = f"update {table_name} set {col_name} = ? where graph_name = ?"
    cursor.execute(query, [val, graph_name])
    conn.commit()
    conn.close()
    return


def get_directed(graph_name):
    """
    Query the statistics table to check whether a graph is directed or not
    :param graph_name:
    :return:
    """
    return single_val_get('directed', 'statistics', graph_name)


def get_volume(graph_name):
    """
    Get the number of edges in a graph
    :param graph_name:
    :return:
    """
    return single_val_get('volume', 'statistics', graph_name)


def get_size(graph_name):
    """
    Get the number of vertices in a graph
    :param graph_name:
    :return:
    """
    return single_val_get('size', 'statistics', graph_name)


def get_n_m(graph_name):
    # TODO get size, volume instead of n, m -- n, m might be empty
    n = get_n(graph_name)
    m = get_m(graph_name)
    if not n:
        num_nodes = single_val_get('size', 'statistics', graph_name)
    else:
        num_nodes = n
    if not m:
        num_edges = single_val_get('volume', 'statistics', graph_name)
    else:
        num_edges = m
    return num_nodes, num_edges


def get_n(graph_name):
    return single_val_get('n', 'n_m', graph_name)


def get_n_vertices(graph_name):
    return single_val_get('size', 'statistics', graph_name)


def get_n_edges(graph_name):
    return single_val_get('volume', 'statistics', graph_name)


def get_pr_struct_size(graph_name):
    return single_val_get('pr_struct_size', 'statistics', graph_name)


def get_category(graph_name):
    return single_val_get('category', 'metadata', graph_name)


def get_m(graph_name):
    return single_val_get('m', 'n_m', graph_name)


def get_cache_stats():
    """
        linux-specific, assuming single-threaded execution
    """

    def convert_to_bytes(s):
        if 'K' in s:
            return int(s.replace('K', '')) * 1024
        elif 'M' in s:
            return float(s.replace('M', '')) * 1024 * 1024

    args = ['lscpu', '-C']
    try:
        res = subprocess.check_output(args)
    except subprocess.CalledProcessError:
        return {}
    # parse the output of the simplify program to get the number of vertices and edges
    # in the compressed, simplified representation
    lines = [l.split() for l in res.decode('ascii').split('\n')]
    df = pd.DataFrame(columns=lines[0], data=lines[1:])
    if 'COHERENCY-SIZE' in df.columns:
        line_size = int(df['COHERENCY-SIZE'].unique()[0])
    else:
        line_size = 64
    l1d_size = df.loc[df['NAME'] == 'L1d']['ONE-SIZE'].values[0]
    l2_size = df.loc[df['NAME'] == 'L2']['ONE-SIZE'].values[0]
    l3_size = df.loc[df['NAME'] == 'L3']['ONE-SIZE'].values[0]

    d = {
        'line_size': line_size,
        'l1d_size': convert_to_bytes(l1d_size),
        'l2_size': convert_to_bytes(l2_size),
        'l3_size': convert_to_bytes(l3_size),
    }
    return d


def get_size_in_memory(n_vertices, n_edges):
    """
    back-of-the-envelope calculation to compute the size in memory (in MB) of the
    structs used in PR experiments

    Vertices:
        3 std::vector<double> of size n_vertices; where size_of(double) == 8 bytes
    Edges:
        1 std::vector<Edge> of size n_edges; where:
        struct Edge {
            uint32_t source;        # size_of(uint32_t) == 4 bytes
            uint32_t destination;   # size_of(uint32_t) == 4 bytes
            uint64_t index;         # size_of(uint64_t) == 8 bytes
        }
    """
    size_bytes = n_edges * 16 + 3 * n_vertices * 8

    return size_bytes


def convert_size(size_bytes):
    if size_bytes == 0:
        return "0B"
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)
    return s, size_name[i]


def valid_pr(rows):
    """
    Verify that our graph PR computation is equivalent to igraph
    """
    for r in rows:
        graph_name = r['graph_name']
        n = get_n(graph_name)
        directed = get_directed(graph_name)
        settings = config.settings
        graph_dir = os.path.join(settings['graphs_dir'], graph_name)
        graph_path = os.path.join(
            graph_dir, f"{settings['compressed_el_file_name']}.net")
        pr_path = os.path.join(graph_dir, 'pr')

        pr = np.loadtxt(pr_path).astype(np.float64)

        # use igraph to compute pagerank
        igraph_pr = np.array(igraph.Graph.Read_Edgelist(
            graph_path, directed=directed).pagerank())

        if directed:
            nx_graph = nx.DiGraph
        else:
            nx_graph = nx.Graph

        g = nx.read_edgelist(graph_path, create_using=nx_graph)
        nx_pr_dict = nx.pagerank(nx.read_edgelist(
            graph_path, create_using=nx_graph))

        nx_pr = np.zeros(n).astype(np.float64)
        for k, v in nx_pr_dict.items():
            nx_pr[int(k)] = v

        # print(f"{np.allclose(pr, igraph_pr)=}")
        # print(f"{np.allclose(pr, nx_pr)=}")

        valid_pagerank = np.allclose(pr, igraph_pr, rtol=0, atol=1e-4) and \
            np.allclose(pr, nx_pr, rtol=0, atol=1e-4)

        # for i in range(n):
        #     print(f"{pr[i]} {igraph_pr[i]} {nx_pr[i]}")

        return valid_pagerank


def get_critical_depth():
    n_threads = psutil.cpu_count()
    n_quads = 1
    critical_depth = 0
    while n_quads < n_threads:
        critical_depth += 1
        n_quads = n_quads * 4
    return critical_depth


def next_largest_multiple(n, d):
    """

    Args:
        n:
        d:

    Returns:

    """
    multiple = np.power(2, d)
    # print(multiple)
    return int((n + multiple - 1) / multiple) * multiple


def get_unimputed_features(graph_names):
    """
    Given a set of graph names, return a set of Konect statistics that have been
    computed for all those graphs (i.e. non-null values)
    Args:
        graph_names:

    Returns: a sorted (in order of ascending PageRank struct size) list of graph names to download
    that contain a substantial (exact definition of 'substantial' tbd) number of features with
    non-null values

    """
    conn = connect()
    conn.row_factory = sqlite3.Row
    seq = ','.join(['?'] * len(graph_names))
    graph_name_seq = ','.join([f"\'{g}\'" for g in graph_names])
    sql = f"select * from statistics where graph_name in ({seq})"

    df = pd \
        .read_sql(sql, conn, params=graph_names) \
        .sort_values(by='pr_struct_size', ascending=True)

    zero_cols = df.columns[df.apply(lambda x: np.all(x == 0))]
    null_cols = df.columns[df.apply(lambda x: np.all(x.isnull()))]
    valid_cols = set(df.columns) \
        .difference(
        set(null_cols).union(set(zero_cols))
    )

    # get the number of zero entries per column
    non_zero_val_dict = df \
        .apply(lambda x: np.count_nonzero(x)) \
        .sort_values(ascending=False) \
        .to_dict()

    # filter for features for which we have at least min_n_data_samples
    min_n_data_samples = config.settings['modelling']['min_n_data_samples']
    non_zero_val_dict = {
        k: v for k, v in non_zero_val_dict.items() if v >= min_n_data_samples}

    # ignore all graphs that have zero entries for the filtered features
    valid_cols_str = ',\n'.join(non_zero_val_dict.keys())

    df = df[non_zero_val_dict.keys()]
    res = df.sort_values(by='pr_struct_size', ascending=True)[
        'graph_name'].values

    sql_query = f"select {valid_cols_str} from statistics where graph_name in ({graph_name_seq}) order by pr_struct_size"
    # print(sql_query)
    # cursor = conn.execute(sql, graph_names)
    # rows = cursor.fetchall()

    return res


def save_ground_truth_pr(edgelist_path, graph_name):
    """Run GapBS' Pagerank on the graph specified
    Sideffect - Saves the Ground-Truth PageRank values in a binary file in 
    the graph's directory (for future verification)

    Args:
        edgelist_path (str): path to the graph's text edge list
        graph_name (str): the graph's (konect) name/label
    """
    settings = config.settings
    pr_exec = settings['pr_executable']
    pr_filename = settings['pagerank_file_name']
    graphs_dir = settings['graphs_dir']
    graph_dir = os.path.join(graphs_dir, graph_name)
    pr_path = os.path.join(graph_dir, pr_filename)

    args = [
        pr_exec,
        '-f', edgelist_path,
        '-c', pr_path

    ]
    logging.info(f"Executing: {' '.join(args)}")
    res = subprocess.check_output(args)

    return


def copy_order_file_to_binary(n, input_path, output_path):
    """Use a cpp utility to convert vertex orderings saved as text files
    to a binary file (using boost::serialization archives)

    Args:
        n (int): number of vertices in the graph
        input_path (string): path to vertex ordering text file
        output_path (string): where to save the binary file
    """
    convert_map_to_bin_executable = config.settings['convert_map_to_bin_executable']
    args = [
        convert_map_to_bin_executable,
        '-n', str(n),
        '-i', input_path,
        '-o', output_path,
    ]
    subprocess.check_output(args)
    return
