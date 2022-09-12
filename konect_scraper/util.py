import sqlite3
import os
import json
from pathlib import Path
import logging
import numpy as np
from konect_scraper import config
import pandas as pd

from konect_scraper.sql import connect


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
            arr = np \
                .nan_to_num(df[col].values) \
                .astype(np.ulonglong)
            df[col] = pd.to_numeric(arr, errors='coerce').astype(np.ulonglong)
            # return
        
        dtype = dtypes[col]
        if dtype == 'Int64':
            arr = np \
                .nan_to_num(df[col].values) \
                .astype(np.int64)
            df[col] = pd.to_numeric(arr, errors='coerce').astype(np.int64)
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
    conn.close()

    return


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


def single_val_numeric_get(col_name, table_name, graph_name):
    conn = connect()
    cursor = conn.cursor()
    query = f"select {col_name} from {table_name} where graph_name = ?"
    cursor.execute(query, [graph_name])
    res = cursor.fetchone()[0]
    conn.close()
    return res


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


def create_log_dir_if_not_exists():
    data_dirs = [config.settings['logging']['log_dir']]

    for data_dir in data_dirs:
        Path(data_dir).mkdir(parents=True, exist_ok=True)


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


def get_all_konect_names():
    settings = config.settings

    datasets_json_path = settings['datasets_json_path']
    with open(datasets_json_path) as f:
        datasets = json.load(f)
    # verify that all graph names exist in datasets.json

    return [d['name'] for d in datasets['datasets']]


def delete_all_rows(table):
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute(f'DELETE FROM {table};', );
    logging.info(f'{c.rowcount} rows deleted from the {table}')

    # commit the changes to db
    conn.commit()
    # close the connection
    conn.close()
    return


def get_all_rows(table):
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path)
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

    logging.basicConfig(filename=log_file_name, encoding='utf-8', level=logging.DEBUG, format=fmt, filemode='a')


def get_query_vals_str(n_vals):
    return '(' + \
           ','.join(['?'] * n_vals) + \
           ')'


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
    conn = sqlite3.connect(db_path)
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
    return single_val_numeric_get('directed', 'statistics', graph_name)


def get_volume(graph_name):
    """
    Get the number of edges in a graph
    :param graph_name:
    :return:
    """
    return single_val_numeric_get('volume', 'statistics', graph_name)


def get_size(graph_name):
    """
    Get the number of vertices in a graph
    :param graph_name:
    :return:
    """
    return single_val_numeric_get('size', 'statistics', graph_name)


def get_n(graph_name):
    return single_val_numeric_get('n', 'statistics', graph_name)


def get_m(graph_name):
    return single_val_numeric_get('m', 'statistics', graph_name)


def set_n(graph_name, val):
    return single_val_numeric_set('n', 'statistics', graph_name, val)


def set_m(graph_name, val):
    return single_val_numeric_set('m', 'statistics', graph_name, val)
