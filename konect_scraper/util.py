import sqlite3
import config
import os
import json

def __init__(self):
    config.init()


def connect():
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path)
    return conn


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
    print(f"Deleting {db_path}..")
    os.remove(db_path)

def verify_graphs_in_json(konect_names):
    settings = config.settings

    datasets_json_path = settings['datasets_json_path']
    with open(datasets_json_path) as f:
        datasets = json.load(f)
    # verify that all graph names exist in datasets.json

    graph_names = set([d['name'] for d in datasets['datasets']])
    return all([name in graph_names for name in konect_names])


def get_all_konect_names():
    settings = config.settings

    datasets_json_path = settings['datasets_json_path']
    with open(datasets_json_path) as f:
        datasets = json.load(f)
    # verify that all graph names exist in datasets.json

    return [d['name'] for d in datasets['datasets']]

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


def is_directed(graph_name):
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
