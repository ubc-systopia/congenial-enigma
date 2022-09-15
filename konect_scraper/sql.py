import sqlite3

import pandas as pd

from konect_scraper import config


def connect():
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path)
    return conn


def get_all_graphs_in_categories(categories):
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    seq = ','.join(['?'] * len(categories))
    sql = f"select * from metadata where category in ({seq})"
    cursor = conn.execute(sql, categories)
    rows = cursor.fetchall()
    return rows


def get_all_unipartite_graphs():
    conn = connect()
    conn.row_factory = sqlite3.Row
    sql = f"select * from metadata where network_format like 'Unipartite%'"
    cursor = conn.execute(sql)
    return cursor.fetchall()

def get_all_unipartite_directed_graphs():
    conn = connect()
    conn.row_factory = sqlite3.Row
    sql = f"select * from metadata where network_format like 'Unipartite, directed%'"
    cursor = conn.execute(sql)
    return cursor.fetchall()
def row_as_dict(r):
    return {k: r[k] for k in r.keys()}


def get_all_downloadable_graphs(graph_names):
    conn = connect()
    conn.row_factory = sqlite3.Row
    seq = ','.join(['?'] * len(graph_names))

    sql = f"select * from konect where graph_name in ({seq}) and data_url != 'none'"

    cursor = conn.execute(sql, graph_names)
    return cursor.fetchall()

def read_table_as_table(table):

    # Create your connection.
    cnx = sqlite3.connect(config.settings['sqlite3']['sqlite3_db_path'])

    df = pd.read_sql_query(f"SELECT * FROM {table}", cnx)
    return df


def get_all_rows_by_graph_names(table, graph_names):
    conn = connect()
    conn.row_factory = sqlite3.Row
    seq = ','.join(['?'] * len(graph_names))

    sql = f"select * from {table} where graph_name in ({seq})"

    cursor = conn.execute(sql, graph_names)
    return cursor.fetchall()


def get_all_graphs_by_graph_names_where_stats_between(stats, mins, maxs, graph_names):
    assert len(stats) == len(mins) == len(maxs)
    conn = connect()
    conn.row_factory = sqlite3.Row
    seq = ','.join(['?'] * len(graph_names))

    sql = f"select * from statistics where "
    for col, mn, mx in zip(stats, mins, maxs):
        # sql += f"({col} >= {mn} and {col} <= {mx}) and "
        sql += f"({col} between {mn} and {mx}) and "

    sql += f"(graph_name in ({seq}))"

    cursor = conn.execute(sql, graph_names)
    return cursor.fetchall()


def get_all_graphs_where_stats_between(stats, mins, maxs):
    assert len(stats) == len(mins) == len(maxs)

    conn = connect()
    conn.row_factory = sqlite3.Row

    sql = f"select * from statistics where "
    for col, mn, mx in zip(stats[:-1], mins[:-1], maxs[:-1]):
        sql += f"{col} between {mn} and {mx} and "
    sql += f"{col} between {mn} and {mx}"
    cursor = conn.execute(sql)
    return cursor.fetchall()


def distinct(column, table):
    conn = connect()
    sql = f"select distinct {column} from {table}"
    cursor = conn.execute(sql)
    rows = cursor.fetchall()
    return rows
