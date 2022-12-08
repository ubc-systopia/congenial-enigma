from ctypes import alignment
import sqlite3

import pandas as pd

from konect_scraper import column_names, config


def connect():
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    timeout = config.settings['sqlite3']['timeout']
    conn = sqlite3.connect(db_path, timeout=timeout)

    return conn


def print_rows(rs):
    """Given an arbitrary list of rows, use the rows' keys to pretty print 
    a table of the rows' contents

    Args:
        rs (list of sqlite3.Row): 
    """
    ks = rs[0].keys()
    width = 30

    heading_str = ""
    for k in ks:
        heading_str += f"{k : <{width }}"

    # heading
    print("-" * width * len(ks))
    print(heading_str)
    print("-" * width * len(ks))

    # rows
    for r in rs:
        r_str = ""
        for k in ks:
            r_str += f"{r[k] : <{width }}"
        print(r_str)

    return


def print_row(r):
    """Given an arbitrary 

    Args:
        r (_type_): _description_
    """
    return


def append_df_to_table(df, table):
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    timeout = config.settings['sqlite3']['timeout']
    conn = sqlite3.connect(db_path, timeout=timeout)

    df.to_sql(table, index=False, con=conn, if_exists='append')

    # col_names_str = ','.join(df.columns)
    # vals_str = ','.join(['?'] * len(df.columns))
    # query=f"insert or replace into {table} ({col_names_str}) values ({vals_str})"
    # conn.executemany(query, df.to_records(index=False))
    conn.commit()


def get_graphs_by_graph_numbers(ns, graph_type):
    """given a min, max graph numbers and a graph type, return all relevant rows

    Args:
        ns (int pair): pair of ints - [min graph number, max graph number)
        graph_type (string): type of graph - one of: directed, undirected, 
                                                     bipartite
    """
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    mn = ns[0]
    mx = ns[1]
    sql = f"select * from {graph_type} where graph_number between {mn} and {mx}"
    cursor = conn.execute(sql)
    rows = cursor.fetchall()

    return rows


def get_all_graphs_in_categories(categories):
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    seq = ','.join(['?'] * len(categories))
    sql = f"select * from metadata where category in ({seq})"
    cursor = conn.execute(sql, categories)
    rows = cursor.fetchall()
    return rows


def is_bipartite(graph_name):
    return "Bipartite" in single_val_get('network_format', 'metadata', graph_name)


def get_all_bipartite_graphs():
    conn = connect()
    conn.row_factory = sqlite3.Row
    sql = f"select * from metadata where network_format like 'Bipartite%'"
    cursor = conn.execute(sql)
    return cursor.fetchall()


def get_all_unipartite_graphs():
    conn = connect()
    conn.row_factory = sqlite3.Row
    sql = f"select * from metadata where network_format like 'Unipartite%'"
    cursor = conn.execute(sql)
    return cursor.fetchall()


def read_precomputed_stats():
    conn = connect()
    additional_stat_cols = [
        'orig_bandwidth',
        'cm_bandwidth',
        'sb_k',
        'par_sb_k',
        'sb_n_iters',
        'par_sb_n_iters',
        'pr_struct_size',
    ]
    try:
        df = read_table_as_table('statistics')
    except:
        return
    return df[additional_stat_cols]


def get_all_unipartite_undirected_graphs():
    conn = connect()
    conn.row_factory = sqlite3.Row
    sql = f"select * from metadata where network_format like 'Unipartite, undirected%'"
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


def get_all_graphs_by_graph_names(graph_names):
    conn = connect()
    conn.row_factory = sqlite3.Row
    seq = ','.join(['?'] * len(graph_names))

    sql = f"select * from konect where "

    sql += f"(graph_name in ({seq}))"

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


def insert_row_if_not_exists(graph_name, table):
    conn = connect()
    cursor = conn.cursor()
    cols = column_names.features_col_names.keys()
    n_cols = len(cols)

    cols = column_names.features_col_names.keys()
    n_cols = len(cols)

    val_str = ','.join(['?'] * n_cols)
    col_str = ','.join(cols)
    sql = f"insert or ignore into {table}(graph_name) VALUES(?)"
    val_str = ','.join(['?'] * n_cols)
    col_str = ','.join(cols)
    sql = f"insert or ignore into {table}({col_str}) VALUES({val_str})"
    res = cursor.execute(sql, [graph_name] + [-1] * (n_cols - 1))
    r = conn.commit()

    return


def single_val_get(col_name, table_name, graph_name):
    conn = connect()
    cursor = conn.cursor()
    query = f"select {col_name} from {table_name} where graph_name = ?"
    cursor.execute(query, [graph_name])
    try:
        res = cursor.fetchone()[0]
    except TypeError:
        return None
    conn.close()
    return res
