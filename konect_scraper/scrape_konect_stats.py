#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ast

import requests
import pandas as pd
import numpy as np
import sqlite3
import os.path
from sqlalchemy.dialects.sqlite import insert

from konect_scraper import column_names, config


def insert_on_duplicate(table, conn, keys, data_iter):
    insert_stmt = insert(table.table).values(list(data_iter))
    on_duplicate_key_stmt = insert_stmt.on_duplicate_key_update(insert_stmt.inserted)
    conn.execute(on_duplicate_key_stmt)


def get_feature_labels(url, meta_df, stats_df, graph_name):
    r = requests.get(url)
    tables = pd.read_html(url)

    meta_df, directed = parse_meta_table(tables[0], meta_df, graph_name)
    stats_df = parse_stats_table(tables[1], stats_df, graph_name, directed)

    return meta_df, stats_df


def parse_numeric(s):
    try:
        return ast.literal_eval(s)
    except:
        if '×' in s:  # string with scientific notation
            l = s.split('×')
            a = ast.literal_eval(''.join(l[0].replace(u'\xa0', u'')))
            b = l[1].strip()
            if '−' in b:
                exp = b.replace('10', '')
                exp = exp.replace('−', '')
                return a * np.float_power(10, -int(exp))
            else:
                exp = int(b.replace('10', ''))
                return a * np.float_power(10, exp)
            return 0
        elif '−' in s:
            s = s.replace('−', '')
            return -ast.literal_eval(s.replace(u'\xa0', u''))
        elif '+' in s:
            s = s.replace('+', '')
            return ast.literal_eval(s.replace(u'\xa0', u''))

        else:
            return ast.literal_eval(''.join(s.replace(u'\xa0', u'')))


def convert(s):
    c = s.lower() \
        .replace(' ', '_') \
        .replace('-', '_') \
        .replace('/', '_')

    # sqlite doesn't like column names starting with numeric - move it to end of string
    if c.split('_')[0].isnumeric():
        pre = c.split('_')[0]
        c = c.replace(pre + '_', '')
        return c + '_' + pre

    if c == 'join':
        return 'network_join'

    return c


def parse_stats_table(table, df, graph_name, directed):
    data = {}
    data['graph_name'] = graph_name
    s1 = set()
    for key, value in zip(table[0].values, table[2].values):
        literal = parse_numeric(value)
        s1.add(convert(key))
        data[convert(key)] = literal

    data['directed'] = directed

    s2 = set(list(column_names.stat_col_names.keys()))

    missing_stat_cols = s1.difference(s2)
    if missing_stat_cols:
        print(f"{missing_stat_cols=}!")

    # get the difference between the expected columns and what is available on konect

    df = pd.concat([df, pd.DataFrame.from_records([data])], ignore_index=True)
    return df


def parse_meta_table(table, df, graph_name):
    data = {}
    data['graph_name'] = graph_name
    directed = 0
    for key, value in zip(table[0].values, table[2].values):
        if key == 'Network format':
            if 'undirected' in value:
                directed = 0
            else:
                directed = 1
        data[convert(key)] = value

    df = pd.concat([df, pd.DataFrame.from_records([data])], ignore_index=True)

    return df, directed


def aggregate_feature_labels():
    return


def write_to_sqlite3(df, table_name, conn):
    columns_str = ','.join(df.columns)

    values_str = '('
    for c in df.columns[:-1]:
        values_str += '?,'
    values_str += '?)'

    query = f''' insert or replace into {table_name} ({columns_str}) values {values_str} '''
    # print(query)
    # # print(df.replace(np.nan, '').to_records(index=False))
    # for row in df.to_records(index=False):
    #     print(list(row))
    #     conn.execute(query, list(row))
    conn.executemany(query, df.to_records(index=False))
    conn.commit()


def main(datasets):
    settings = config.settings

    dataframes_dir = settings['dataframes_dir']

    stats_df = pd.DataFrame(columns=list(column_names.stat_col_names.keys()))
    meta_df = pd.DataFrame(columns=column_names.meta_col_names.keys())
    # meta_df.set_index('graph_name', inplace=True)
    # stats_df.set_index('graph_name', inplace=True)
    # get all konect url for datasets in datasets.json
    for dataset in datasets:
        url = dataset['konect-url']
        graph_name = dataset['name']
        # aggregate all possible feature labels available for datasets on konect
        meta_df, stats_df = get_feature_labels(url, meta_df, stats_df, graph_name)

    # populate a dataframe with feature values for datasets
    meta_df_path = os.path.join(dataframes_dir, "meta.csv")
    stats_df_path = os.path.join(dataframes_dir, "stats.csv")

    meta_df.to_csv(meta_df_path)
    stats_df.to_csv(stats_df_path)

    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path)

    write_to_sqlite3(meta_df, 'metadata', conn)
    write_to_sqlite3(stats_df, 'statistics', conn)

    # write_to_sqlite3(meta_df, 'metadata')

    # persist meta info and stats dataframes to sqlite3 database
    conn.close()
    return


if __name__ == '__main__':
    # column_names.init()
    # config.init()
    main()
