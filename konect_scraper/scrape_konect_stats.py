#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ast
import logging
import os.path
import sqlite3
import unicodedata
import urllib.request

import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from sqlalchemy.dialects.sqlite import insert

from konect_scraper import column_names, config
from konect_scraper.sql import read_precomputed_stats
from konect_scraper.stats import verify_stat
from konect_scraper.util import get_query_vals_str, cast_np_dtypes, get_size_in_memory


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


def get_data_url(konect_url):
    def is_unavailable(soup):
        # find metadata table
        h2 = soup.find('h2', text='Metadata')
        table = h2.find_next_sibling('table')
        rows = table.find('tr')
        for row in rows:
            cells = row.find_all("td")
            rn = cells[0].get_text()
            if "Availability" in rn:
                for line in rn.split('\n'):
                    if 'Availability' in line:
                        return line.replace('Availability', '') == 'Dataset is not available for download'

    tables = pd.read_html(konect_url)
    html_page = urllib.request.urlopen(konect_url)
    soup = BeautifulSoup(html_page, "html.parser")
    tables = soup.find(lambda tag: tag.name == 'table')
    hrefs = tables.findAll('a', href=True)
    konect_files_url = "http://konect.cc/files"
    files_relative_path = '../../files'

    if is_unavailable(soup):
        url = "none"
    else:
        for a in hrefs:
            if a.text == 'Dataset is available for download':
                href = a['href']
                if href.startswith(files_relative_path):
                    url = os.path.join(konect_files_url + href.replace(files_relative_path, ''))
                else:
                    url = href

    return url


def fill_konect_table():
    settings = config.settings
    all_networks_url = settings['all_networks_url']
    tables = pd.read_html(all_networks_url)
    html_page = urllib.request.urlopen(all_networks_url)
    soup = BeautifulSoup(html_page, "html.parser")
    urls = []
    i = 1
    links = soup.findAll('a')[2:-2]
    assert (len(links) - 1) / 2 == tables[0].shape[0]
    for link in links:
        if i % 2 == 0:
            urls.append(all_networks_url + link.get('href'))
        i += 1
    rows = []
    for url, (i, row) in zip(urls, tables[0][1:].iterrows()):
        internal_name = os.path.basename(os.path.normpath(url))
        logging.info(f"{internal_name} : {url}")

        internal_code = row[0]
        name = row[1]
        n = int(row[3])
        m = int(row[4])

        rows.append((
            # graph_name, name, code, n, m, konect_url, data_url
            internal_name, name, internal_code, n, m, url, get_data_url(url)
        ))
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    query_vals_str = get_query_vals_str(len(column_names.konect_col_names))
    query = f"insert into konect values {query_vals_str}"
    cursor.executemany(query, rows)
    conn.commit()
    conn.close()
    return


def parse_numeric(s):
    utimes = '\u00D7'
    uminus = '\u2212'
    uplus = '\u002B'

    def parse_scientific_notation_string(string):
        line = string.split(utimes)
        coef = ''.join(line[0].replace(u'\xa0', u''))

        if uminus in coef:
            coef = ast.literal_eval(coef.replace(uminus, ''))
            coef = -int(coef)
        else:
            coef = ast.literal_eval(coef)

        base = line[1].strip()
        exp = base.replace('10', '', 1)
        if uminus in exp:
            exp = exp.replace(uminus, '')
            exp = -int(exp)
        else:
            exp = int(exp)

        return coef * np.float_power(10, exp)

    s = unicodedata.normalize('NFC', s)
    try:
        return ast.literal_eval(s)
    except:
        if utimes in s:  # string with scientific notation
            return parse_scientific_notation_string(s)
        elif uminus in s:
            s = s.replace(uminus, '')
            return -ast.literal_eval(s.replace(u'\xa0', u''))
        elif uplus in s:
            s = s.replace(uplus, '')
            return ast.literal_eval(s.replace(u'\xa0', u''))
        elif 'Inf' in s:
            return np.inf
        elif '-Inf' in s:
            return np.NINF
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

def update_to_sqlite3(stats, table_name, conn):
    sql = f'UPDATE {table_name} SET '

    sql += ', '.join([f'{k} = ?' for k in stats.keys() if k != 'graph_name'])
    sql += ' WHERE graph_name = ?'
    vals = stats.values()
    cur = conn.cursor()
    cur.execute(sql, list(stats.values()))
    conn.commit()

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


def main(rows):
    settings = config.settings

    dataframes_dir = settings['dataframes_dir']

    stats_df = pd.DataFrame(columns=list(column_names.stat_col_names.keys()))
    meta_df = pd.DataFrame(columns=column_names.meta_col_names.keys())
    # meta_df.set_index('graph_name', inplace=True)
    # stats_df.set_index('graph_name', inplace=True)
    # get all konect url for datasets in datasets.json
    logging.info(f"Scraping {len(rows)} graphs from konect.cc/networks..")
    # for row in rows:
    #
    #     url = row['konect_url']
    #     graph_name = row['graph_name']
    #     logging.info(f"Scraping {graph_name} at {url}")
    #     # aggregate all possible feature labels available for datasets on konect
    #     meta_df, stats_df = get_feature_labels(url, meta_df, stats_df, graph_name)
    #
    # # populate a dataframe with feature values for datasets
    meta_df_path = os.path.join(dataframes_dir, "meta.csv")
    stats_df_path = os.path.join(dataframes_dir, "stats.csv")
    #
    # meta_df.to_csv(meta_df_path, index=False)
    # stats_df.to_csv(stats_df_path, index=False)

    meta_df = pd.read_csv(meta_df_path, index_col=0)
    stats_df = pd.read_csv(stats_df_path, index_col=0)

    # add precomputed (e.g. slashburn) stats if they already exist
    precomputed_stats_df = read_precomputed_stats()
    # return
    # use the (uncompressed) size and volume of the graphs to compute (roughly)
    # calculate the sizes of the PageRank computation structs
    stats_df['pr_struct_size'] = stats_df.apply(
        lambda row: get_size_in_memory(row['size'], row['volume']),
        axis=1
    ).astype(np.int64)

    meta_dtypes = {
        k: column_names.sql_to_np_dtypes[column_names.meta_col_names[k]] for k in
        meta_df.columns
    }

    stats_dtypes = {
        k: column_names.sql_to_np_dtypes[column_names.stat_col_names[k]] for k in
        stats_df.columns
    }

    cast_np_dtypes(meta_df, meta_dtypes)
    cast_np_dtypes(stats_df, stats_dtypes)

    meta_df = meta_df.astype(meta_dtypes)
    for col in stats_df.columns:
        stats_df[col] = stats_df[col].astype(stats_dtypes[col])
    row = stats_df.loc[stats_df['graph_name'] == 'zhishi-baidu-internallink']
    # some konect stats are either not computed or incorrect - recalculate them
    stat_cols_to_verify = [
        # 'fill',
    ]

    for col in stat_cols_to_verify:
        logging.info(f"Verifying {col}..")
        stat_col = []
        for graph_name in stats_df['graph_name'].values:
            stat_col.append(verify_stat(col, graph_name))
        print(stat_col)

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
