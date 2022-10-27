import logging
import os
import sqlite3
from datetime import datetime

import pandas as pd

from konect_scraper import config, column_names, scrape_konect_stats
from konect_scraper.util import create_data_dirs_if_not_exists, create_sql_table, init_logger, \
    get_all_rows, update_table_schema, create_pr_expt_table, create_log_dir_if_not_exists


def update_ns_ms():
    settings = config.settings
    graphs_dir = settings['graphs_dir']
    sub_dirs = os.listdir(graphs_dir)
    for graph_name in sub_dirs:
        sb_ord_path = os.path.join(graphs_dir, graph_name, "sb")
        with open(sb_ord_path, 'r') as f:
            n = int(f.readline().strip())
            m = int(f.readline().strip())
            print(graph_name, n, m)


def main():
    """
    This script:
    - Initializes a sqlite3 db with the following tables:
        1. konect     - stores the download urls for all graphs on konect
        2. preproc    - the preprocessing and reordering times for respective graphs
        3. metadata   - metadata (format, category, data source) for each graph
        4. statistics - numerical values that characterizes a network (from http://konect.cc/statistics/)

        parser.add_argument('-i', '--initialize',
                        action=argparse.BooleanOptionalAction, required=True,
                        help='Whether to initialize the sqlite database or not. '
                             'If true, existing sqlite db in the path specified by '
                             'config.settings.sqlite3.sqlite3_db_path will be deleted')
    """
    config.init()
    column_names.init()
    settings = config.settings

    log_dir = settings['logging']['log_dir']
    create_log_dir_if_not_exists()
    curr_time = datetime.now().strftime("%H_%M_%d_%m_%Y")
    log_path = f"init_db_{curr_time}"
    log_file_name = os.path.join(log_dir, log_path + '.' + 'log')
    init_logger(log_file_name)

    logging.getLogger("requests").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    logging.getLogger("charset_normalizer").setLevel(logging.WARNING)

    db_path = settings['sqlite3']['sqlite3_db_path']
    repo_root = settings['repo_root']

    # remove any existing sqlite3 db TODO add argument for conditional first exection
    # delete_graphs_db()

    # create data dirs (if not exists)
    create_data_dirs_if_not_exists()

    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path)

    tables = ['metadata', 'statistics', 'preproc', 'konect', 'n_m', 'directed',
              'undirected', 'bipartite']
    columns = [
        column_names.meta_col_names,
        column_names.stat_col_names,
        column_names.preproc_col_names,
        column_names.konect_col_names,
        column_names.n_m_col_names,
    ] + [column_names.graph_dataframe_col_names] * 3

    for table, cols in zip(tables, columns):
        create_sql_table(conn, table, cols)
        update_table_schema(table, cols)
    db_path = config.settings['sqlite3']['sqlite3_db_path']

    graph_table_names = ['directed', 'undirected', 'bipartite']
    for table_name in graph_table_names:
        df = pd.read_csv(os.path.join(settings['dataframes_dir'], f"{table_name}.csv"))
        conn = sqlite3.connect(db_path)
        df.to_sql(table_name, con=conn, if_exists='replace')

    create_pr_expt_table()
    # scrape_konect_stats.fill_konect_table() TODO add argument for conditional first exection

    # remove any existing rows scraped from konect
    # for table in ['metadata', 'statistics', 'preproc']:
    #     delete_all_rows(table)

    # optionally, write konect to csv
    conn = sqlite3.connect(db_path)
    # db_df = pd.read_sql_query("select * from konect", conn)
    # db_df.to_csv('./konect.csv', index=False)

    # on first execution populate konect table from konect.csv
    df = pd.read_csv(
        os.path.join(
            config.settings['dataframes_dir'],
            'konect.csv'
        )
    )
    df.to_sql('konect', con=conn, if_exists='replace')
    # get all rows from konect table and scrape konect for those rows
    rows = get_all_rows("konect")
    scrape_konect_stats.main(rows)

    # add existing n, m values to now reset db
    # update_ns_ms()
    return


if __name__ == '__main__':
    main()
