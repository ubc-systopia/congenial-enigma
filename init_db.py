import logging
import os
import sqlite3
from datetime import datetime

from konect_scraper import config, column_names, scrape_konect_stats
from konect_scraper.util import delete_graphs_db, create_data_dirs_if_not_exists, create_sql_table, init_logger, \
    get_datasets, get_all_rows, delete_all_rows, column_exists, add_column_if_not_exists


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

    tables = ['metadata', 'statistics', 'preproc', 'konect']
    columns = [
        column_names.meta_col_names,
        column_names.stat_col_names,
        column_names.preproc_col_names,
        column_names.konect_col_names,
    ]


    for table, cols in zip(tables, columns):
        create_sql_table(conn, table, cols)


    for col in column_names.meta_col_names:
        dtype = column_names.meta_col_names[col]
        add_column_if_not_exists(col, 'metadata', dtype)

    for col in column_names.stat_col_names:
        dtype = column_names.stat_col_names[col]
        add_column_if_not_exists(col, 'statistics', dtype)

    # scrape_konect_stats.fill_konect_table() TODO add argument for conditional first exection

    # remove any existing rows scraped from konect
    for table in ['metadata', 'statistics', 'preproc']:
        delete_all_rows(table)

    # get all rows from konect table and scrape konect for those rows
    rows = get_all_rows("konect")
    scrape_konect_stats.main(rows)
    return


if __name__ == '__main__':
    main()
