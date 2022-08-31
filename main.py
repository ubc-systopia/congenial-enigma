import argparse
from argparse import RawTextHelpFormatter
import konect_scraper.config as config
import sqlite3
import konect_scraper.column_names as column_names
from konect_scraper import scrape_konect_stats, download_and_extract
from konect_scraper.util import \
    create_sql_table, delete_graphs_db, verify_graphs_in_json, get_datasets, create_data_dirs_if_not_exists


def main(args):
    init = args.initialize
    plot = args.plot
    preprocess = args.preprocess
    download = args.download
    konect_internal_names = args.graph_names

    if konect_internal_names:
        assert verify_graphs_in_json(konect_internal_names)
        datasets = get_datasets(konect_internal_names)
    else:
        datasets = get_datasets()

    if init:
        # create data dirs (if not exists)
        create_data_dirs_if_not_exists()

        # remove any existing sqlite3 db
        delete_graphs_db()

        db_path = config.settings['sqlite3']['sqlite3_db_path']
        conn = sqlite3.connect(db_path)

        create_sql_table(conn, 'metadata', column_names.meta_col_names)
        create_sql_table(conn, 'statistics', column_names.stat_col_names)
        create_sql_table(conn, 'preproc', column_names.preproc_col_names)

        scrape_konect_stats.main(datasets)

    if download:
        download_and_extract.main(datasets)

    if preprocess:
        # verify that the requested ordering to compute are supported
        print(preprocess)

    return


if __name__ == '__main__':
    config.init()
    column_names.init()
    argparse_desc = """
    Konect Scraper: an automatic scraper that:
    1.  Downloads either:
        1.1.    all graphs in datasets.json or
        1.2.    a single graph specified by the konect `Internal Name` (must exist in datasets.json)
        Specify single graph using -g argument (all graphs will be downloaded by default).\n
    2.  Scrapes all available stats for the graph(s) from konect and populates a sqlite3 database\n
    3.  Computes orderings using dbg, rabbit, cuthill-mckee, and slashburn
        3.1.    The runtime of each reordering operation is persisted\n
    4.  Optionally, plot the adjacency matrix of each isomorphism. (spy plots should only be plotted for relatively small graphs - < 10,000 vertices). 
        4.1.    Specify optional plotting using -p argument.\n
    5.  Iterates over the edges of the isomorphisms using either Row/Column-major or Hilbert order to compute the PageRank (PR)
        5.1.    The runtime of each experiment is persisted.
    """
    parser = argparse.ArgumentParser(description=argparse_desc, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--initialize',
                        action=argparse.BooleanOptionalAction, required=True,
                        help='Whether to initialize the sqlite database or not. '
                             'If true, existing sqlite db in the path specified by '
                             'config.settings.sqlite3.sqlite3_db_path will be deleted')

    parser.add_argument('-g', '--graph-names', nargs='+',
                        help='If specified, only download and scrape these graphs '
                             'from dataset (must exist in datasets.json). '
                             'Otherwise, download ALL graphs in datasets.json.', )

    parser.add_argument('-d', '--download',
                        action=argparse.BooleanOptionalAction, required=True,
                        help='Whether to download the graph from konect or not.')

    parser.add_argument('-p', '--preprocess', nargs='+',
                        help='A list of orderings to compute, if any.\n'
                             'Possible Orders: {\n'
                             '  `rnd`: Random, \n'
                             '  `rbt`: Rabbit, \n'
                             '  `sb`:  Slashburn, \n'
                             '  `cm`:  Cuthill-McKee, \n'
                             '  `srt`: Descending Degree Sort, \n'
                             '  `hc`:  HubCluster, \n'
                             '  `hs`:  HubSort, \n'
                             '  `dbg`: Degree Based Grouping,\n'
                             '}')

    parser.add_argument('-l', '--plot', action=argparse.BooleanOptionalAction,
                        help='Whether to plot the adjacency matrices or not.', required=True)

    main(parser.parse_args())

