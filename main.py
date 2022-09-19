import argparse
from argparse import RawTextHelpFormatter

import numpy as np

import konect_scraper.config as config
import sqlite3
import konect_scraper.column_names as column_names
from konect_scraper import scrape_konect_stats, download_and_extract, reorder, pr_experiments
from konect_scraper.sql import distinct, get_all_graphs_in_categories, get_all_graphs_where_stats_between, \
    get_all_unipartite_graphs, get_all_graphs_by_graph_names_where_stats_between, get_all_rows_by_graph_names, \
    get_all_downloadable_graphs, row_as_dict, get_all_unipartite_directed_graphs
from konect_scraper.util import \
    create_sql_table, delete_graphs_db, verify_graphs_in_json, get_datasets, create_data_dirs_if_not_exists, \
    create_log_dir_if_not_exists, init_logger, valid_orderings, valid_pr, get_category, get_pr_struct_size
from datetime import datetime
import logging
import konect_scraper.plot as plotting
import os
from konect_scraper.config import IOMode


def main(args):
    create_log_dir_if_not_exists()
    plot = args.plot
    orders = args.reorder
    io_modes = args.io_modes
    download = args.download
    run_pr_expts = args.run_pr_expts
    debug = args.debug

    config.settings['debug'] = debug

    log_dir = config.settings['logging']['log_dir']
    curr_time = datetime.now().strftime("%H_%d_%m_%Y")
    log_path = f"{config.settings['app_name']}_{curr_time}"
    log_file_name = os.path.join(log_dir, log_path + '.' + 'log')
    init_logger(log_file_name)

    categories = [
        'Citation network',
        'Online social network'
    ]

    rows = get_all_graphs_in_categories(categories)
    graphs = [r['graph_name'] for r in rows]

    # get all where 50 < size < 100 and 100 < volume < 1000
    rows = get_all_graphs_where_stats_between(
        stats=['size', 'volume'],
        mins=[50, 100],
        maxs=[100, 1000]
    )

    # rows = get_all_unipartite_graphs()
    rows = get_all_unipartite_directed_graphs()
    graph_names = [r['graph_name'] for r in rows]
    print(f"{len(graph_names)} unipartite graphs in dataset.")

    # get all graphs that are 10x-?x as big as the l3 cache
    # 1000 * l3_cache_size
    # np.inf
    # l3_cache_size = config.settings['cpu-info']['cache-sizes']['l3_size']
    # rows = get_all_graphs_by_graph_names_where_stats_between(
    #     stats=['pr_struct_size', ],
    #     mins=[5 * l3_cache_size, ],
    #     maxs=[2**64 - 1, ],
    #     graph_names=graph_names
    # )

    # get all where 50 < size < 100 and 100 < volume < 1000
    rows = get_all_graphs_by_graph_names_where_stats_between(
        stats=['size', ],
        mins=[1000, ],
        maxs=[2000, ],
        graph_names=graph_names
    )

    # a selection of relevant graph_names have been identified,
    # download, reorder, and plot them

    # if io mode is unspecified, use text as default
    if not io_modes:
        io_modes = [IOMode.text]
    else:
        modes = []
        for mode in io_modes:
            match mode:
                case 'binary':
                    modes.append(IOMode.binary)
                case 'text':
                    modes.append(IOMode.text)
                case _:
                    logging.error(f"{mode}: Unsupported IO mode!")
        io_modes = modes
    graph_names = [r['graph_name'] for r in rows]
    rows = get_all_downloadable_graphs(graph_names)[:]

    rows = sorted(rows, key=lambda r: get_pr_struct_size(r['graph_name']), reverse=False)[:]
    # print([r['graph_name'] for r in rows])
    graph_name_start_idx = -2
    graph_name_end_idx =-1
    rows = rows[graph_name_start_idx:graph_name_end_idx]

    for i, row in enumerate(rows):
        print(f"{i : <5} {row['graph_name'] : <40}"
              f"{get_pr_struct_size(row['graph_name']): <40}"
              f"{get_category(row['graph_name']): <40}")
    # return
    if download:
        download_and_extract.main(rows, io_modes)

    if orders:
        if orders == ['all']:
            orders = config.settings['orderings'].keys()
        # verify that the requested ordering to compute are supported
        assert valid_orderings(orders)
        reorder.main(rows, orders)

    if plot:
        plotting.main(rows, orders)


    if run_pr_expts:
        pr_experiments.main(rows, list(orders) + ['orig'])
        # DEBUG - plot the edge orderings
        if debug:
            for vorder_str in orders:
                plotting.plot_edge_orderings(rows, vorder_str)
            assert valid_pr(rows)

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

    parser.add_argument('-g', '--graph-names', nargs='+',
                        help='If specified, only download and scrape these graphs '
                             'from dataset (must exist in datasets.json). '
                             'Otherwise, download ALL graphs in datasets.json.', )

    parser.add_argument('-m', '--io-modes', nargs='+',
                        help='The IO mode that the graphs and isomorphisms will be saved as.'
                             'At least one of [binary, text] should be specified', )

    parser.add_argument('-d', '--download',
                        action=argparse.BooleanOptionalAction, required=True,
                        help='Whether to download the graph from konect or not.')

    parser.add_argument('-r', '--reorder', nargs='+',
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
                             '  `all:  All the above orderings.\n'
                             '}')

    parser.add_argument('-l', '--plot', action=argparse.BooleanOptionalAction,
                        help='Whether to plot the adjacency matrices or not.', required=True)

    parser.add_argument('-e', '--run-pr-expts',
                        action=argparse.BooleanOptionalAction, required=True,
                        help='Whether to run Edge-Centric PageRank computation experiments or not.')

    parser.add_argument('--debug',
                        action=argparse.BooleanOptionalAction, required=True,
                        help='Run in Debug Mode.'
                        )

    main(parser.parse_args())
