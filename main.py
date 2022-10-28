import argparse
import json
import logging
import os
from argparse import RawTextHelpFormatter
from datetime import datetime

import konect_scraper.column_names as column_names
import konect_scraper.config as config
import konect_scraper.plot as plotting
from konect_scraper import download_and_extract, reorder, pr_experiments
from konect_scraper.config import IOMode
from konect_scraper.sql import get_all_graphs_by_graph_names_where_stats_between, get_all_downloadable_graphs, \
    get_all_unipartite_directed_graphs, get_all_graphs_by_graph_names, get_all_unipartite_graphs, \
    get_all_unipartite_undirected_graphs
from konect_scraper.util import \
    create_log_dir_if_not_exists, init_logger, valid_orderings, valid_pr, get_category, get_pr_struct_size, \
    get_unimputed_features, get_directed, get_n, get_m, get_n_vertices, get_n_edges, convert_size

from konect_scraper.cluster import execute as cc_download

def get_io_modes(io_modes):
    # if io mode is unspecified, use both binary and text as default
    if not io_modes:
        io_modes = [IOMode.binary, IOMode.text]
    else:
        modes = []
        for io_mode in io_modes:
            match io_mode:
                case 'binary':
                    modes.append(IOMode.binary)
                case 'text':
                    modes.append(IOMode.text)
                case _:
                    logging.error(f"{mode}: Unsupported IO mode!")
        io_modes = modes
    return io_modes

def main(args):
    create_log_dir_if_not_exists()
    plot = args.plot
    orders = args.reorder
    io_modes = args.io_modes
    directed = args.directed
    run_pr_expts = args.run_pr_expts
    json_args_path = args.json_args
    debug = args.debug
    graph_ns = list(map(int, args.graph_numbers))
    exec_mode = args.mode
    environment = args.environment

    io_modes = get_io_modes(io_modes)

    config.settings['debug'] = debug

    log_dir = config.settings['logging']['log_dir']
    curr_time = datetime.now().strftime("%H_%d_%m_%Y")
    log_path = f"{config.settings['app_name']}_{curr_time}"
    log_file_name = os.path.join(log_dir, log_path + '.' + 'log')
    init_logger(log_file_name)

    if directed:
        graph_type = 'directed'
    else:
        graph_type = 'undirected'

    
    directed = True
    if directed:
        rows = get_all_unipartite_directed_graphs()
    else:
        rows = get_all_unipartite_undirected_graphs()
    n_unipartite_graphs = len(rows)
    graph_names = [r['graph_name'] for r in rows]
    print(f"{len(graph_names)} unipartite graphs in dataset.")

    # a selection of relevant graph_names have been identified,
    # download, reorder, and plot them

    graph_names = [r['graph_name'] for r in rows]
    rows = get_all_downloadable_graphs(graph_names)[:]
    graph_names = get_unimputed_features([r['graph_name'] for r in rows])
    if json_args_path:
        with open(json_args_path, 'r') as j:
            json_args = json.loads(j.read())
        graph_names = json_args['graph_names']
    rows = get_all_graphs_by_graph_names(graph_names)
    # print(rows)
    graph_name_start_idx = 283
    graph_name_end_idx = -2
    # graph_name_end_idx = len(rows)
    rows = sorted(rows, key=lambda r: get_pr_struct_size(
        r['graph_name']), reverse=False)

    rows = rows[graph_name_start_idx:graph_name_end_idx]

    # print heading
    print(f"{'Index' : <5} {'Graph Name' : <40}"
          f"{'|V|': <40}"
          f"{'|E|': <40}"
          f"{'PageRank Struct Size ': <40}"
          f"{'Graph Category': <40}"
          f"{'Is Directed?': <40}")
    print('-' * 245)
    for i, row in enumerate(rows):
        # if int(get_directed(row['graph_name'])) == 1:
        #     continue
        print(f"{i + graph_name_start_idx: <5} {row['graph_name'] : <40}"
              f"{get_n_vertices(row['graph_name']): <40}"
              f"{get_n_edges(row['graph_name']): <40}"
              f"{str(size) + size_string: <40}"
              f"{get_category(row['graph_name']): <40}"
              f"{get_directed(row['graph_name']): <40}")
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
    Konect Scraper: a driver that performs the following functions:
    1. Downloads graphs from konect.cc
    2. Preprocesses downloaded graphs
        - removes duplicate edges, self-loops
        - maps the vertex ids to a dense space: [0, |V| - 1)
    3. Reorders preprocessed graphs 
        - each graph directory will store the isomorphism map for the computed
          ordering
    4. Plots computed isomorphism (i.e. vertex orderings)
    5. Runs PageRank Experiments
        - each experiment consists of iterating over all the edges of the graph
          using a pair of <vertex order, edge order> combination
    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=RawTextHelpFormatter)

    parser.add_argument('-v', '--environment', choices={'local', 'cluster'},
                        required=True, help="Whether to run konect scraper and"
                        "preprocess locally or on Compute Canada")

    exec_modes = {'download', 'preprocess', 'reorder', 'plot', 'pr_expt'}
    parser.add_argument('-m', '--mode', choices=exec_modes, required=True,
                        help="Specify the execution mode. One of: "
                        "\{'download', 'preprocess', 'reorder', 'plot', "
                        "'pr_expt'\}")

    parser.add_argument('-d', '--directed',
                        action=argparse.BooleanOptionalAction,
                        required=True,
                        help='Whether to download and preprocess directed/undirected graphs.\n'
                        '(Bipartite graphs currently unsupported).')

    parser.add_argument('-g', '--graph-numbers', nargs='+', required=True,
                        help='If specified, only download and scrape these graphs\n'
                        'e.g. `--directed -g 0 100` would download all directed graphs whose graph number'
                        'is [0, 100) ')

    parser.add_argument('-i', '--io-modes', nargs='+',
                        help='The IO mode that the graphs and isomorphisms will be saved as.'
                             'At least one of [binary, text] should be specified', )

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
                        help='Whether to plot the adjacency matrices or not.',)

    parser.add_argument('-e', '--run-pr-expts',
                        action=argparse.BooleanOptionalAction,
                        help='Whether to run Edge-Centric PageRank computation experiments or not.')

    parser.add_argument('--debug',
                        action=argparse.BooleanOptionalAction,
                        help='Run in Debug Mode.'
                        )
    parser.add_argument('--json-args',
                        help="Path to json file containing arguments listing which graphs"
                             "to preprocess, run experiments on etc.")

    main(parser.parse_args())


"""Examples

Only download and preprocess (i.e. simplify graph):
$ python main.py --download --io-modes binary text --no-plot --no-run-pr-expts --no-debug

No download, just reorder

$ python main.py --no-download --io-modes binary text --reorder all --no-plot --no-run-pr-expts --no-debug

No download, just experiments (reorder all detects that vertex orderings have been computed, and skips recomputation)

$ python main.py --no-download --io-modes binary text --reorder all --no-plot --run-pr-expts --no-debug


Download, preprocess, reorder, and plot; No PR experiments;
$ python main.py --download --io-modes binary text --reorder all --plot --no-run-pr-expts --no-debug

$ python main.py --download --io-modes binary text --reorder all --no-plot --run-pr-expts --no-debug


    1.  Download either:
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
