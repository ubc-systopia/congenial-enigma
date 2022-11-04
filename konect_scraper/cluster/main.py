import argparse
import logging
from datetime import datetime

from ..config import *
from ..column_names import *
from ..util import *
from . import execute
from konect_scraper.cluster.util import verify_vertex_orders


def get_io_modes(io_modes):
    # if io mode is unspecified, use both binary and text as default
    if not io_modes:
        io_modes = [config.IOMode.binary, config.IOMode.text]
    else:
        modes = []
        for io_mode in io_modes:
            match io_mode:
                case 'binary':
                    modes.append(config.IOMode.binary)
                case 'text':
                    modes.append(config.IOMode.text)
                case _:
                    logging.error(f"{io_mode}: Unsupported IO mode!")
        io_modes = modes
    return io_modes


def parse_and_init_data_dir():
    prsr = argparse.ArgumentParser()
    prsr.add_argument('config_path')
    prsr.add_argument('config_idx')
    prsr.add_argument('data_dir')
    return prsr.parse_args()


def main(args):

    create_log_dir_if_not_exists()
    plot = args.plot
    overwrite = args.overwrite
    orders = args.reorder
    io_modes = args.io_modes
    directed = args.directed
    json_args_path = args.json_args
    debug = args.debug
    graph_ns = list(map(int, args.graph_numbers))
    exec_mode = args.mode

    io_modes = get_io_modes(io_modes)

    config.settings['debug'] = debug

    log_dir = config.settings['logging']['log_dir']
    curr_time = datetime.now().strftime("%H_%d_%m_%Y")
    log_path = f"{config.settings['app_name']}_{curr_time}"
    log_file_name = os.path.join(log_dir, log_path + '.' + 'log')
    init_logger(log_file_name)

    slurm_params = {
        'time': args.time,
        'mem': args.mem,
        'cpus-per-task': args.cpus_per_task,
    }

    if args.constraint:
        slurm_params['constraint'] = args.constraint

    if directed:
        graph_type = 'directed'
    else:
        graph_type = 'undirected'

    match exec_mode:
        case 'download':
            print(
                f"Downloading {graph_ns[1] - graph_ns[0]} graphs to {config.settings['graphs_dir']}")
            execute.main(graph_type, graph_ns, slurm_params, 'download')
            return
        case 'preprocess':
            execute.main(graph_type, graph_ns, slurm_params, 'preprocess')
            return
        case 'reorder':
            orders = verify_vertex_orders(orders, config.settings)
            execute.main(graph_type, graph_ns, slurm_params,
                         'reorder', orders, overwrite)
            return
        case 'plot':
            # todo
            print("Unimplemented")
            return
        case 'pr-expt':
            orders = verify_vertex_orders(orders, config.settings)
            execute.main(graph_type, graph_ns, slurm_params, 'pr_expt',
                         list(orders) + ['orig'])
            return
        case _:
            print(f"Unsupported execution mode on cluster: {exec_mode}")
            return


if __name__ == '__main__':

    argparse_desc = """
    Konect Scraper: a driver that performs the following functions on:
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
    Meant to be executed on a cluster using Slurm - 1-5 are executed using 
    `sbatch`
    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=argparse.RawTextHelpFormatter)

    exec_modes = {'download', 'preprocess', 'reorder', 'plot', 'pr-expt'}
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

    parser.add_argument('-a', '--data-dir', help='Absolute path to the '
                        'directory that should store all graphs, orderings, plots, etc.\n'
                        'If supplied, will overwrite path given by config.py')

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

    parser.add_argument('-o', '--overwrite', action='store_true',
                        help='If true, previous graphs/compressed graphs/vertex'
                        'orders will be redownloaded/recompressed/re-reordered.')

    parser.add_argument('-l', '--plot', action=argparse.BooleanOptionalAction,
                        help='Whether to plot the adjacency matrices or not.',)

    parser.add_argument('--debug',
                        action=argparse.BooleanOptionalAction,
                        help='Run in Debug Mode.'
                        )
    parser.add_argument('--json-args',
                        help="Path to json file containing arguments listing which graphs"
                             "to preprocess, run experiments on etc.")

    # required slurm params
    parser.add_argument('--time', required=True,
                        help='maximum time given to each sbatch job. DD-HH:MM:SS')

    parser.add_argument('--mem', required=True,
                        help='minimum memory required for each sbatch job')

    parser.add_argument('--cpus-per-task', required=True,
                        help='Since task are exclusively run on single nodes'
                        ', e.g. --nodes=1-1, cpus-per-task dictates the number'
                        'of cores used.')

    parser.add_argument('--constraint',
                        help='One of broadwell, cascade, skylake. If any are sufficient, can be specified as e.g. '
                        '[skylake|cascade]')

    args = parser.parse_args()

    # overwrite user supplied data directory in config
    data_dir = args.data_dir
    config.init(data_dir)
    column_names.init()

    main(args)
