import os.path
import subprocess
import sys
from pathlib import Path
from enum import Enum

import psutil
from PyQt5.QtWidgets import QApplication

from konect_scraper.util import get_cache_stats


def get_monitors_dpi():
    print(os.environ['DISPLAY'])
    app = QApplication(sys.argv)
    screens = app.screens()
    dpis = [screen.physicalDotsPerInch() for screen in screens]
    app.quit()
    return dpis


def get_n_threads():
    return psutil.cpu_count()


class IOMode(Enum):
    binary = 1
    text = 2


def init(input_data_dir=None):
    global settings

    debug = False

    # TODO update repo home
    repo_root = Path(os.path.dirname(os.path.realpath(__file__))).parent

    repo_name = "congenial-enigma"
    app_name = "graph_preprocess"
    all_networks_url = "http://konect.cc/networks/"
    repo_home = os.path.join(repo_root, "konect_scraper")
    rabbit_home = os.path.join(repo_root, "rabbit_order")

    pbrcm_home = os.path.join(repo_root, "ParallelBatchRCM")

    dbg_home = os.path.join(repo_root, "dbg")  # todo replace

    dbg_apps_dir = os.path.join(dbg_home, "apps")
    dbg_convert_dir = os.path.join(dbg_home, "graph-convert-utils")
    dbg_clean_el_executable = os.path.join(
        dbg_convert_dir, "clean_edgelist.py")
    dbg_convert_script = os.path.join(dbg_convert_dir, "convert.sh")
    dbg_convert_executable = os.path.join(dbg_convert_dir, "convert")
    dbg_datasets_dir = os.path.join(dbg_home, "datasets")

    if input_data_dir:
        data_dir = input_data_dir
    else:
        data_dir = os.path.join(repo_home, "data")
    # data_dir = '/media/atrostan/patterson_backup/data/'

    sqlite3_db_path = os.path.join(data_dir, "graphs.db")
    # sqlite3_db_path = "/home/atrostan/workspace/data/graphs.db"# todo replace
    datasets_json_path = os.path.join(repo_home, "datasets.json")
    dataframes_dir = os.path.join(repo_root, "konect_dataframes")
    graphs_dir = os.path.join(data_dir, "graphs")
    plots_dir = os.path.join(data_dir, "plots")
    results_dir = os.path.join(data_dir, "results")

    edgelist_file_suffix = 'net'
    scipy_csr_suffix = 'npz'

    orig_el_file_name = f"orig.{edgelist_file_suffix}"
    compressed_el_file_name = "comp"
    pagerank_file_name = "pr"

    cmake_build_dir = "cmake-build-debug"
    graph_preprocess_dir = os.path.join(repo_root, "graph_preprocess")

    # EXECUTABLES
    graph_preprocess_executable = os.path.join(
        graph_preprocess_dir, cmake_build_dir, "graph_preprocess")
    slashburn_executable = os.path.join(
        graph_preprocess_dir, cmake_build_dir, "slashburn")
    cuthill_mckee_executable = os.path.join(
        graph_preprocess_dir, cmake_build_dir, "cuthill_mckee")
    convert_map_to_bin_executable = os.path.join(
        graph_preprocess_dir, cmake_build_dir, "convert_map_to_binary")
    compute_ccs_executable = os.path.join(
        graph_preprocess_dir, cmake_build_dir, "compute_ccs")
    stats_executable = os.path.join(
        graph_preprocess_dir, cmake_build_dir, "stats")

    parallel_batch_rcm_executable = os.path.join(
        pbrcm_home, 'build', 'CuthillMcKee')

    rabbit_cmake_build_dir = os.path.join(rabbit_home, "demo", cmake_build_dir)
    rabbit_order_executable = os.path.join(rabbit_cmake_build_dir, "reorder")

    pr_experiments_executable = os.path.join(
        graph_preprocess_dir, cmake_build_dir, "pr_experiments")

    # abseil and parallel slashburn
    par_slashburn_dir = os.path.join(repo_root, 'par_slashburn')
    abseil_repo_dir = os.path.join(par_slashburn_dir, 'abseil-cpp')
    abseil_install_include_dir = os.path.join(
        par_slashburn_dir, 'install', 'include')
    par_slashburn_executable = os.path.join(
        par_slashburn_dir, cmake_build_dir, 'par_slashburn')
    pr_executable = os.path.join(par_slashburn_dir, cmake_build_dir, 'pr')

    webgraph_dir = os.path.join(repo_root, 'webgraph')

    prgrn_dir = os.path.join(repo_root, 'peregrine')
    prgrn_bin_dir = os.path.join(prgrn_dir, 'bin')
    prgrn_convert_data_executable = os.path.join(prgrn_bin_dir, 'convert_data')
    prgrn_count_executable = os.path.join(prgrn_bin_dir, 'count')

    # LOGGING
    log_dir = os.path.join(repo_root, "logs")
    # log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    log_format = "[%(asctime)s %(filename)s->%(funcName)s():%(lineno)s]%(levelname)s: %(message)s"

    # HYPERPARAMETERS
    slashburn_percent = 0.005

    # PLOTTING
    marker = ','
    markersize = .5
    plot_format = 'png'
    try:
        dpis = get_monitors_dpi()
        dpi = max(list(map(int, get_monitors_dpi())))
    except:
        dpi = 200

    ax_size = 5  # sidelength of an ax in a matrix of plots; used to calculate the total figure size

    # CPU INFO
    cache_stats = get_cache_stats()  # todo parse larger table if lscpu -C doesn't work
    if not cache_stats:
        cache_stats['line_size'] = -1
        cache_stats['l1d_size'] = -1
        cache_stats['l2_size'] = -1
        cache_stats['l3_size'] = -1

    settings = {
        "debug": debug,
        "repo_root": repo_root,
        "repo_home": repo_home,
        "repo_name": repo_name,
        "app_name": app_name,

        "rabbit_home": rabbit_home,

        "dbg_home": dbg_home,
        "dbg_apps_dir": dbg_apps_dir,
        "dbg_convert_dir": dbg_convert_dir,
        "dbg_clean_el_executable": dbg_clean_el_executable,
        "dbg_convert_script": dbg_convert_script,
        "dbg_convert_executable": dbg_convert_executable,
        "dbg_datasets_dir": dbg_datasets_dir,

        "pbrcm_home": pbrcm_home,

        "make_executable": "make",
        "cmake_executable": "cmake",

        # cmake
        "graph_preprocess_dir": graph_preprocess_dir,
        "cmake_build_dir": cmake_build_dir,
        "rabbit_cmake_build_dir": rabbit_cmake_build_dir,
        "cmake_build_type": "Debug",
        "cmake_make_program": "ninja",

        "all_networks_url": all_networks_url,
        "sqlite3": {
            # "tables": {
            #     ["metadata", "statistics"]
            # },
            "sqlite3_db_path": sqlite3_db_path,
            "timeout": 60,
        },
        "datasets_json_path": datasets_json_path,
        "data_dir": data_dir,
        "dataframes_dir": dataframes_dir,
        "graphs_dir": graphs_dir,
        "plots_dir": plots_dir,
        "results_dir": results_dir,
        "orig_el_file_name": orig_el_file_name,
        "compressed_el_file_name": compressed_el_file_name,
        "pagerank_file_name": pagerank_file_name,

        # file suffixes
        "edgelist_file_suffix": edgelist_file_suffix,
        "scipy_csr_suffix": scipy_csr_suffix,
        "binary_suffix": "bin",

        # abseil
        "abseil_install_include_dir": abseil_install_include_dir,
        "abseil_repo_url": "https://github.com/abseil/abseil-cpp.git",
        "abseil_repo_dir": abseil_repo_dir,
        # parallel slashburn
        "par_slashburn_dir": par_slashburn_dir,

        "webgraph_dir": webgraph_dir,

        # Executables
        "graph_preprocess_executable": graph_preprocess_executable,
        "slashburn_executable": slashburn_executable,
        "par_slashburn_executable": par_slashburn_executable,
        "cuthill_mckee_executable": cuthill_mckee_executable,
        "parallel_batch_rcm_executable": parallel_batch_rcm_executable,
        "rabbit_order_executable": rabbit_order_executable,
        "pr_experiments_executable": pr_experiments_executable,
        "pr_executable": pr_executable,
        "convert_map_to_bin_executable": convert_map_to_bin_executable,
        'compute_ccs_executable': compute_ccs_executable,
        'stats_executable': stats_executable,

        "comment_strings": ["%", "#"],
        "n_threads": get_n_threads(),
        "plot": {
            "marker": marker,
            "markersize": markersize,
            "format": plot_format,
            "dpi": dpi,
            "adj_mat_format": "spy",  # or matshow,
            "ax_size": ax_size,
            "bbox_inches": 'tight',
            "pad_inches": 0,
            # the number of graphs to show per aggregated spy plot,
            "max_rows_per_agg_spy_plot": 5,
            "max_n": 100_000,  # the largest graph size that is plottable as as adjacency matrix
        },
        "orderings": {
            'rnd': "random",
            'rbt': "rabbit",
            'sb': "slashburn",
            'parsb': "par_slashburn",
            'cm': "cuthill-mckee",
            'rev_cm': "reverse-cuthill-mckee",
            'srt': "sort",
            'hc': "hubcluster",
            'hs': "hubsort",
            'dbg': "degree-based-grouping",
            'orig': "original",
        },
        "edge_orderings": {
            'row': "Row",
            'column': "Column",
            'hilbert': "Hilbert",
            # 'fgf': 'FastGeneralForm',
        },
        "logging": {
            "log_dir": log_dir,
            "slurm_log_dir": os.path.join(log_dir, 'slurm'),
            "log_format": log_format
        },
        "hyperparameters": {
            "slashburn": {
                "percent": slashburn_percent
            },
            "cuthill-mckeee": {

            },
            "pr-experiments": {
                "damping_factor": 0.85,
                "edge_orderings": [
                    'row',
                    'column',
                    'hilbert'
                ],
                "num_iters": 20,
                "num_expts": 5,
            }
        },
        "dbg": {
            "order_idx_dict": {
                'Random': 1,
                'Sort': 2,
                'HubSort': 3,
                'HubCluster': 4,
                'DBG': 5,
            },
            "order_str_dict": {
                'rnd': 'Random',
                'srt': 'Sort',
                'hs': 'HubSort',
                'hc': 'HubCluster',
                'dbg': 'DBG',
            },
            "degree_used_for_reordering": 0,
            # only run 1 iteration of PR - we're interested in the ordering,
            # not the PR
            'max_iters': 1,
            'edgelist_file_suffix': 'el',
        },
        "cpu-info": {
            "cache-sizes": {
                'line_size': cache_stats['line_size'],
                'l1d_size': cache_stats['l1d_size'],
                'l2_size': cache_stats['l2_size'],
                'l3_size': cache_stats['l3_size'],
            }
        },
        "modelling": {
            "proportion": 0.5,
            # need at least this many data samples to build a dataset for the PR
            "min_n_data_samples": 50,
            # expts; this value will be used to decide which features we'll use
            # to train a predictive model of vertex+edge ordering performance
        },
        'compute_canada': {
            'job_array_dir': os.path.join(repo_home, 'cluster', 'csvs'),
            'scripts_dir': os.path.join(repo_home, 'cluster', 'scripts'),
            'repo_root': repo_root,
            'data_dir': data_dir,
            'account': 'def-mseltzer',
            'user': 'atrostan',
            'image': '/home/atrostan/singularity-images/congenial_enigma.sif',
            'execution_modes': [
                'download',
                'preprocess',
                'reorder',
                'plot',
                'pr_expt',
                'all'
            ]
        },
        'webgraph': {
            'hyperball': {
                # the logarithm of the number of registers per counter.
                'log2m': 10, 
                'gc': 'CMS', # G1 (newer) or CMS (as in docs)
            }
        },
        'peregrine': {
            'prgrn_dir': prgrn_dir,
            'prgrn_bin_dir': prgrn_bin_dir,
            'prgrn_convert_data_executable': prgrn_convert_data_executable,
            'prgrn_count_executable': prgrn_count_executable,
            'motifs': {
                # a map from 3-motifs peregrine pattern strings to their label
                '3': {
                    '[1-3](1~2)[2-3]': 'vertex_induced_wedge',
                    '[1-2][1-3][2-3]': 'triangle',
                    '[1-3][2-3]': 'edge_induced_wedge'
                },
                # a map from 4-motifs peregrine pattern strings to their label
                '4': {
                    '[1-2][1-4](1~3)[2-3](2~4)(3~4)': 'path_4',
                    '[1-2][1-3][1-4](2~3)(2~4)(3~4)': 'star_4',
                    '[1-2][1-3][1-4][2-3](2~4)(3~4)': 'tailed_triangle',
                    '[1-2][1-4](1~3)[2-3](2~4)[3-4]': 'cycle_4',
                    '[1-2][1-3][1-4][2-3](2~4)[3-4]': 'chordal_cycle_4',
                    '[1-2][1-3][1-4][2-3][2-4][3-4]': 'clique_4',
                },
                '5': {
                    '[1-2][1-3][1-4][1-5][2-3][2-4][2-5][3-4][3-5][4-5]': 'clique_5',
                }
            }
        },
        'plfit': {
            'stats': {
                'moments': {
                    'mean': 'mean',
                    'variance': 'variance',
                    'std.dev.': 'stddev',
                    'skewness': 'skewness',
                    'kurtosis': 'kurtosis',
                },
                'mle': {
                    'alpha': 'alpha',
                    'xmin': 'xmin',
                    'L': 'log_likelihood',
                    'D': 'ks_stat',
                    'p': 'p_value'
                }
            }
        }
    }
