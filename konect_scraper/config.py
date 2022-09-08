import os.path
from pathlib import Path
from enum import Enum


class IOMode(Enum):
    binary = 1
    text = 2


def init():
    global settings

    # TODO update repo home
    repo_root = Path(os.path.dirname(os.path.realpath(__file__))).parent
    app_name = "graph_preprocess"
    all_networks_url = "http://konect.cc/networks/"
    repo_home = os.path.join(repo_root, "konect_scraper")
    rabbit_home = os.path.join(repo_root, "rabbit_order")
    dbg_home = os.path.join(repo_root, "dbg")

    data_dir = os.path.join(repo_home, "data")
    sqlite3_db_path = os.path.join(data_dir, "graphs.db")
    datasets_json_path = os.path.join(repo_home, "datasets.json")
    dataframes_dir = os.path.join(data_dir, "dataframes")
    graphs_dir = os.path.join(data_dir, "graphs")
    plots_dir = os.path.join(data_dir, "plots")
    orig_el_file_name = "orig.net"
    compressed_el_file_name = "comp"
    cmake_build_dir = "cmake-build-debug"
    graph_preprocess_dir = os.path.join(repo_root, "graph_preprocess")
    # EXECUTABLES
    graph_preprocess_executable = os.path.join(graph_preprocess_dir, cmake_build_dir, "graph_preprocess")
    slashburn_executable = os.path.join(graph_preprocess_dir, cmake_build_dir, "slashburn")
    cuthill_mckee_executable = os.path.join(graph_preprocess_dir, cmake_build_dir, "cuthill_mckee")
    rabbit_order_executable = os.path.join(rabbit_home, "demo", "reorder")

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
    dpi = 200
    ax_size = 5  # sidelength of an ax in a matrix of plots; used to calculate the total figure size

    settings = {
        "repo_root": repo_root,
        "repo_home": repo_home,
        "app_name": app_name,
        "all_networks_url": all_networks_url,
        "sqlite3": {
            # "tables": {
            #     ["metadata", "statistics"]
            # },
            "sqlite3_db_path": sqlite3_db_path
        },
        "datasets_json_path": datasets_json_path,
        "data_dir": data_dir,
        "dataframes_dir": dataframes_dir,
        "graphs_dir": graphs_dir,
        "plots_dir": plots_dir,
        "orig_el_file_name": orig_el_file_name,
        "compressed_el_file_name": compressed_el_file_name,

        # Executables
        "graph_preprocess_executable": graph_preprocess_executable,
        "slashburn_executable": slashburn_executable,
        "cuthill_mckee_executable": cuthill_mckee_executable,
        "rabbit_order_executable": rabbit_order_executable,

        "comment_strings": ["%", "#"],
        "plot": {
            "marker": marker,
            "markersize": markersize,
            "format": plot_format,
            "dpi": dpi,
            "adj_mat_format": "spy",  # or matshow,
            "ax_size": ax_size,
            "bbox_inches": 'tight',
            "pad_inches": 0,
            "max_rows_per_agg_spy_plot": 4,  # the number of graphs to show per aggregated spy plot,
            "max_n": 20_000,  # the largest graph size that is plottable as as adjacency matrix
        },
        "orderings": {
            'rnd': "random",
            'rbt': "rabbit",
            'sb': "slashburn",
            'cm': "cuthill-mckee",
            'srt': "sort",
            'hc': "hubcluster",
            'hs': "hubsort",
            'dbg': "degree-based-grouping"
        },
        "logging": {
            "log_dir": log_dir,
            "log_format": log_format
        },
        "hyperparameters": {
            "slashburn": {
                "percent": slashburn_percent
            },
            "cuthill-mckeee": {

            }
        }
    }
