import os.path


def init():
    global settings
    # TODO update repo home
    repo_root = "/home/atrostan/Workspace/repos/congenial-enigma/"
    app_name = "graph_preprocess"
    repo_home = os.path.join(repo_root, "konect_scraper")
    data_dir = os.path.join(repo_home, "data")
    sqlite3_db_path = os.path.join(data_dir, "graphs.db")
    datasets_json_path = os.path.join(repo_home, "datasets.json")
    dataframes_dir = os.path.join(data_dir, "dataframes")
    graphs_dir = os.path.join(data_dir, "graphs")
    plots_dir = os.path.join(data_dir, "plots")
    orig_el_file_name = "orig.net"
    compressed_el_file_name = "comp.net"
    cmake_build_dir = "cmake-build-debug"
    graph_preprocess_executable = os.path.join(repo_root, "graph_preprocess", cmake_build_dir, "graph_preprocess")

    log_dir = os.path.join(repo_root, "logs")

    settings = {
        "repo_home": repo_home,
        "app_name": app_name,
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
        "graph_preprocess_executable": graph_preprocess_executable,
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
            "log_dir": log_dir
        }
    }
