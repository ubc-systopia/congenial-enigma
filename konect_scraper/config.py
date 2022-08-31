import os.path


def init():
    global settings
    repo_home = "./"
    repo_home = "/home/atrostan/PycharmProjects/tmp/konect_scraper"
    data_dir = os.path.join(repo_home, "data")
    sqlite3_db_path = os.path.join(data_dir, "graphs.db")
    datasets_json_path = os.path.join(repo_home, "datasets.json")
    dataframes_dir = os.path.join(data_dir, "dataframes")
    graphs_dir = os.path.join(data_dir, "graphs")
    orig_el_file_name = "orig.net"
    compressed_el_file_name = "comp.net"
    graph_simplify_executable = "/home/atrostan/CLionProjects/graph-simplify/cmake-build-debug/graph_simplify"
    settings = {
        "repo_home": repo_home,
        "sqlite3": {
            # "tables": {
            #     ["metadata", "statistics"]
            # },
            "sqlite3_db_path": sqlite3_db_path
        },
        "datasets_json_path": datasets_json_path,
        "dataframes_dir": dataframes_dir,
        "graphs_dir": graphs_dir,
        "orig_el_file_name": orig_el_file_name,
        "compressed_el_file_name": compressed_el_file_name,
        "graph_simplify_executable": graph_simplify_executable,
        "orderings": {
            'rnd': "random",
            'rbt': "rabbit",
            'sb': "slashburn",
            'cm': "cuthill-mckee",
            'srt': "sort",
            'hc': "hubcluster",
            'hs': "hubsort",
            'dbg': "degree-based-grouping"
        }
    }
