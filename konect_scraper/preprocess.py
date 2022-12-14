import logging
import os
from pathlib import Path
from konect_scraper import config
from konect_scraper.download_and_extract import compress
from konect_scraper.util import get_directed, get_size, get_size_in_memory, get_volume, save_connected_component, \
    save_ground_truth_pr, set_n_m, single_val_numeric_set, save_peregrine, save_webgraph
import networkx as nx
import igraph as ig
import numpy as np
import numpy as np

from scipy.linalg import eig, eigh

from scipy.sparse.linalg import eigs, eigsh
def process(graph_name, io_modes, overwrite):

    settings = config.settings
    graphs_dir = settings['graphs_dir']
    graph_dir = os.path.join(graphs_dir, graph_name)

    compressed_edge_list_path = os.path.join(
        graph_dir, f"{settings['compressed_el_file_name']}.net")

    compressed_graph_path = os.path.join(
        graph_dir, settings["compressed_el_file_name"])

    compressed_graph_path_with_extension = \
        os.path.join(
            f"{compressed_graph_path}.{settings['edgelist_file_suffix']}")

    if Path(compressed_edge_list_path).is_file():
        if overwrite:
            logging.info(f"{graph_name} already compressed; Overwriting.")
        else:
            logging.info(f"{graph_name} already compressed; skipping.")
            return

    # compress the graph's edgelist
    directed = bool(get_directed(graph_name))
    sz = get_size(graph_name); vol = get_volume(graph_name)
    n, m = compress(graph_dir, directed, sz, vol, io_modes, graph_name)
    set_n_m(graph_name, n, m)

    # update the PageRank experiments structs size in db
    get_size_in_memory(n, m)

    # update the text file size in the database
    single_val_numeric_set('compressed_txt_file_size', 'metadata', graph_name,
                            os.path.getsize(compressed_graph_path + ".net"))

    # save the compressed edgelist's PageRank for verification of
    logging.info(f"Computing {graph_name}'s PageRank")
    save_ground_truth_pr(compressed_graph_path_with_extension, graph_name)

    logging.info(f"Preprocessing {graph_name} for peregrine..")
    save_peregrine(compressed_graph_path_with_extension, graph_dir)

    logging.info(f"Preprocessing {graph_name} for webgraph..")
    save_webgraph(compressed_graph_path_with_extension, graph_dir)


def main(rows, io_modes, overwrite):
    """
    1. Compress a graph's edge list
        - removes 
            - vertex, edge duplicates
            - self directed edges
        - compresses the vertex ID space to [0, N) (where N = num_vertices)
    2. Preprocess the graph's edge list using DBG utils
    3. Computes and persists the graph's PageRank (for future verification)
    4. Saves the graph's edgelist as a scipy.sparse.csr_matrix
        - to be used for:
            - computing linalg features of graph (e.g. eigenvalues)

    Args:
        rows (list[sqlite3.Row]): rows of konect table listing graphs
    """

    settings = config.settings
    graphs_dir = settings['graphs_dir']

    for row in rows:
        graph_name = row['graph_name']
        process(graph_name, io_modes, overwrite)

        

    return
