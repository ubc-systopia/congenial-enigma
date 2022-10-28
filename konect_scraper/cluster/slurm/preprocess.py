import argparse
import logging
from pathlib import Path
import pandas as pd
import os

from konect_scraper.cluster.main import parse_and_init_data_dir
from konect_scraper.cluster.slurm.utils import is_graph_compressed
from ... import config
import os
from ... import download_and_extract
from ...util import get_directed, get_size, get_size_in_memory, get_volume, set_n_m, single_val_numeric_set


def main(config_path, config_idx):
    settings = config.settings
    if is_graph_compressed(settings, config_path, config_idx):
        return

    df = pd.read_csv(config_path)
    row = df.iloc[int(config_idx)]
    graph_name = row['graph_name']
    graph_dir = os.path.join(settings['graphs_dir'], graph_name)

    # compress the graph's edgelist
    directed = bool(get_directed(graph_name))
    sz = get_size(graph_name)
    vol = get_volume(graph_name)
    io_modes = [config.IOMode.binary, config.IOMode.text]
    n, m = download_and_extract.compress(
        graph_dir, directed, sz, vol, io_modes, graph_name)
    set_n_m(graph_name, n, m)
    # update the PageRank experiments structs size in db
    get_size_in_memory(n, m)

    # update the text file size in the database
    compressed_graph_path = os.path.join(
        graph_dir, settings["compressed_el_file_name"])
    single_val_numeric_set('compressed_txt_file_size', 'metadata', graph_name,
                           os.path.getsize(compressed_graph_path + ".net"))

    return


if __name__ == '__main__':
    args = parse_and_init_data_dir()
    config.init(args.data_dir)
    main(args.config_path, args.config_idx)
