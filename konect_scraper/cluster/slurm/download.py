

import argparse
import logging
from pathlib import Path
import pandas as pd

from ... import config
import os
from ... import download_and_extract
from ...util import single_val_numeric_set 

def main(config_path, config_idx):
    settings = config.settings
    print(config_path, config_idx)
    df = pd.read_csv(config_path)
    row = df.iloc[int(config_idx)]
    graph_name = row['graph_name']
    konect_url = row['konect_url']
    data_url = row['data_url']
    graph_dir = os.path.join(settings['graphs_dir'], graph_name)
    compressed_edge_list_path = os.path.join(graph_dir, f"{settings['compressed_el_file_name']}.net")
    if Path(compressed_edge_list_path).is_file():
        logging.info(f"{graph_name} already compressed; skipping.")
        return 

    # create directory for the graph (if not exists)
    Path(graph_dir).mkdir(parents=True, exist_ok=True)
    download_and_extract.download_graph(data_url, graph_dir)
    # update the graph's metadata to point to its directory
    download_and_extract.update_meta_dir(graph_dir, graph_name)

    # update the text file size in the database
    graph_path = os.path.join(graph_dir, settings["orig_el_file_name"])
    single_val_numeric_set('txt_file_size', 'metadata', graph_name, os.path.getsize(graph_path))

    return 

if __name__ == '__main__':
    config.init()
    parser = argparse.ArgumentParser()
    parser.add_argument('config_path')
    parser.add_argument('config_idx')
    args = parser.parse_args()
    main(args.config_path, args.config_idx)
    