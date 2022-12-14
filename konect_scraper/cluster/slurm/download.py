import argparse
import logging
from pathlib import Path
import pandas as pd
from konect_scraper.cluster.main import parse_and_init_data_dir
from konect_scraper.cluster.slurm.utils import is_graph_compressed


from ... import config
import os
from ... import download_and_extract
from ...util import single_val_numeric_set


def main(config_path, config_idx):
    settings = config.settings
    if is_graph_compressed(settings, config_path, config_idx):
        return

    df = pd.read_csv(config_path)
    row = df.iloc[int(config_idx)]
    graph_name = row['graph_name']
    data_url = row['data_url']
    graph_dir = os.path.join(settings['graphs_dir'], graph_name)

    # create directory for the graph (if not exists)
    Path(graph_dir).mkdir(parents=True, exist_ok=True)
    download_and_extract.download_graph(data_url, graph_dir)
    # update the graph's metadata to point to its directory
    download_and_extract.update_meta_dir(graph_dir, graph_name)

    # update the text file size in the database
    graph_path = os.path.join(graph_dir, settings["orig_el_file_name"])
    single_val_numeric_set('txt_file_size', 'metadata',
                           graph_name, os.path.getsize(graph_path))

    return


if __name__ == '__main__':
    args = parse_and_init_data_dir()
    config.init(args.data_dir)
    main(args.config_path, args.config_idx)
