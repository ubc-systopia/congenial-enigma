import logging
from pathlib import Path
import pandas as pd
import os

def is_graph_compressed(settings, config_path, config_idx):
    df = pd.read_csv(config_path)
    row = df.iloc[int(config_idx)]
    graph_name = row['graph_name']
    konect_url = row['konect_url']
    data_url = row['data_url']
    graph_dir = os.path.join(settings['graphs_dir'], graph_name)
    compressed_edge_list_path = os.path.join(
        graph_dir, f"{settings['compressed_el_file_name']}.net")
    if Path(compressed_edge_list_path).is_file():
        logging.info(f"{graph_name} already compressed; skipping.")
        return True
    return False
