import argparse
import logging
from pathlib import Path
import pandas as pd

import os
from konect_scraper import config
from konect_scraper.cluster.main import parse_and_init_data_dir
from konect_scraper.reorder import compute_ordering

def main(config_path, config_idx):
    df = pd.read_csv(config_path)
    row = df.iloc[int(config_idx)]
    graph_name = row['graph_name']
    vertex_order = row['vertex_order']
    overwrite = bool(row['overwrite'])
    compute_ordering(graph_name, vertex_order, overwrite)
    return 


if __name__ == '__main__':
    args = parse_and_init_data_dir()
    config.init(args.data_dir)
    main(args.config_path, args.config_idx)
