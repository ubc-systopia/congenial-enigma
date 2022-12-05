import json
import logging
import sqlite3
import subprocess
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from konect_scraper import config, column_names
from konect_scraper.io import find
from konect_scraper.sql import connect, append_df_to_table
from konect_scraper.stats import compute_deg_stats, \
    compute_plfit_stats, compute_scipy_stats, compute_radii, percentile_effective_diameter, compute_distance_stats, \
    hyperball, compute_motif_stats
from konect_scraper.scrape_konect_stats import write_to_sqlite3
import numpy as np
import sys
from scipy.special import comb

from konect_scraper.util import get_n


def compute_stats(graph_name):
    n = get_n(graph_name)
    stats = {}
    logging.info(f"Computing {graph_name}'s algebraic stats..")
    stats.update(compute_scipy_stats(graph_name))
    logging.info(f"Computing {graph_name}'s degree stats..")
    stats.update(compute_deg_stats(graph_name))
    logging.info(f"Computing {graph_name}'s Powerlaw stats..")
    stats.update(compute_plfit_stats(graph_name))
    logging.info(f"Computing {graph_name}'s distance stats..")
    stats.update(compute_distance_stats(graph_name))
    logging.info(f"Computing {graph_name}'s motif stats..")
    stats.update(compute_motif_stats(graph_name))
    # TODO in gini coefficient

    # for k, v in stats.items():
        # print(k, v)

    diff = set(stats.keys()).difference(set(column_names.features_col_names))
    return stats

def main(rows):
    settings = config.settings
    graphs_dir = settings['graphs_dir']
    compressed_fname = settings['compressed_el_file_name']
    edgelist_file_suffix = settings['edgelist_file_suffix']
    settings['scipy_csr_suffix']
    # compute the given orders for each of the datasets
    conn = connect()
    # df = pd.concat([df, pd.DataFrame([d])], ignore_index=True)

    for _, row in enumerate(rows):
        feats_df = pd.DataFrame(columns=column_names.features_col_names.keys())

        graph_name = row['graph_name']
        d = compute_stats(graph_name)
        d['graph_name'] = graph_name
        feats_df = pd.concat([feats_df, pd.DataFrame([d])], ignore_index=True)

    
    
        db_path = config.settings['sqlite3']['sqlite3_db_path']
        timeout = config.settings['sqlite3']['timeout']
        conn = sqlite3.connect(db_path, timeout=timeout)
        write_to_sqlite3(feats_df, 'features', conn)
    # append_df_to_table(feats_df, 'features')
    # print(f"{feats_df=}")
    # feats_df.to_csv("tmp_feats.csv")
    return
