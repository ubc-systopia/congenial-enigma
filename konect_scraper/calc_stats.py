import json
import logging
import sqlite3
import subprocess
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from konect_scraper import config, column_names
from konect_scraper.io import find
from konect_scraper.sql import connect, append_df_to_table, insert_row_if_not_exists
from konect_scraper.stats import compute_deg_stats, \
    compute_plfit_stats, compute_scipy_stats, compute_radii, percentile_effective_diameter, compute_distance_stats, \
    hyperball, compute_motif_stats
from konect_scraper.scrape_konect_stats import write_to_sqlite3, update_to_sqlite3
import numpy as np
import sys
from scipy.special import comb
import os
from konect_scraper.util import get_n


def compute_stats(graph_name):
    n = get_n(graph_name)
    stats = {} 

    # create an empty row for the graph in the features sqlite3 table
    insert_row_if_not_exists(graph_name, 'features')

    db_path = config.settings['sqlite3']['sqlite3_db_path']
    graphs_dir = config.settings['graphs_dir']

    # cpp: compute connected components stats and write degree arrays if not 
    # exist
    logging.info(f"Computing {graph_name}'s cc stats and writing deg arrays..")

    graph_dir = os.path.join(graphs_dir, graph_name)
    cc_exec = config.settings['compute_ccs_executable']
    args = [
        cc_exec, 
        '-f', os.path.join(graph_dir, "comp.net"),
        '-o', os.path.join(graph_dir, "lcc.net"),
        '-s', 
        '-d', db_path
    ]
    logging.info(" ".join(args))
    res = subprocess.check_output(args)

    args = [
        cc_exec, 
        '-f', os.path.join(graph_dir, "comp.net"),
        '-o', os.path.join(graph_dir, "lscc.net"),
        '-d', db_path
    ]
    logging.info(" ".join(args))
    res = subprocess.check_output(args)
    # cpp: eigen stats

    logging.info(f"Computing {graph_name}'s eigen stats..")
    stats_executable = config.settings['stats_executable']
    args = [
        stats_executable, '-n', str(n), '-g', graph_dir, '-d', db_path
    ]
    logging.info(" ".join(args))
    res = subprocess.check_output(args)

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
        # return 
        feats_df = pd.concat([feats_df, pd.DataFrame([d])], ignore_index=True)

        d['graph_name'] = graph_name

        db_path = config.settings['sqlite3']['sqlite3_db_path']
        timeout = config.settings['sqlite3']['timeout']
        conn = sqlite3.connect(db_path, timeout=timeout)
        diff = set(column_names.features_col_names).difference(set(d.keys()))
        # some values have not been computed so update instead
        if (len(diff) > 0):
            update_to_sqlite3(d, 'features', conn)
        else:
            write_to_sqlite3(feats_df, 'features', conn)
    # append_df_to_table(feats_df, 'features')
    # print(f"{feats_df=}")
    # feats_df.to_csv("tmp_feats.csv")
    return
