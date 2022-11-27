import logging
import os
import subprocess
from konect_scraper.util import create_dir_if_not_exists

import konect_scraper.config as config
from konect_scraper.util import get_directed, get_n, get_m
from konect_scraper.sql import append_df_to_table
import pandas as pd
import datetime


def run_pr_expt(graph_name, order_str, edge_order_str):
    settings = config.settings
    debug = settings['debug']
    directed = get_directed(graph_name)
    n = get_n(graph_name)
    m = get_m(graph_name)
    num_iters = settings['hyperparameters']['pr-experiments']['num_iters']
    num_expts = settings['hyperparameters']['pr-experiments']['num_expts']
    damping_factor = settings['hyperparameters']['pr-experiments']['damping_factor']
    graphs_dir = settings['graphs_dir']
    graph_dir = os.path.join(graphs_dir, graph_name)
    results_dir = os.path.join(settings['results_dir'], graph_name)
    create_dir_if_not_exists(results_dir)

    date_str = datetime.datetime.utcnow().strftime("%a_%b_%d_%H_%M_%S_UTC_%Y")
    results_path = os.path.join(
        results_dir,
        f'{order_str}_{edge_order_str}_{date_str}.csv'
    )
    sqlite3_db_path = settings['sqlite3']['sqlite3_db_path']
    pr_experiments_executable = settings['pr_experiments_executable']
    args = [pr_experiments_executable]
    if directed:
        args += ['-d']
    if debug:
        args += ['-e']
    args += [
        '-n', str(n),
        '-m', str(m),
        '-i', str(num_iters),
        '-x', str(num_expts),
        '-p', str(damping_factor),
        '-g', str(graph_dir),
        '-b', str(sqlite3_db_path),
        '-o', str(order_str),
        '-r', str(results_path),
        '-s', edge_order_str
    ]

    logging.info(f"Executing: " + ' '.join(args))

    res = subprocess.check_output(args)

    # after execution, read the results csv and append to sqlite3 table
    # append_df_to_table(
    #     pd.read_csv(results_path),
    #     'pr_expts'
    # )
    # appending results of each experiment may fail due to locking of sqlite3
    # table instead, after completion of _all_ pr-expts submitted using sbatch
    # arr call:
    # `python konect_scraper.utilities.append_results_to_db.py \
    #    --results-dir {path to results} \
    #    --data-dir {path to dir with graphs.db}`
    # to append all results contained in results dir, and then erase them


    return


def main(rows, orders):
    """"""
    settings = config.settings
    edge_orders = settings['edge_orderings']
    # compute the given orders for each of the datasets
    for row in rows:
        graph_name = row['graph_name']
        for vertex_order in orders:
            for edge_order in edge_orders:
                logging.info(
                    f"Running PageRank experiments on {graph_name}-{vertex_order}-{edge_order}..")
                run_pr_expt(graph_name, vertex_order, edge_order)

    return
