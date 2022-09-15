import logging
import os
import subprocess

import konect_scraper.config as config
from konect_scraper.util import get_directed, get_n, get_m


def run_pr_expt(graph_name, order_str):
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
        '-o', str(order_str)
    ]

    logging.info(f"Executing: " + ' '.join(args))

    res = subprocess.check_output(args)

    return


def main(rows, orders):
    """"""
    settings = config.settings

    # compute the given orders for each of the datasets
    for row in rows:
        graph_name = row['graph_name']
        for order in orders:
            logging.info(f"Running PageRank experiments on {graph_name}-{order}..")
            run_pr_expt(graph_name, order)

    return


if __name__ == '__main__':
    main()
