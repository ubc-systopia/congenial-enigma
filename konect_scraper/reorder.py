from pathlib import Path

from konect_scraper import config
import os
import logging
import subprocess

from konect_scraper.config import IOMode
from konect_scraper.util import get_directed, get_n, get_m


def compute_random(graph_path):
    return

def compute_rabbit(graph_path):
    return

def compute_slashburn(graph_path, order_path, directed, n, m):
    settings = config.settings
    percent = settings['hyperparameters']['slashburn']['percent']
    sqlite3_db_path = settings['sqlite3']['sqlite3_db_path']

    executable = settings['slashburn_executable']
    args = [executable]

    if directed:
        args += ['-d']

    # for io_mode in io_modes:
    #     match io_mode:
    #         case IOMode.text:
    #             args += ['-t']
    #         case IOMode.binary:
    #             args += ['-i']

    args += [
        '-n', str(n),
        '-m', str(m),
        '-p', str(percent),
        '-g', graph_path,
        '-b', sqlite3_db_path,
        '-o', order_path,
    ]

    logging.info(f"Executing: " + ' '.join(args))

    res = subprocess.check_output(args)

    return
def compute_cuthill_mckee(graph_path, n, m):
    settings = config.settings
    percent = settings['hyperparameters']['slashburn']['percent']
    sqlite3_db_path = settings['sqlite3']['sqlite3_db_path']

    executable = settings['cuthill_mckee_executable']
    args = [executable]

    args += [
        '-n', str(n),
        '-m', str(m),
        '-g', graph_path,
        '-b', sqlite3_db_path,
    ]

    logging.info(f"Executing: " + ' '.join(args))

    res = subprocess.check_output(args)

    return
def compute_ordering(graph_name, order):
    settings = config.settings

    orderings = settings['orderings']
    graphs_dir = settings['graphs_dir']
    graph_dir = os.path.join(graphs_dir, graph_name)

    extension = ".bin"

    comp_graph_path = os.path.join(graph_dir, settings['compressed_el_file_name'] + extension)
    directed = bool(get_directed(graph_name))
    n = get_n(graph_name)
    m = get_m(graph_name)

    order_str = orderings[order]

    order_path = os.path.join(graph_dir, order)

    if Path(order_path).is_file():  # if already computed, skip
        logging.info(f"{graph_name}-{order_str} already computed; skipping.")
        return
    match order:
        case "rnd":
            compute_random(comp_graph_path)
        case "rbt":
            compute_rabbit(comp_graph_path)
        case "sb":
            compute_slashburn(comp_graph_path, order_path, directed, n, m)
        case "cm":
            compute_cuthill_mckee(comp_graph_path, n, m)
        case _:
            logging.error(f"{order}: Unsupported Ordering!")
        # case ""

def main(rows, orders):
    settings = config.settings


    # compute the given orders for each of the datasets
    for row in rows:
        graph_name = row['graph_name']

        for order in orders:
            compute_ordering(graph_name, order)
    return

if __name__ == '__main__':
    main()
