import logging
import os
import shutil
import subprocess
from pathlib import Path

from konect_scraper import config
from konect_scraper.util import get_directed, get_n, get_m


def compute_dbg_order(graph_name, order_str):
    settings = config.settings
    order_str_dict = settings['dbg']['order_str_dict']
    order_idx_dict = settings['dbg']['order_idx_dict']
    degree_used_for_reordering = settings['dbg']['degree_used_for_reordering']
    dbg_order_idx = order_idx_dict[order_str_dict[order_str]]
    dbg_apps_dir = settings['dbg_apps_dir']
    dbg_home = settings['dbg_home']
    dbg_datasets_dir = settings['dbg_datasets_dir']
    make_executable = settings['make_executable']
    max_iters = settings['dbg']['max_iters']
    sqlite3_db_path = settings['sqlite3']['sqlite3_db_path']
    order_file = os.path.join(dbg_datasets_dir, f"{graph_name}.{dbg_order_idx}.map")

    graphs_dir = settings['graphs_dir']
    graph_dir = os.path.join(graphs_dir, graph_name)
    order_file = os.path.join(graph_dir, f"{dbg_order_idx}.map")
    # order_file = os.path.join(dbg_datasets_dir, f"{graph_name}.{dbg_order_idx}.map")

    n = get_n(graph_name)
    m = get_m(graph_name)
    args = [
        make_executable,
        f"REORDERING_ALGO={dbg_order_idx}",
        f"DEGREE_USED_FOR_REORDERING={degree_used_for_reordering}",
        f"DATASET={graph_name}",
        f"MAXITERS={max_iters}",
        f"ORDER_FILE={order_file}",
        f"NUM_VERTICES={n}",
        f"NUM_EDGES={m}",
        f"GRAPH_NAME={graph_name}",
        f"SQLITE_DB_PATH={sqlite3_db_path}",
        f"GPATH={graph_dir}",
        f"DBG_ROOT={dbg_home}",
        f"COMPRESSED_EDGE_LIST_FNAME={settings['compressed_el_file_name']}",
        "run-PageRank",
    ]
    env = {
        **os.environ,
        # "DBG_ROOT": dbg_home,
    }
    print(f'{dbg_home=}')
    print(f'{dbg_apps_dir=}')
    logging.info(f"Executing: " + ' '.join(args))
    print(f"Executing: " + ' '.join(args))
    res = subprocess.check_output(args, cwd=dbg_apps_dir, env=env)
    print(res.decode('utf-8'))

    # copy the map from dbg dataset dir
    # reprocess map file to match format used by app
    shutil.copy(
        order_file,
        os.path.join(
            settings['graphs_dir'],
            graph_name,
            order_str
        )
    )

    # remove the duplicate file from the dbg directory
    try:
        os.remove(order_file)
    except OSError:
        pass

    # print(' '.join(args))

    return


def compute_rabbit(graph_path, order_path, graph_name, comms_path):
    settings = config.settings
    sqlite3_db_path = settings['sqlite3']['sqlite3_db_path']

    executable = settings['rabbit_order_executable']
    args = [executable]
    args += [graph_path, order_path, graph_name, sqlite3_db_path, comms_path]

    logging.info(f"Executing: " + ' '.join(args))

    res = subprocess.check_output(args)

    return

def compute_par_slashburn(graph_path, order_path):
    settings = config.settings
    percent = settings['hyperparameters']['slashburn']['percent']
    sqlite3_db_path = settings['sqlite3']['sqlite3_db_path']

    executable = settings['par_slashburn_executable']
    args = [executable]
    args += [
        '-f', graph_path,
        '-s',
        '-o', order_path,
        '-d', sqlite3_db_path,
        '-p', str(percent)
    ]
    logging.info(f"Executing: " + ' '.join(args))

    res = subprocess.check_output(args)
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


def compute_parallel_batch_cm(graph_path, n, m, directed, cm_path, rcm_path):
    settings = config.settings
    percent = settings['hyperparameters']['slashburn']['percent']
    sqlite3_db_path = settings['sqlite3']['sqlite3_db_path']

    executable = settings['parallel_batch_rcm_executable']
    args = [executable]

    args += [
        graph_path,
        '-i', 'CPU_BATCH',
        '-t', str(config.settings['n_threads']),
        '-o', cm_path,
        '-n', str(n),
        '-m', str(m),
        '--rcm', rcm_path,
        '--sqlite3_db', sqlite3_db_path,
    ]
    if directed:
        args += ['-d']
    # print(f"Executing: " + ' '.join(args))

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


def compute_ordering(graph_name, order, ovewrite):
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

    if Path(order_path).is_file() and not ovewrite:  # if already computed, skip
        logging.info(f"{graph_name}-{order_str} already computed; skipping.")
        return

    match order:

        case 'rnd' | 'srt' | 'hs' | 'hc' | 'dbg':
            compute_dbg_order(graph_name, order)

        case "rbt":
            extension = ".net"

            comp_graph_path = os.path.join(graph_dir, settings['compressed_el_file_name'] + extension)
            comms_path = os.path.join(graph_dir, "comms")
            compute_rabbit(comp_graph_path, order_path, graph_name, comms_path)

        case "parsb":
            comp_graph_path = os.path.join(graph_dir, settings['compressed_el_file_name'] + ".net")
            compute_par_slashburn(comp_graph_path, order_path)

        case "sb":
            # comp_graph_path = os.path.join(graph_dir, settings['compressed_el_file_name'] + ".net")
            # compute_par_slashburn(comp_graph_path, order_path)
            compute_slashburn(comp_graph_path, order_path, directed, n, m)

        case "cm":
            extension = ".net"

            comp_graph_path = os.path.join(graph_dir, settings['compressed_el_file_name'] + extension)
            cm_order_path = order_path
            rcm_order_path = os.path.join(graph_dir, 'rev_cm')
            # compute_cuthill_mckee(comp_graph_path, n, m)
            compute_parallel_batch_cm(comp_graph_path, n, m, directed, cm_order_path, rcm_order_path)
        case "rev_cm":
            return

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
