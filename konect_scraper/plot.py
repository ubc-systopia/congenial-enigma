import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess
import os
import gc

import multiprocessing

from matplotlib import patches, cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as plticker
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import sys
from PyQt5.QtWidgets import QApplication
from konect_scraper.io import get_adj_mat_from_edge_list, read_iso_map, save_spy_plots, read_image
from konect_scraper import config
from konect_scraper.util import get_directed, get_n, get_m, create_plot_dirs_if_not_exists, get_graph_dir, get_plot_dir, \
    translate_adj_mat, chunks, single_val_get, get_critical_depth, next_largest_multiple
import matplotlib.colors


def __init__(self):
    config.init()


def adj_mat_spy(canvas, arr, path):
    settings = config.settings
    marker = settings['plot']['marker']
    markersize = settings['plot']['markersize']
    fmt = settings['plot']['format']
    dpi = settings['plot']['dpi']
    canvas.spy(arr, marker=marker, markersize=markersize)
    canvas.savefig(f"{path}.{fmt}", dpi=dpi)
    plt.cla()


def adj_matshow(canvas, arr, path):
    settings = config.settings
    fmt = settings['plot']['format']
    dpi = settings['plot']['dpi']
    canvas.matshow(arr)
    canvas.savefig(f"{path}.{fmt}", dpi=dpi)
    plt.cla()


def plot_edge_list_as_adj_mat(canvas, graph_name, directed, label_str, el_file_name):
    settings = config.settings

    graph_dir = get_graph_dir(graph_name)
    mat_plot_dir = get_plot_dir(graph_name, "adj_mat")
    spy_plot_dir = get_plot_dir(graph_name, "spy")
    graph_path = os.path.join(graph_dir, settings[el_file_name])

    adj_mat = get_adj_mat_from_edge_list(graph_path, directed)
    logging.info(f"Plotting {label_str} {graph_name}'s Adjacency Matrix..")
    adj_matshow(canvas, adj_mat, os.path.join(mat_plot_dir, label_str))
    logging.info(f"Plotting {label_str} {graph_name}'s Spy Plot..")
    adj_mat_spy(canvas, adj_mat, os.path.join(spy_plot_dir, label_str))


def plot_compressed(canvas, graph_name, directed):
    plot_edge_list_as_adj_mat(canvas, graph_name, directed, "comp", "compressed_el_file_name")

    return


def plot_orig(canvas, graph_name, directed):
    plot_edge_list_as_adj_mat(canvas, graph_name, directed, "orig", "orig_el_file_name")
    return


# ax_plot_fn(ax, graph_name, directed, plot_type, label_str, markersize, label_str, el_file_name, )


def plot_ordering(canvas, graph_name, directed, order):
    settings = config.settings

    orderings = settings['orderings']
    graphs_dir = settings['graphs_dir']
    graph_dir = os.path.join(graphs_dir, graph_name)

    mat_plot_dir = get_plot_dir(graph_name, "adj_mat")
    spy_plot_dir = get_plot_dir(graph_name, "spy")
    graph_path = os.path.join(graph_dir, settings["compressed_el_file_name"])

    adj_mat = get_adj_mat_from_edge_list(graph_path, directed)
    iso_path = os.path.join(graph_dir, order)
    iso_map = read_iso_map(iso_path)

    map_mat = translate_adj_mat(adj_mat, iso_map)

    logging.info(f"Plotting {order} {graph_name}'s Adjacency Matrix..")
    adj_matshow(canvas, map_mat, os.path.join(mat_plot_dir, order))
    logging.info(f"Plotting {order} {graph_name}'s Spy Plot..")
    adj_mat_spy(canvas, map_mat, os.path.join(spy_plot_dir, order))


def get_markersize(n):
    if n > 10000:
        return 1

    if n > 2000:
        return 1

    return 1


def split_aggregate_plots(graph_names, orders, fmt):
    settings = config.settings
    max_rows_per_agg_spy_plot = settings['plot']['max_rows_per_agg_spy_plot']
    plot_num = 0
    plots_dir = settings['plots_dir']

    if len(graph_names) > max_rows_per_agg_spy_plot:
        chnks = chunks(graph_names, max_rows_per_agg_spy_plot)
        for graph_chunk in chnks:
            all_spy_path = os.path.join(plots_dir, f"all_spy_{plot_num}.{fmt}")
            aggregate_plots(graph_chunk, orders, all_spy_path)
            plot_num += 1
    else:
        all_spy_path = os.path.join(plots_dir, f"all_spy_{plot_num}.{fmt}")
        aggregate_plots(graph_names, orders, all_spy_path)


def aggregate_plots(graph_names, orders, all_spy_path):
    settings = config.settings
    orderings = settings['orderings']
    bbox_inches = settings['plot']['bbox_inches']
    pad_inches = settings['plot']['pad_inches']
    order_names = ["Original", "Compressed"] + [orderings[o] for o in orders]
    adj_mat_format = settings['plot']['adj_mat_format']
    plots_dir = settings['plots_dir']
    fmt = settings['plot']['format']
    dpi = settings['plot']['dpi']
    ax_size = settings['plot']['ax_size']
    n_rows = len(graph_names)
    n_cols = len(orders) + 2  # include the original and compressed
    fig_size = (n_cols * ax_size, n_rows * ax_size,)

    image_paths = np.zeros((n_rows, n_cols), dtype=object)

    image_paths[:, 0] = [
        os.path.join(plots_dir, graph_name, adj_mat_format, f'orig.{fmt}') for graph_name in graph_names
    ]

    image_paths[:, 1] = [
        os.path.join(plots_dir, graph_name, adj_mat_format, f'comp.{fmt}') for graph_name in graph_names
    ]

    for row_idx, graph_name in enumerate(graph_names):
        for ord_idx, order in enumerate(orders):
            col_idx = ord_idx + 2
            image_paths[row_idx][col_idx] = os.path.join(
                plots_dir, graph_name, adj_mat_format, f'{order}.{fmt}')

    fig, axs = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=fig_size)

    for row_idx, graph_name in enumerate(graph_names):
        n = get_n(graph_name)
        for col_idx in range(n_cols):
            im_path = image_paths[row_idx, col_idx]
            if len(axs.shape) == 1:  # 1d, reshape to 2d
                axs = axs.reshape((1, axs.shape[0]))
            ax = axs[row_idx, col_idx]

            img = read_image(im_path, fmt)
            ax.imshow(img, interpolation='none', extent=[0, n - 1, n - 1, 0])
            ax.xaxis.tick_top()

    for ax, col in zip(axs[0], order_names):
        ax.set_title(col)

    for ax, row in zip(axs[:, 0], graph_names):
        ax.set_ylabel(row, rotation=45, size='large')
        # graph names in y label of axes overwrites some figures - adjust it
        ax.yaxis.set_label_coords(-.2, 0.5)
    fig.tight_layout()

    plt.autoscale()

    fig.savefig(all_spy_path, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)
    return


def ax_plot_order(ax, graph_name, directed, plot_type, markersize, order, ):
    settings = config.settings
    graphs_dir = settings['graphs_dir']

    graph_dir = os.path.join(graphs_dir, graph_name)
    graph_path = os.path.join(graph_dir, settings["compressed_el_file_name"])

    graph_path += ".net"

    adj_mat = get_adj_mat_from_edge_list(graph_path, directed)
    iso_path = os.path.join(graph_dir, order)
    iso_map = read_iso_map(iso_path)

    map_mat = translate_adj_mat(adj_mat, iso_map)
    if plot_type == "adj_mat":
        ax.matshow(map_mat)
    else:
        ax.spy(map_mat, markersize=markersize)

    # highlight the dense region defined by sb_k x sb_num_iters
    if order == 'sb':
        sb_k = single_val_get('sb_k', 'statistics', graph_name)
        sb_n_iters = single_val_get('sb_n_iters', 'statistics', graph_name)
        if sb_k and sb_n_iters:
            sqr_width = sb_k * sb_n_iters
            rect = patches.Rectangle((0, sqr_width), sqr_width, -sqr_width,
                                     linewidth=1, edgecolor='r', alpha=0.3, zorder=2, facecolor='r')
            ax.add_patch(rect)

    if order == 'rbt':
        # read the original vertex id community assignment
        comms_path = os.path.join(graph_dir, "comms")
        comms = np.loadtxt(comms_path).astype(np.uint32)
        mapped_comms = np.zeros(adj_mat.shape[0])
        np.set_printoptions(threshold=sys.maxsize)
        for i, c in enumerate(comms):
            mapped_comms[iso_map[i]] = c

        plot_rbt_comms(ax, mapped_comms)

    return

def plot_rbt_comms(ax, comms):
    n_comms = len(np.unique(comms))
    color = iter(cm.rainbow(np.linspace(0, 1, n_comms)))
    for comm in np.unique(comms):

        start = np.argwhere(comms == comm)[0][0]
        end = np.argwhere(comms == comm)[-1][0]
        comm_size = end - start

        # plot the updated region size
        critical_depth = get_critical_depth()

        comm_size = next_largest_multiple(comm_size, critical_depth + 2)

        c = next(color)
        rect = patches.Rectangle((start, end), comm_size, -comm_size,
                                 linewidth=1, edgecolor=c, alpha=0.3, zorder=2, facecolor=c)
        ax.add_patch(rect)

def ax_plot_adj_mat(ax, graph_name, directed, plot_type, markersize, label_str, el_file_name):
    settings = config.settings

    graph_dir = get_graph_dir(graph_name)
    graph_path = os.path.join(graph_dir, settings[el_file_name])

    if settings[el_file_name] == "comp":
        graph_path += ".net"

    adj_mat = get_adj_mat_from_edge_list(graph_path, directed)
    if plot_type == "adj_mat":
        ax.matshow(adj_mat)
    else:
        ax.spy(adj_mat, markersize=markersize)
    return


def ax_plot_and_clear(graph_name, directed, label_str, el_file_name, fmt, markersize, figsize, plots_dir, dpi,
                      ax_plot_fn, plot_type):
    settings = config.settings
    ax_path = os.path.join(plots_dir, graph_name, plot_type, f"{label_str}.{fmt}")

    if Path(ax_path).is_file():  # if already computed, skip
        logging.info(f"{graph_name}-{label_str} already plotted; skipping.")
        return
    fig, ax = plt.subplots(figsize=figsize)
    logging.info(f"Plotting {graph_name}-{label_str}")

    bbox_inches = settings['plot']['bbox_inches']
    pad_inches = settings['plot']['pad_inches']

    match ax_plot_fn.__name__:
        case 'ax_plot_order':
            ax_plot_fn(ax, graph_name, directed, plot_type, markersize, label_str, )
        case 'ax_plot_adj_mat':
            ax_plot_fn(ax, graph_name, directed, plot_type, markersize, label_str, el_file_name, )
    ax.axis('off')
    plt.tight_layout()
    fig.savefig(ax_path, bbox_inches=bbox_inches, pad_inches=pad_inches, dpi=dpi)
    fig.clear()
    plt.close(fig)


def get_figsize(ax_size, n):
    if n < 5_000:
        return (ax_size, ax_size)
    side_len = int(n * 0.000_2) * ax_size

    return side_len, side_len


def plot_edge_orderings(rows, vorder_str):
    edge_orderings = config.settings['hyperparameters']['pr-experiments']['edge_orderings']
    for row in rows:
        graph_name = row['graph_name']
        for eorder_str in edge_orderings:
            plot_edge_ordering(graph_name, vorder_str, eorder_str)

    return


def plot_edge_ordering(graph_name, vorder_str, eorder_str):
    # colour map the adjacency matrix

    fig, ax = plt.subplots()

    n = get_n(graph_name)
    graphs_dir = config.settings['graphs_dir']
    plots_dir = config.settings['plots_dir']
    graph_dir = os.path.join(graphs_dir, graph_name)
    plot_dir = os.path.join(plots_dir, graph_name)
    plot_format = config.settings['plot']['format']
    dpi = config.settings['plot']['dpi']
    adj_mat_format = config.settings['plot']['adj_mat_format']
    edgelist_path = os.path.join(graph_dir, f"{vorder_str}.{eorder_str}")
    plot_path = os.path.join(plot_dir, f"{vorder_str}_{eorder_str}.{plot_format}")
    print(plot_path)

    sorted_edges = np.loadtxt(edgelist_path).astype(np.uint32)
    adj_mat = np.zeros((n, n))
    c = np.zeros(sorted_edges.shape[0])
    eid = 1
    for e in sorted_edges:
        adj_mat[e[0], e[1]] = eid
        c[eid - 1] = eid;

        eid += 1

    adj_mat = np.ma.masked_where(adj_mat == 0, adj_mat)

    cmap = mpl.cm.get_cmap("plasma").copy()
    # cmap.set_bad(color='white')
    # norm = plt.Normalize(1, eid)
    # plt.spy(adj_mat, cmap=cmap, interpolation='none')
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.scatter(sorted_edges[:, 1], sorted_edges[:, 0], c=c, cmap=cmap, marker='.')
    # ax.set_extent([0, n - 1, n - 1, 0])
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    plt.colorbar()
    # fig.colorbar(im, cax=cax)
    plt.savefig(plot_path, figsize=(5, 5), dpi=dpi, )
    plt.close()

    return

def plot_graph(graph_name, orders, settings):
    logging.info(f"Plotting {graph_name}..")
    create_plot_dirs_if_not_exists(graph_name)

    directed = bool(get_directed(graph_name))
    m = get_m(graph_name)

    n = get_n(graph_name)
    ax_size = settings['plot']['ax_size']
    max_n = settings['plot']['max_n']
    dpi = settings['plot']['dpi']

    # if n > max_n:
    #     logging.error(f"Not plotting {graph_name} - {n} is too many nodes!")
    # plot_orig(plt, graph_name, directed)
    # plot_compressed(plt, graph_name, directed)
    markersize = get_markersize(n)
    plots_dir = settings['plots_dir']
    plot_type = "spy"

    figsize = get_figsize(ax_size, n)
    fmt = settings['plot']['format']

    # plot the original adjacency matrix
    ax_plot_and_clear(graph_name, directed, "orig", "orig_el_file_name", fmt, markersize, figsize, plots_dir, dpi,
                        ax_plot_adj_mat, plot_type)
    # plot the compressed adjacency matrix
    ax_plot_and_clear(graph_name, directed, "comp", "compressed_el_file_name", fmt, markersize, figsize, plots_dir,
                        dpi,
                        ax_plot_adj_mat, plot_type)

    # adj_mat_fig, adj_mat_axs = plt.subplots(nrows, ncols, figsize=figsize)
    # spy_fig, spy_axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    # all_axs = [spy_axs]
    plot_types = ["adj_mat", "spy"]
    plot_types = ["spy"]

    for order_idx, order in enumerate(orders):
        # plot_ordering(plt, graph_name, directed, order)
        for plot_type in plot_types:
            ax_plot_and_clear(graph_name, directed, order, None, fmt, markersize, figsize,
                                plots_dir,
                                dpi,
                                ax_plot_order, plot_type)

    gc.collect()

def main(rows, orders):
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)
    plot_module = sys.modules[__name__]
    settings = config.settings
    adj_mat_format = settings['plot']['adj_mat_format']
    dpi = settings['plot']['dpi']
    if adj_mat_format == "spy":
        plot_fn = getattr(plot_module, "adj_mat_spy")
    else:
        plot_fn = getattr(plot_module, "adj_matshow")
    fmt = settings['plot']['format']

    nrows = len(rows)
    ncols = len(orders) + 2  # plot all orders + original and compressed isomorphisms
    ax_size = settings['plot']['ax_size']
    max_n = settings['plot']['max_n']
    figsize = (ncols * ax_size, nrows * ax_size,)

    from multiprocessing import Process
    procs = []
    # compute the given orders for each of the datasets
    for row_idx, row in enumerate(rows):

        graph_name = row['graph_name']
        args = [graph_name, orders, settings]
        args = {
            'graph_name': graph_name,
            'orders': orders,
            'settings': settings
        }
        print(args)
        proc = Process(target=plot_graph, kwargs=args)
        procs.append(proc)
        proc.start()
        # pool.map(plot_graph, args=(args,))
        # plot_graph(graph_name, orders, settings)
    # complete the processes
    for proc in procs:
        proc.join()

    graph_names = [row['graph_name'] for row in rows]

    logging.info(f"Plotting an aggregate plots of {len(graph_names)} graphs and {len(orders)} orders..")
    split_aggregate_plots(graph_names, orders, fmt)
    return


if __name__ == '__main__':
    main()
