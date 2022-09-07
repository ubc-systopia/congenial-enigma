import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os
import sys
import gc

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from konect_scraper.io import get_adj_mat_from_edge_list, read_iso_map, save_spy_plots, read_image
from konect_scraper import config
from konect_scraper.util import get_directed, get_n, get_m, create_plot_dirs_if_not_exists, get_graph_dir, get_plot_dir, \
    translate_adj_mat


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


def ax_plot_adj_mat(ax, graph_name, directed, label_str, el_file_name, plot_format, markersize):
    settings = config.settings

    graph_dir = get_graph_dir(graph_name)
    graph_path = os.path.join(graph_dir, settings[el_file_name])

    if settings[el_file_name] == "comp":
        graph_path += ".net"

    adj_mat = get_adj_mat_from_edge_list(graph_path, directed)
    if plot_format == "adj_mat":
        ax.matshow(adj_mat)
    else:
        ax.spy(adj_mat, markersize=markersize)


def plot_compressed(canvas, graph_name, directed):
    plot_edge_list_as_adj_mat(canvas, graph_name, directed, "comp", "compressed_el_file_name")

    return


def plot_orig(canvas, graph_name, directed):
    plot_edge_list_as_adj_mat(canvas, graph_name, directed, "orig", "orig_el_file_name")
    return


def ax_plot_order(ax, graph_name, directed, order, plot_format, markersize):
    settings = config.settings
    graphs_dir = settings['graphs_dir']

    graph_dir = os.path.join(graphs_dir, graph_name)
    graph_path = os.path.join(graph_dir, settings["compressed_el_file_name"])

    graph_path += ".net"

    adj_mat = get_adj_mat_from_edge_list(graph_path, directed)
    iso_path = os.path.join(graph_dir, order)
    iso_map = read_iso_map(iso_path)

    map_mat = translate_adj_mat(adj_mat, iso_map)
    if plot_format == "adj_mat":
        ax.matshow(map_mat)
    else:
        ax.spy(map_mat, markersize=markersize)

    return


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
        return 3

    if n > 2000:
        return 1

    return 1


def aggregate_plots(graph_names, orders):
    settings = config.settings
    orderings = settings['orderings']

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

    for row_idx in range(n_rows):
        for col_idx in range(n_cols):
            im_path = image_paths[row_idx, col_idx]
            ax = axs[row_idx, col_idx]
            img = read_image(im_path, fmt)
            ax.imshow(img)

    all_spy_path = os.path.join(plots_dir, f"all_spy.{fmt}")

    for ax, col in zip(axs[0], order_names):
        ax.set_title(col)

    for ax, row in zip(axs[:, 0], graph_names):
        ax.set_ylabel(row, rotation=45, size='large')
    fig.tight_layout()
    fig.savefig(all_spy_path, )
    return


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

    nrows = len(rows)
    ncols = len(orders) + 2  # plot all orders + original and compressed isomorphisms
    ax_size = settings['plot']['ax_size']
    figsize = (ncols * ax_size, nrows * ax_size,)

    adj_mat_fig, adj_mat_axs = plt.subplots(nrows, ncols, figsize=figsize)
    spy_fig, spy_axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    all_axs = [spy_axs]
    plot_formats = ["adj_mat", "spy"]
    plot_formats = ["spy"]
    # compute the given orders for each of the datasets
    for row_idx, row in enumerate(rows):

        graph_name = row['graph_name']
        logging.info(f"Plotting {graph_name}..")
        create_plot_dirs_if_not_exists(graph_name)

        directed = bool(get_directed(graph_name))
        m = get_m(graph_name)

        n = get_n(graph_name)

        if n > 10_000:
            logging.error(f"Not plotting {graph_name} - {n} is too many nodes!")
        # plot_orig(plt, graph_name, directed)
        # plot_compressed(plt, graph_name, directed)
        markersize = get_markersize(n)
        plots_dir = settings['plots_dir']
        plot_format = "spy"
        # plot the original adjacency matrix
        figsize = (ax_size, ax_size)
        fig, ax = plt.subplots(figsize=figsize)
        ax_plot_adj_mat(ax, graph_name, directed, "orig", "orig_el_file_name", plot_format, markersize)
        fmt = settings['plot']['format']
        ax_path = os.path.join(plots_dir, graph_name, plot_format, f"orig.{fmt}")
        ax.axis('off')
        plt.tight_layout()
        fig.savefig(ax_path, dpi=dpi)
        fig.clear()
        plt.close(fig)

        # plot the compressed adjacency matrix
        fig, ax = plt.subplots(figsize=figsize)
        ax_plot_adj_mat(ax, graph_name, directed, "comp", "compressed_el_file_name", plot_format, markersize)
        ax_path = os.path.join(plots_dir, graph_name, plot_format, f"comp.{fmt}")
        ax.axis('off')
        plt.tight_layout()
        fig.savefig(ax_path, dpi=dpi)
        fig.clear()

        plt.close(fig)


        for order_idx, order in enumerate(orders):
            # plot_ordering(plt, graph_name, directed, order)
            for plot_format, axs in zip(plot_formats, all_axs):
                fig, ax = plt.subplots(figsize=figsize)

                ax_plot_order(ax, graph_name, directed, order, plot_format, markersize)
                ax_path = os.path.join(plots_dir, graph_name, plot_format, f"{order}.{fmt}")
                ax.axis('off')
                plt.tight_layout()
                fig.savefig(ax_path, dpi=dpi)
                fig.clear()

                plt.close(fig)
        gc.collect()

    graph_names = [row['graph_name'] for row in rows]

    logging.info(f"Plotting an aggregate plots of {len(graph_names)} graphs and {len(orders)} orders..")
    aggregate_plots(graph_names, orders)

    # all_spy_path = os.path.join(plots_dir, f"all_spy.{fmt}")
    # save_spy_plots(spy_fig, spy_axs, graph_names, order_names, all_spy_path)
    return


if __name__ == '__main__':
    main()
