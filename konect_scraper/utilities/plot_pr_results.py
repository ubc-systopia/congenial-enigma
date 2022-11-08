import argparse
import os
from pathlib import Path

from konect_scraper import config, column_names
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os import listdir
from os.path import isfile, join
from dateutil import parser
from konect_scraper.sql import read_table_as_table
from konect_scraper.util import get_n, get_m, get_pr_struct_size, convert_size, single_val_get, get_additional_stats, \
    get_stat_dtype, get_n_m


def plot_pr_results(graph_name, df):
    num_nodes, num_edges = get_n_m(graph_name)
    vorder_strs = sorted(df.vertex_order.unique())
    eorder_strs = sorted(df.edge_order.unique())
    vorder_map = config.settings['orderings']
    vorder_map['orig'] = 'original'
    eorder_map = config.settings['edge_orderings']

    pr_struct_size = get_pr_struct_size(graph_name)
    size, size_str = convert_size(pr_struct_size)

    results_dir = config.settings['results_dir']

    # make results dir if note exists    
    Path(results_dir).mkdir(parents=True, exist_ok=True)

    dpi = config.settings['plot']['dpi']
    plot_format = config.settings['plot']['format']
    n = len(vorder_strs) * len(eorder_strs)
    w = .25
    x = np.arange(0, len(vorder_strs))
    means = []
    stds = []
    for eorder_str in eorder_strs:
        mns = []
        sds = []
        for vorder_str in vorder_strs:
            slc = df.loc[(
                    (df['edge_order'] == eorder_str) &
                    (df['vertex_order'] == vorder_str)
            )].runtime.values
            mn = slc.mean()
            std = slc.std()
            mns.append(mn)
            sds.append(std)
        means.append(mns)
        stds.append(sds)
    for i, value in enumerate(means):
        std = stds[i]
        eorder_str = eorder_strs[i]
        position = x + (w * (1 - n) / 2) + i * w
        plt.bar(position, value, yerr=std, width=w, label=eorder_str)

    vorder_strs = [vorder_map[i] for i in vorder_strs]
    offset = 3.0
    plt.xticks(x - offset, vorder_strs, rotation=45)
    plt.ylabel("Runtime (ms)")
    plt.legend()
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.rc('text.latex', preamble=r'\usepackage{amssymb} \usepackage{amsmath}')
    stat_cols = ['fill']
    stat_dict = get_additional_stats(stat_cols, graph_name)
    title = r"{} - ".format(graph_name) + \
            r"{} nodes, {} edges, $\approx$ {} {}, ".format(
                "{:,}".format(num_nodes),
                "{:,}".format(num_edges),
                size,
                size_str
            )

    for k, v in stat_dict.items():
        dtype = get_stat_dtype(k)
        title += r"{} - ".format(k)
        if dtype == "Float64":
            title += r"{}".format(f"{v:.4E}")
    plt.title(title)
    plot_path = os.path.join(results_dir, f"{graph_name}.{plot_format}")
    plt.tight_layout()

    plt.savefig(plot_path, dpi=dpi)
    plt.close()
    return


def main(args):
    settings = config.settings
    pr_df = read_table_as_table('pr_expts')
    pr_df['datetime'] = pd.to_datetime(pr_df['datetime'], format='%d-%m-%Y %H-%M-%S')

    filter_datetime = args.min_date  # only look at experiments after this date
    pr_df = pr_df[
        (pr_df['datetime'] >= pd.to_datetime(filter_datetime, format='%d-%m-%Y %H-%M-%S'))
    ]
    pr_df['pr_struct_size'] = pr_df.apply(
        lambda row: get_pr_struct_size(row['graph_name']),
        axis=1
    ).astype(np.int64)
    pr_df = pr_df.sort_values('pr_struct_size', ascending=True)

    vorders = settings['orderings']
    eorders = settings['edge_orderings']
    n_expts = settings['hyperparameters']['pr-experiments']['num_expts']
    reqd_num_expts = (len(vorders) + 1) * len(eorders) * n_expts

    gdfs = [[name, gdf] for name, gdf in pr_df.groupby('graph_name')]

    gdfs = sorted(gdfs, key=lambda x: get_pr_struct_size(x[0]), reverse=True)


    for name, gdf in gdfs:
        # if gdf.shape[0] != reqd_num_expts:
        #     print(name)
        #     continue
        plot_pr_results(name, gdf)
    return


if __name__ == '__main__':
    argparse_desc = """
    A module that plots the results of PageRank experiments using a 
    multicoloured bar plot
    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--data-dir', required=True,
                        help="Path to directory that should store sqlite3.db")
    parser.add_argument('-m', '--min-date', required=True,
                        help="Only plot experimental results after this date. "
                        "Expected format: 'd-m-Y H-M-S' "
                        "e.g. 18-10-2022 15-16-45" )
    args = parser.parse_args()
    config.init(args.data_dir)
    column_names.init()
    main(args)
