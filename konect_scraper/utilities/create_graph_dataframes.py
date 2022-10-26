
import os
from .. import config, column_names
from ..sql import (get_all_unipartite_directed_graphs,
                   get_all_unipartite_undirected_graphs,
                   get_all_bipartite_graphs,
                   get_all_downloadable_graphs,
                   get_all_graphs_by_graph_names)
from ..util import (get_n_vertices, get_n_edges, convert_size, get_pr_struct_size,
                    get_category, get_directed)

import pandas as pd

"""
Given all graphs on konect.cc, persist them into 3 tables:
Directed (unipartite), Undirected (unipartite), and Bipartite.

Save these tables into a sqlite.db and csv
"""


def populate_df(rows, df):
    for i, row in enumerate(rows):
        size, size_string = convert_size(get_pr_struct_size(row['graph_name']))
        d = {
            'index': i,
            'graph_name': row['graph_name'],
            'num_vertices': get_n_vertices(row['graph_name']),
            'num_edges': get_n_edges(row['graph_name']),
            'pr_struct_size': str(size) + size_string,
            'category': get_category(row['graph_name']),
        }
        df = pd.concat([df, pd.DataFrame([d])], ignore_index=True)
    return df


def sort_downloadable_graphs(rows):
    graph_names = [r['graph_name'] for r in rows]
    rows = get_all_downloadable_graphs(graph_names)
    graph_names = [r['graph_name'] for r in rows]
    rows = get_all_graphs_by_graph_names(graph_names)
    rows = sorted(rows, key=lambda r: get_pr_struct_size(
        r['graph_name']), reverse=False)
    return rows

def main():
    cols = ['index', 'graph_name', 'num_vertices', 'num_edges', 'pr_struct_size',
            'category']

    config.init()
    settings = config.settings

    dir_rows = sort_downloadable_graphs(get_all_unipartite_directed_graphs())
    undir_rows = sort_downloadable_graphs(
        get_all_unipartite_undirected_graphs())
    bip_rows = sort_downloadable_graphs(get_all_bipartite_graphs())

    dir_df = pd.DataFrame(columns=cols)
    undir_df = pd.DataFrame(columns=cols)
    bip_df = pd.DataFrame(columns=cols)

    dir_df = populate_df(dir_rows, dir_df)
    undir_df = populate_df(undir_rows, undir_df)
    bip_df = populate_df(bip_rows, bip_df)

    dfs_dir = settings['dataframes_dir']

    dir_df.to_csv(os.path.join(dfs_dir, 'directed.csv'))
    undir_df.to_csv(os.path.join(dfs_dir, 'undirected.csv'))
    bip_df.to_csv(os.path.join(dfs_dir, 'bipartite.csv'))

    return


if __name__ == '__main__':
    main()
