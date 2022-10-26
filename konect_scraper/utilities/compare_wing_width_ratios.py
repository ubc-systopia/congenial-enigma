import sqlite3

import pandas as pd

from konect_scraper import config


def main():
    """
    After computing the slashburn/parslashburn orderings for each graph,
    identify the graphs with lowest wing width ratio
    Returns:

    """
    cnx = sqlite3.connect(config.settings['sqlite3']['sqlite3_db_path'])

    df = pd.read_sql_query(f"SELECT * FROM 'statistics'", cnx)
    nm_df = pd.read_sql_query(f"SELECT * FROM 'n_m'", cnx)
    print(nm_df)
    additional_stat_cols = [
        # 'orig_bandwidth',
        # 'cm_bandwidth',
        'sb_k',
        'par_sb_k',
        'sb_n_iters',
        'par_sb_n_iters',
        # 'pr_struct_size',
    ]
    df = df \
        .drop(columns=['n', 'm']) \
        .set_index('graph_name'). \
        join(nm_df.set_index('graph_name'), how='left')
    df = df[df['par_sb_n_iters'].notna()]
    df = df.loc[:, (df != 0).any(axis=0)]
    df = df.dropna(axis=1, how='all')

    df['par_wwr'] = df['par_sb_n_iters'] * df['par_sb_k'] / df['n']
    df['wwr'] = df['sb_n_iters'] * df['sb_k'] / df['n']
    df = df.sort_values(by='par_wwr', ascending=True)
    print(df.shape)
    print(f"{df=}")
    num_zeros = (df == 0).astype(int).sum(axis=0)
    num_nulls = (df.isnull()).astype(int).sum(axis=0)
    num_zeros = num_zeros.sort_values()
    for i, r in enumerate(num_zeros):
        print(i, r, num_zeros.index[i])
    df.to_csv('wing_width_ratio_comparison.csv')
    return


if __name__ == '__main__':
    config.init()

    main()
