

import sqlite3
from konect_scraper import column_names, config
from konect_scraper.calc_stats import compute_stats
from konect_scraper.cluster.main import parse_and_init_data_dir
import pandas as pd
from konect_scraper.scrape_konect_stats import write_to_sqlite3

from konect_scraper.sql import append_df_to_table


def main(config_path, config_idx):
    global settings
    column_names.init()
    settings = config.settings
    df = pd.read_csv(config_path)
    row = df.iloc[int(config_idx)]
    graph_name = row['graph_name']
    mem = row['mem']
    settings['slurm_params'] = {}
    settings['slurm_params']['mem'] = mem
    feats_df = pd.DataFrame(columns=column_names.features_col_names.keys())

    d = compute_stats(graph_name)
    d['graph_name'] = graph_name
    feats_df = pd.concat([feats_df, pd.DataFrame([d])], ignore_index=True)
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    timeout = config.settings['sqlite3']['timeout']
    conn = sqlite3.connect(db_path, timeout=timeout)
    write_to_sqlite3(feats_df, 'features', conn)
    # append_df_to_table(feats_df, 'features')
    return 


if __name__ == '__main__':
    args = parse_and_init_data_dir()
    config.init(args.data_dir)
    main(args.config_path, args.config_idx)