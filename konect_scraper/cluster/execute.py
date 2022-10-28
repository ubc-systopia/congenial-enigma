import logging
from ..sql import *
from .. import config
from .util import *
import os
from pathlib import Path
import subprocess
"""
A wrapper that given graph indices and type of graph type, will submit 
sbatch array jobs to download the relevant in parallel

Each sbatch will download a single graph
So, if given a 100 graphs to download, assume 100 sbatch jobs will be submitted
"""

def write_sbatch_array_csv(rows, name):
    # each row corresponds to `graph_name, konect_url, data_url`
    csvs_dir = config.settings['compute_canada']['job_array_dir']
    Path(csvs_dir).mkdir(parents=True, exist_ok=True)
    with open(os.path.join(csvs_dir, f'{name}.csv'), 'w') as f:
        for k in rows[0].keys()[:-1]:
            f.write(str(k) +',')
        f.write(rows[0].keys()[-1])
        f.write('\n')
        for r in rows:
            for k in rows[0].keys()[:-1]:
                f.write(str(r[k]) + ',')
            f.write(str(r[rows[0].keys()[-1]]))
            f.write('\n')
            


def main(graph_type, graph_ns, mode_str):
    settings = config.settings
    rows = get_graphs_by_graph_numbers(graph_ns, graph_type)
    df = rows_to_df(rows)
    print(df['graph_name'].values)
    rows = get_all_graphs_by_graph_names(df['graph_name'].values)
    write_sbatch_array_csv(rows, mode_str)
    scripts_dir = config.settings['compute_canada']['scripts_dir']
    csvs_dir = config.settings['compute_canada']['job_array_dir']
    print(f"{config.settings['compute_canada']['repo_root']}")
    args=[
        os.path.join(scripts_dir, mode_str, 'singularity-start.sh'),
        os.path.join(csvs_dir, f'{mode_str}.csv'), # CONFIG_FILE
        config.settings['logging']['slurm_log_dir'], # LOGDIR
        config.settings['compute_canada']['image'], # IMAGE
        str(config.settings['compute_canada']['repo_root']), # REPO_HOME
        config.settings['compute_canada']['data_dir'], # DATA_DIR
    ]
    print(" ".join(args))
    # run the sbatch command 
    res = subprocess.check_output(args)
    print(res.decode('utf-8'))
    return 
