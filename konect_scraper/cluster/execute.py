import logging

from konect_scraper.util import remove_all_files_in_directory
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


def write_sbatch_array_csv(rows, name, vertex_orders=None, overwrite=False):
    # each row corresponds to `graph_name, konect_url, data_url`
    csvs_dir = config.settings['compute_canada']['job_array_dir']
    Path(csvs_dir).mkdir(parents=True, exist_ok=True)

    graph_column_names = ['graph_name', 'konect_url', 'data_url']
    column_names = graph_column_names.copy()
    lines = []
    if vertex_orders:
        column_names += ['vertex_order']
    if overwrite:
        column_names += ['overwrite']

    lines.append(",".join(column_names))  # header

    for row in rows:
        line = ",".join([str(row[k]) for k in graph_column_names])
        if vertex_orders:
            for vertex_order in vertex_orders:
                line += f',{vertex_order}'
                if overwrite:
                    line += f',1'
                else:
                    line += f',0'
                lines.append(line)

        else:
            lines.append(line)
    with open(os.path.join(csvs_dir, f'{name}.csv'), 'w') as f:
        for l in lines:
            f.write(f'{l}\n')


def main(graph_type, graph_ns, slurm_params, mode_str,
         vertex_orders=None, overwrite=False):
    settings = config.settings
    rows = get_graphs_by_graph_numbers(graph_ns, graph_type)
    df = rows_to_df(rows)
    rows = get_all_graphs_by_graph_names(df['graph_name'].values)

    if mode_str == 'reorder' or mode_str == 'pr_expt':
        write_sbatch_array_csv(
            rows, mode_str, vertex_orders, overwrite=overwrite)
    else:
        write_sbatch_array_csv(rows, mode_str, overwrite=overwrite)

    scripts_dir = config.settings['compute_canada']['scripts_dir']
    csvs_dir = config.settings['compute_canada']['job_array_dir']

    local_config = os.path.join(csvs_dir, f'{mode_str}.csv')
    # CONFIG_FILE will be locally mounted to /
    mounted_config = os.path.join(
        '/', settings['repo_name'], 'konect_scraper',
        'cluster', 'csvs', f'{mode_str}.csv',
    )

    # clean slurm logs from previous execution
    log_dir = os.path.join(
        config.settings['logging']['slurm_log_dir'],
        mode_str
    )
    remove_all_files_in_directory(log_dir)

    args = [
        os.path.join(scripts_dir, 'singularity-start.sh'),
        local_config,
        mounted_config,
        log_dir,
        config.settings['compute_canada']['image'],
        str(config.settings['compute_canada']['repo_root']),
        config.settings['compute_canada']['data_dir'],
        mode_str,
        slurm_params['time'],
        slurm_params['mem'],
        slurm_params['cpus-per-task'],
    ]
    print("CALLING: " + " ".join(args))
    # run the sbatch command
    res = subprocess.check_output(args)
    print(res.decode('utf-8'))
    return
