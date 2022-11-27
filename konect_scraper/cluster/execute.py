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


def write_sbatch_array_csv(rows, name, vertex_orders=None, edge_orders=None, overwrite=False):
    # each row corresponds to `graph_name, konect_url, data_url`
    csvs_dir = config.settings['compute_canada']['job_array_dir']
    Path(csvs_dir).mkdir(parents=True, exist_ok=True)
    print(f'{overwrite=}')
    graph_column_names = ['graph_name', 'konect_url', 'data_url']
    column_names = graph_column_names.copy()
    lines = []
    if vertex_orders:
        column_names += ['vertex_order']
    if edge_orders:
        column_names += ['edge_order']
    column_names += ['overwrite']

    lines.append(",".join(column_names))  # header

    for row in rows:
        line = ",".join([str(row[k]) for k in graph_column_names])
        if vertex_orders:
            for vertex_order in vertex_orders:
                if edge_orders:
                    for edge_order in edge_orders:
                        if overwrite:
                            lines.append(
                                line + f',{vertex_order}' + f',{edge_order}' + f',1')
                        else:
                            lines.append(
                                line + f',{vertex_order}' + f',{edge_order}' + f',0')
                else:
                    if overwrite:
                        lines.append(line + f',{vertex_order}' + f',1')
                    else:
                        lines.append(line + f',{vertex_order}' + f',0')
            
        else:
            if overwrite:
                lines.append(line + f',1')
            else:
                lines.append(line + f',0')
    
    with open(os.path.join(csvs_dir, f'{name}.csv'), 'w') as f:
        for l in lines:
            f.write(f'{l}\n')

def prep_sbatch_array_submit(graph_type, graph_ns, slurm_params, mode_str,
         vertex_orders=None, overwrite=False):

    settings = config.settings
    rows = get_graphs_by_graph_numbers(graph_ns, graph_type)
    df = rows_to_df(rows)
    rows = get_all_graphs_by_graph_names(df['graph_name'].values)
    edge_orders = list(settings['edge_orderings'].keys())
    config_filename = f'{mode_str}_{graph_ns[0]}_{graph_ns[-1]}'

    if mode_str == 'reorder':
        write_sbatch_array_csv(
            rows, config_filename, vertex_orders=vertex_orders, overwrite=overwrite)
    elif mode_str == 'pr_expt':
        write_sbatch_array_csv(
            rows, config_filename, vertex_orders=vertex_orders, 
            edge_orders=edge_orders, overwrite=overwrite)
    else:
        write_sbatch_array_csv(rows, config_filename, overwrite=overwrite)

    return 


def main(graph_type, graph_ns, slurm_params, mode_str,
         vertex_orders=None, overwrite=False):
    
    prep_sbatch_array_submit(graph_type, graph_ns, slurm_params, mode_str,
         vertex_orders=None, overwrite=False)
    settings = config.settings
    config_filename = f'{mode_str}_{graph_ns[0]}_{graph_ns[-1]}'
    scripts_dir = config.settings['compute_canada']['scripts_dir']
    csvs_dir = config.settings['compute_canada']['job_array_dir']

    local_config = os.path.join(csvs_dir, f'{config_filename}.csv')
    # CONFIG_FILE will be locally mounted to /
    mounted_config = os.path.join(
        '/', settings['repo_name'], 'konect_scraper',
        'cluster', 'csvs', f'{config_filename}.csv',
    )

    # clean slurm logs from previous execution
    log_dir = os.path.join(
        config.settings['logging']['slurm_log_dir'],
        mode_str
    )

    remove_all_files_in_directory(log_dir)

    args = [
        # os.path.join(scripts_dir, 'singularity-start.sh'),
        'python -m '
        '--local-config', local_config,
        '--mounted-config', mounted_config,
        '--singularity-image', config.settings['compute_canada']['image'],
        '--repo-root', str(config.settings['compute_canada']['repo_root']),
        '--data-dir', config.settings['compute_canada']['data_dir'],
        '--mode-str', mode_str,
        '--time', slurm_params['time'],
        '--mem', slurm_params['mem'],
        '--cpus-per-task', slurm_params['cpus-per-task'],
        '--constraint', slurm_params['constraint']
    ]

    print("CALLING: " + " ".join(args))
    # run the sbatch command
    res = subprocess.check_output(args)
    print(res.decode('utf-8'))
    return
