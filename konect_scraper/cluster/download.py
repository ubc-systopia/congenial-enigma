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


def get_max_array_jobs():
    res = subprocess.check_output("scontrol show config | grep MaxArraySize | cut -d= -f2")
    print(res)
def download_and_extract(graph_name, konect_url, data_url):
    return 

def write_sbatch_array_csv(rows):
    # each row corresponds to `graph_name, konect_url, data_url`
    csvs_dir = config.settings['compute_canada']['job_array_dir']
    Path(csvs_dir).mkdir(parents=True, exist_ok=True)
    with open(os.path.join(csvs_dir, 'download.csv'), 'w') as f:
        for k in rows[0].keys()[:-1]:
            f.write(str(k) +',')
        f.write(rows[0].keys()[-1])
        f.write('\n')
        for r in rows:
            for k in rows[0].keys()[:-1]:
                f.write(str(r[k]) + ',')
            f.write(str(r[rows[0].keys()[-1]]))
            f.write('\n')
            


def main(graph_type, graph_ns):
    settings = config.settings
    rows = get_graphs_by_graph_numbers(graph_ns, graph_type)
    print(f"Downloading the following {graph_type} graphs to {settings['graphs_dir']}")
    df = rows_to_df(rows)
    print(df['graph_name'].values)
    rows = get_all_graphs_by_graph_names(df['graph_name'].values)
    write_sbatch_array_csv(rows)

    # run the sbatch command on 
    return 
