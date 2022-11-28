import os
from konect_scraper import config
from konect_scraper.cluster.execute import prep_sbatch_array_submit
from konect_scraper.cluster.util import parse_sacct
from konect_scraper.util import remove_all_files_in_directory
from slurmpy import Slurm


SBATCH_SCRIPT="""
module load singularity/3.7

CONFIG_ID=${SLURM_ARRAY_TASK_ID}
offset_line_number=$((CONFIG_ID + 2))
line=$(sed "${offset_line_number}q;d" ${LOCAL_CFG})
echo ${line}
SCRIPTS_DIR=/congenial-enigma/konect_scraper/cluster/scripts/

echo "singularity exec --bind ${DATA_DIR}:/data,${REPO_HOME}:/congenial-enigma \
        ${IMAGE} \
        ${SCRIPTS_DIR}singularity-exec.sh ${CFG_FILE} ${CONFIG_ID} ${MODE}"

singularity exec --bind ${DATA_DIR}:/data,${REPO_HOME}:/congenial-enigma \
        ${IMAGE} \
        ${SCRIPTS_DIR}singularity-exec.sh ${CFG_FILE} ${CONFIG_ID} ${MODE}
"""

def main(graph_type, graph_ns, slurm_params, mode_str,
         vertex_orders=None, overwrite=False, depends_on=None):
    # ignore rev_cm, it is computed when computing cm
    if vertex_orders:
        skip_vorders = ['rev_cm', 'orig']
        vertex_orders = [o for o in vertex_orders if o not in skip_vorders]
    prep_sbatch_array_submit(graph_type, graph_ns, slurm_params, mode_str,
         vertex_orders, overwrite)
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
    
    singularity_image_path = settings['compute_canada']['image']
    remove_all_files_in_directory(log_dir)
    min_graph_n = graph_ns[0]
    max_graph_n = graph_ns[-1]
    n_graphs = max_graph_n - min_graph_n
    n_array_jobs = n_graphs 
    edge_orders = list(settings['edge_orderings'].keys())
    if vertex_orders:
        n_array_jobs += 1
        n_array_jobs *= len(vertex_orders)
        if mode_str == 'pr_expt':
            n_array_jobs *= len(edge_orders)
    slurm_init_params = {
        'account': settings['compute_canada']['account'],
        # 'user': 'atrostan',
        'time': slurm_params['time'],
        'mem': slurm_params['mem'],
        'cpus-per-task': slurm_params['cpus-per-task'],
        'constraint': slurm_params['constraint'],
        'nodes': '1-1',
        'array': f'0-{n_array_jobs}',
    }
    cmd_kwargs = {
        'CFG_FILE': mounted_config,
        'IMAGE': singularity_image_path,
        'REPO_HOME': str(config.settings['compute_canada']['repo_root']),
        'DATA_DIR': config.settings['compute_canada']['data_dir'],
        'MODE': mode_str,
        'LOCAL_CFG': local_config,
        'SCRIPTS_DIR': scripts_dir,
    }
    print(cmd_kwargs)

    s = Slurm(mode_str, slurm_init_params)
    
    # parse_sacct(['50743422_0'])
    # A task of this job array can begin execution after the corresponding task 
    # ID in the specified job has completed successfully 
    depends_how='aftercorr' 
    if depends_on:
        job_id = s.run(command=SBATCH_SCRIPT, cmd_kwargs=cmd_kwargs,
        depends_on=depends_on, depends_how=depends_how)
    else:
        job_id = s.run(command=SBATCH_SCRIPT, cmd_kwargs=cmd_kwargs)

    # # submit a job that will execute once SBATCH_SCRIPT completes
    # # to record all sacct -j metrics 

    # array_job_ids = [f'{job_id}_{i}' for i in range(n_graphs + 1)]
    # print(f"Submitted: {config_filename}")
    # print(f"Job ID: {job_id}")
    return job_id

   



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Use slurmpy to spawn '
        'compute canada jobs with dependencies')
    parser.add_argument('--local-config', required=True,
                        help='path to local configuration file')
    
    parser.add_argument('--mounted-config', required=True,
                        help='path to mounted configuration file')
    
    parser.add_argument('--singularity-image', required=True,
                        help='absolute path to singularity image')
    
    parser.add_argument('--repo-root', required=True,
                        help="absolute path to this repository's root")
    
    parser.add_argument('--data-dir', required=True,
                        help='absolute path to graphs.db and graphs')
    
    parser.add_argument('--mode-str', required=True,
                        help='name of execution mode')
    
    parser.add_argument('--time', required=True,
                        help='sbatch time param')
    parser.add_argument('--mem', required=True,
                        help='sbatch mem param')
    parser.add_argument('--cpus-per-task', required=True,
                        help='sbatch cpus-per-task param')
    
    parser.add_argument('--constraint', required=True,
                        help='sbatch constraint param')
    
    args = parser.parse_args()
    print(args)
    # return 
    # main(args)


"""
delete all 'core.*' files in directory and subdirectories
find . -regex "^.*/core.[0-9]*$" -type f -delete
"""