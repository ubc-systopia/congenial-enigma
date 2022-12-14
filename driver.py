import argparse
import json
import subprocess
import yaml
from konect_scraper.utilities import batch_submit

def construct_main_args_from_config(args):
    config_path = args.config_path
    slurm = args.slurm
    data = {}
    with open(config_path, "r") as stream:
        try:
            data = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    args = []
    
    if slurm:
        args = ['python', '-m', 'konect_scraper.cluster.main']
    else:
        args = ['python', 'main.py']

    val_args = ['mode', 'data-dir']
    if slurm:
        val_args += ['time', 'mem', 'cpus-per-task', 'constraint']

    bool_args = ['directed', 'overwrite', 'debug']
    for val_arg in val_args:
        args += [f'--{val_arg}', str(data[val_arg])]

    if data['reorder'] == 'all':
        args += [f'--reorder', 'all']
    else:
        args += [f'--reorder'] 
        for order in data['reorder']:
            args += [order]

    for bool_arg in bool_args:
        if data[bool_arg]:
            args += [f'--{bool_arg}']

    mn_graph_n = data['graph-numbers']['min']
    mx_graph_n = data['graph-numbers']['max']
    args += [f'--graph-numbers', str(mn_graph_n), str(mx_graph_n)]

    return args


def main(args):

    main_args = construct_main_args_from_config(args)
    print(f' '.join(main_args))

    # batch_submit.main(main_args)
    res = subprocess.check_output(main_args)
    print(res.decode('utf-8'))

    return


if __name__ == '__main__':
    argparse_desc = """
    A wrapper that parses konect_scraper.cluster configuration, and calls 
    konect_scraper.cluster.main 
    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-c', '--config-path', required=True,
                        help='absolute path to config yaml file required for '
                        'slurm execution of konect_scraper')

    parser.add_argument('-s', '--slurm', action='store_true',
                        help='If set, submits an array job on slurm using '
                        'konect_scraper.cluster.main. Otherwise, runs locally using main.py.')

    args = parser.parse_args()
    # print(args)
    main(args)
