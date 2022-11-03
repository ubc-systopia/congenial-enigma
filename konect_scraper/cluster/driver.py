import argparse
import json
import subprocess
import yaml


def construct_main_args_from_config(config_path):
    data = {}
    with open(config_path, "r") as stream:
        try:
            data = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    args = ['python', '-m', 'konect_scraper.cluster.main']
    val_args = ['mode', 'reorder', 'data-dir', 'time', 'mem', 'cpus-per-task']
    bool_args = ['directed', 'overwrite']
    for val_arg in val_args:
        args += [f'--{val_arg}', str(data[val_arg])]
    for bool_arg in bool_args:
        if data[bool_arg]:
            args += [f'--{bool_arg}']

    mn_graph_n = data['graph-numbers']['min']
    mx_graph_n = data['graph-numbers']['max']
    args += [f'--graph-numbers', str(mn_graph_n), str(mx_graph_n)]

    return args


def main(args):

    main_args = construct_main_args_from_config(args.config_path)
    print(f' '.join(main_args))
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

    args = parser.parse_args()
    main(args)
