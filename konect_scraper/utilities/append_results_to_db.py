


import argparse
import os, shutil
import pandas as pd
from konect_scraper import config
from konect_scraper.sql import append_df_to_table


def delete_contents_of_folder(folder):
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

def main(args):
    for path, subdirs, files in os.walk(args.results_dir):
        for name in files:
            results_path = os.path.join(path, name)
            append_df_to_table(pd.read_csv(results_path), 'pr_expts')

    # delete after append
    delete_contents_of_folder(args.results_dir)

    return 

if __name__ == '__main__':
    argparse_desc = """
        A utility module that traverses a directory and appends experimental
        results to a sqlite3 database
    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-r', '--results-dir', required=True,
                        help="""Absolute path to directory containing results
                        of experiments. Contains subdirectories for each graph
                        experiments were run on. 
                        """)

    parser.add_argument('-d', '--data-dir', required=True,
                        help="""Location of graphs.db (sqlite3 table to append
                        results to)
                        """)

    args = parser.parse_args()
    config.init(args.data_dir)
    # print(args)
    main(args)
