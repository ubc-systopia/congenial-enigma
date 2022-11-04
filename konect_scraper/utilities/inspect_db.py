import argparse
import os
import sqlite3
import pandas as pd


def main(args):
    data_dir = args.data_dir
    print(f'{data_dir}')
    db_path = os.path.join(data_dir, 'graphs.db')
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT * FROM n_m", conn)
    print(df.to_string())

if __name__ == '__main__':

    argparse_desc = """
    Module to interface and view contents of graphs.db
    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-a', '--data-dir', required=True, 
                        help="Absolute path to the data directory containing graphs.db and all graph files.")

    args = parser.parse_args()

    main(args)