import os.path
from os import listdir
from os.path import isfile, join
import config
import json
from pathlib import Path
import urllib.request
from urllib.parse import urlparse
import tarfile
import shutil
import os
import sqlite3
from util import *
import subprocess
import gzip

def get_edge_list_filename(directory):
    fs = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            fs.append(os.path.join(path, name))

    sorted_files = sorted(fs, key=os.path.getsize)

    return sorted_files[-1]


def is_comment(s):
    b = "%" in s or "#" in s
    return b


def get_n_comments_at_top_of_file(path):
    # count how many comment lines are in the top of the file
    source_file = open(path, 'r')
    i = 0
    while True:
        line = source_file.readline().strip()
        i += 1
        if not is_comment(line):
            break

    return i


def remove_comments_from_edgelist(path):
    """
    The first lines of edge list files typically contain comments
    Delete these commented lines and overwrite the original edge list to strictly contain the edges of the graph
    :param path:
    :return:
    """
    n_comment_lines = get_n_comments_at_top_of_file(path)
    i = 0
    source_file = open(path, 'r')

    while i < n_comment_lines - 1:
        source_file.readline()
        i += 1

    tmp_path = path + '.tmp'
    # this will truncate the file, so need to use a different file name:
    target_file = open(tmp_path, 'w')

    shutil.copyfileobj(source_file, target_file)

    # delete the original file (contains comments) and copy the file (without comments) in its place
    os.remove(path)
    shutil.copyfile(tmp_path, path)
    os.remove(tmp_path)

    return


def download_graph(url, directory):
    """
    Download a graph, extract it, store it in dir, and update the db 'directory' column to point to the directory
    :param url:
    :param directory:
    :return:
    """
    settings = config.settings
    graph_path = os.path.join(directory, settings["orig_el_file_name"])
    fname = os.path.basename(urlparse(url).path)
    tmp_file = os.path.join(directory, fname)

    print(f"Downloading edge list from {url} to {tmp_file}..")

    urllib.request.urlretrieve(url, tmp_file)
    suffixes = Path(tmp_file).suffixes

    mode = ''
    tmp_dir = os.path.join(directory, 'tmp')

    if '.tar' in suffixes:
        if '.bz2' in suffixes:
            mode = 'r:bz2'
        if '.gz' in suffixes:
            mode = 'r:gz'
        print(f"Extracting {tmp_file} - mode: {mode}..")

        tar = tarfile.open(tmp_file, mode)
        tar.extractall(tmp_dir)
        tar.close()
        edge_list_filename = get_edge_list_filename(tmp_dir)
        shutil.copyfile(edge_list_filename, graph_path)

    if '.txt' in suffixes and '.gz' in suffixes:
        with gzip.open(tmp_file, 'rb') as f_in:
            tmp_txt = os.path.splitext(tmp_file)[0]
            with open(tmp_txt, 'wb') as f_out:
                print(f"Extracting {tmp_file} into {tmp_txt} mode: gzip..")
                shutil.copyfileobj(f_in, f_out)
        edge_list_filename = tmp_txt
        shutil.copyfile(edge_list_filename, graph_path)
        os.remove(tmp_txt)

    # copy the edgelist from the extracted directory (this should be the largest file in the directory)

    # remove all first lines that contain a comment
    remove_comments_from_edgelist(graph_path)

    # delete the archive and extracted file (only care about the edgelist)
    os.remove(tmp_file)
    shutil.rmtree(tmp_dir, ignore_errors=True)

    # update the graph's edgelist size (raw text file) in db
    return


def update_meta_dir(directory, name):
    """
    Find the graph in the metadata table and update its directory path
    :param directory:
    :param name:
    :return:
    """
    db_path = config.settings['sqlite3']['sqlite3_db_path']
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    query = """update metadata set directory = ? where graph_name = ?"""
    cursor.execute(query, (directory, name))
    conn.commit()
    conn.close()
    return


def run_graph_simplify(input_path, output_path, directed, n, m):
    settings = config.settings
    executable = settings['graph_simplify_executable']
    sqlite3_db_path = settings['sqlite3']['sqlite3_db_path']
    args = [executable]
    if directed:
        args += ['-d']
    args += [
        '-g', input_path,
        '-b', sqlite3_db_path,
        '-n', str(n),
        '-m', str(m),
        '-o', output_path
    ]

    print(' '.join(args))

    res = subprocess.check_output(args)
    # parse the output of the simplify program to get the number of vertices and edges
    # in the compressed, simplified representation

    lines = res.decode('ascii').split('\n')
    n = int(lines[0].replace('n: ', ''))
    m = int(lines[1].replace('m: ', ''))

    return n, m


def compress(directory, directed, n, m):
    """
    use igraph_simplify to compress a graph's edgelist in directory
    :param directory:
    :return:
    """
    settings = config.settings
    graph_path = os.path.join(directory, settings["orig_el_file_name"])
    compressed_graph_path = os.path.join(directory, settings["compressed_el_file_name"])
    comp_n, comp_m = run_graph_simplify(graph_path, compressed_graph_path, directed, n, m)

    # TODO record the compressed raw text file in db

    return comp_n, comp_m


def main(datasets):
    """
    1. Download graph from urls specified in datasets.json
    2. Extract graph file into graph's directory
    3. Simplify edge list
    :return:
    """

    settings = config.settings

    datasets_json_path = settings['datasets_json_path']
    graphs_dir = settings['graphs_dir']

    for dataset in datasets:
        graph_name = dataset['name']
        konect_url = dataset['konect-url']
        data_url = dataset['data-url']

        graph_dir = os.path.join(graphs_dir, graph_name)
        # create directory for the graph (if not exists)
        Path(graph_dir).mkdir(parents=True, exist_ok=True)
        download_graph(data_url, graph_dir)
        # update the graph's metadata to point to its directory
        update_meta_dir(graph_dir, graph_name)

        # update the text file size in the database
        graph_path = os.path.join(graph_dir, settings["orig_el_file_name"])
        single_val_numeric_set('txt_file_size', 'metadata', graph_name, os.path.getsize(graph_path))

        # compress the graph's edgelist
        directed = bool(is_directed(graph_name))
        sz = get_size(graph_name)
        vol = get_volume(graph_name)
        print(f"{sz=}")
        print(f"{vol=}")
        print(f"{directed=}")
        n, m = compress(graph_dir, directed, sz, vol)
        set_n(graph_name, n)
        set_m(graph_name, m)

        # update the text file size in the database
        compressed_graph_path = os.path.join(graph_dir, settings["compressed_el_file_name"])
        single_val_numeric_set('compressed_txt_file_size', 'metadata', graph_name,
                               os.path.getsize(compressed_graph_path))
    return


if __name__ == '__main__':
    main()
