import os.path
from os import listdir
from os.path import isfile, join
import json
from pathlib import Path
import urllib.request
from urllib.parse import urlparse
import tarfile
import shutil
import os
import sqlite3
import subprocess
import gzip

import numpy as np

from konect_scraper import config
from konect_scraper.config import IOMode
from konect_scraper.util import single_val_numeric_set, get_size, get_volume, get_directed, set_n, set_m, \
    single_val_get, get_size_in_memory, set_n_m, save_ground_truth_pr
import logging

def get_edge_list_filename(directory):
    fs = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            fs.append(os.path.join(path, name))

    sorted_files = sorted(fs, key=os.path.getsize)
    for f in sorted_files:
        if os.path.basename(f).split('.')[0] == 'out':
            return f

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
    with open(path, 'r') as source_file:
        while i < n_comment_lines - 1:
            source_file.readline()
            i += 1

        tmp_path = path + '.tmp'
        # this will truncate the file, so need to use a different file name:
        with open(tmp_path, 'w') as target_file:
            shutil.copyfileobj(source_file, target_file)
    # delete the original file (contains comments) and copy the file (without comments) in its place
    os.remove(path)
    shutil.copyfile(tmp_path, path)
    os.remove(tmp_path)

    return


def get_n_columns(path):
    n_comment_lines = get_n_comments_at_top_of_file(path)
    i = 0
    with open(path, 'r') as source_file:
        while i < n_comment_lines + 1:
            line = source_file.readline()
            i += 1
    if ' ' in line:
        return len(line.split())
    elif '\t' in line:
        return len(line.split('\t'))


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

    logging.info(f"Downloading edge list from {url} to {tmp_file}..")

    urllib.request.urlretrieve(url, tmp_file)
    suffixes = Path(tmp_file).suffixes

    mode = ''
    tmp_dir = os.path.join(directory, 'tmp')

    if '.tar' in suffixes:
        if '.bz2' in suffixes:
            mode = 'r:bz2'
        if '.gz' in suffixes:
            mode = 'r:gz'
        logging.info(f"Extracting {tmp_file} - mode: {mode}..")

        tar = tarfile.open(tmp_file, mode)
        tar.extractall(tmp_dir)
        tar.close()
        edge_list_filename = get_edge_list_filename(tmp_dir)
        shutil.copyfile(edge_list_filename, graph_path)

    if '.txt' in suffixes and '.gz' in suffixes:
        with gzip.open(tmp_file, 'rb') as f_in:
            tmp_txt = os.path.splitext(tmp_file)[0]
            with open(tmp_txt, 'wb') as f_out:
                logging.info(
                    f"Extracting {tmp_file} into {tmp_txt} mode: gzip..")
                shutil.copyfileobj(f_in, f_out)
        edge_list_filename = tmp_txt
        shutil.copyfile(edge_list_filename, graph_path)
        os.remove(tmp_txt)

    # copy the edgelist from the extracted directory (this should be the largest file in the directory)
    # remove all first lines that contain a comment
    remove_comments_from_edgelist(graph_path)
    n_columns = get_n_columns(graph_path)
    if n_columns > 2:
        # todo: use np to remove unneeded cols - hacky
        comments = settings['comment_strings']
        arr = np.loadtxt(graph_path, comments=comments).astype(np.uint32)
        np.savetxt(graph_path, arr[:, :2], fmt='%i', )
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


def run_graph_preprocess(input_path, output_path, directed, n, m, io_modes, graph_name):
    settings = config.settings
    executable = settings['graph_preprocess_executable']
    sqlite3_db_path = settings['sqlite3']['sqlite3_db_path']
    args = [executable]
    # TODO directedness not required for ingest?
    # if directed:
    #     args += ['-d']

    for io_mode in io_modes:
        match io_mode:
            case IOMode.text:
                args += ['-t']
            case IOMode.binary:
                args += ['-i']

    args += [
        '-g', input_path,
        '-b', sqlite3_db_path,
        '-m', str(m),
        '-o', output_path
    ]

    logging.info(f"Executing: " + ' '.join(args))

    res = subprocess.check_output(args)
    # parse the output of the simplify program to get the number of vertices and edges
    # in the compressed, simplified representation
    lines = res.decode('ascii').split('\n')
    n = int(lines[0].replace('n: ', ''))
    m = int(lines[1].replace('m: ', ''))

    # convert the graph to a format dbg requires
    dbg_clean_el_executable = settings['dbg_clean_el_executable']
    dbg_convert_script = settings['dbg_convert_script']
    dbg_datasets_dir = settings['dbg_datasets_dir']
    dbg_home = settings['dbg_home']

    edgelist_file_suffix = settings['edgelist_file_suffix']
    dbg_edgelist_file_suffix = settings['dbg']['edgelist_file_suffix']

    env = {
        **os.environ,
        "DBG_ROOT": dbg_home,
    }
    args = [
        dbg_clean_el_executable,
        f"{output_path}.{edgelist_file_suffix}",
        f"{output_path}.{dbg_edgelist_file_suffix}"  # save to same dir
    ]
    logging.info(f"Executing: " + ' '.join(args))
    res = subprocess.check_output(args, env=env)

    args = [
        dbg_convert_script,
        # os.path.join(dbg_datasets_dir, graph_name),
        output_path
    ]
    logging.info(f"Executing: " + ' '.join(args))
    res = subprocess.check_output(args, env=env)

    return n, m


def compress(directory, directed, n, m, io_modes, graph_name):
    """
    use igraph_simplify to compress a graph's edgelist in directory
    :param directory:
    :return:
    """
    settings = config.settings
    graph_path = os.path.join(directory, settings["orig_el_file_name"])
    compressed_graph_path = os.path.join(
        directory, settings["compressed_el_file_name"])
    comp_n, comp_m = run_graph_preprocess(
        graph_path, compressed_graph_path, directed, n, m, io_modes, graph_name)

    # TODO record the compressed raw text file in db
    return comp_n, comp_m


def main(rows, io_modes):
    """
    1. Download graph from urls specified in datasets.json
    2. Extract graph file into graph's directory
    :return:
    """

    settings = config.settings

    datasets_json_path = settings['datasets_json_path']
    graphs_dir = settings['graphs_dir']

    for row in rows:
        graph_name = row['graph_name']
        konect_url = row['konect_url']
        data_url = row['data_url']
        if data_url == "none":
            continue
        comp_size = single_val_get(
            'compressed_txt_file_size', 'metadata', graph_name)
        # if comp_size > 0:  # dataset has been downloaded already
        #     continue

        graph_dir = os.path.join(graphs_dir, graph_name)
        compressed_edge_list_path = os.path.join(
            graph_dir, f"{settings['compressed_el_file_name']}.net")

        # create directory for the graph (if not exists)
        Path(graph_dir).mkdir(parents=True, exist_ok=True)
        download_graph(data_url, graph_dir)
        # update the graph's metadata to point to its directory
        update_meta_dir(graph_dir, graph_name)

        # update the text file size in the database
        graph_path = os.path.join(graph_dir, settings["orig_el_file_name"])
        single_val_numeric_set('txt_file_size', 'metadata',
                               graph_name, os.path.getsize(graph_path))

        
    return
