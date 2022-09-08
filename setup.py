import logging
import os
import sqlite3
from datetime import datetime
import subprocess
from konect_scraper import config, column_names
import multiprocessing
import psutil


def main():
    config.init()
    settings = config.settings

    cmake_build_dir = settings['cmake_build_dir']
    cmake_build_type = settings['cmake_build_type']
    cmake_make_program = settings['cmake_make_program']  # ninja

    graph_preprocess_dir = settings["graph_preprocess_dir"]
    rabbit_home = settings["rabbit_home"]

    cmake_executable = "cmake"
    make_executable = "make"

    # build graph_preprocess
    args = [
        cmake_executable,
        f"-DCMAKE_BUILD_TYPE={cmake_build_type}",
        f"-DCMAKE_MAKE_PROGRAM={cmake_make_program}",
        "-G", "Ninja",
        "-S", graph_preprocess_dir,
        "-B", cmake_build_dir,
    ]
    print(" ".join(args))
    res = subprocess.check_output(args)
    print(res.decode('ascii'))

    # build rabbit submodule
    args = [
        make_executable,
    ]
    print(" ".join(args))
    subprocess.check_output(args, cwd=os.path.join(rabbit_home, 'demo'))
    print(res.decode('ascii'))

    n_threads = psutil.cpu_count()

    # install graph_preprocess
    args = [
        cmake_executable,
        "--build", cmake_build_dir,
        "--target", "graph_preprocess",
        "-j", str(n_threads)
    ]
    print(" ".join(args))
    res = subprocess.check_output(args)
    print(res.decode('ascii'))

    # install slashburn
    args = [
        cmake_executable,
        "--build", cmake_build_dir,
        "--target", "slashburn",
        "-j", str(n_threads)
    ]
    print(" ".join(args))
    res = subprocess.check_output(args)
    print(res.decode('ascii'))


if __name__ == "__main__":
    main()
