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
    rabbit_cmake_build_dir = settings['rabbit_cmake_build_dir']
    cmake_build_type = settings['cmake_build_type']
    cmake_make_program = settings['cmake_make_program']  # ninja

    graph_preprocess_dir = settings["graph_preprocess_dir"]
    rabbit_home = settings["rabbit_home"]

    cmake_executable = settings["cmake_executable"]
    make_executable = settings["make_executable"]

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

# /home/atrostan/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/222.3739.54/bin/cmake/linux/bin/cmake --build /home/atrostan/Workspace/repos/congenial-enigma/rabbit_order/demo/cmake-build-debug --target all -j 9


    # build rabbit submodule
    args = [
        cmake_executable,
        f"-DCMAKE_BUILD_TYPE={cmake_build_type}",
        f"-DCMAKE_MAKE_PROGRAM={cmake_make_program}",
        "-G", "Ninja",
        "-S", os.path.join(rabbit_home, "demo"),
        "-B", rabbit_cmake_build_dir,
    ]
    print(" ".join(args))
    subprocess.check_output(args, cwd=os.path.join(rabbit_home, 'demo'))
    print(res.decode('ascii'))

    n_threads = config.settings['n_threads']

    # install rabbit submodule
    args = [
        cmake_executable,
        "--build", rabbit_cmake_build_dir,
        "--target", "all",
        "-j", str(n_threads)
    ]

    print(" ".join(args))
    res = subprocess.check_output(args)
    print(res.decode('utf-8'))

    # install graph_preprocess
    args = [
        cmake_executable,
        "--build", cmake_build_dir,
        "--target", "graph_preprocess",
        "-j", str(n_threads)
    ]
    print(" ".join(args))
    res = subprocess.check_output(args)
    print(res.decode('utf-8'))

    # install slashburn
    args = [
        cmake_executable,
        "--build", cmake_build_dir,
        "--target", "slashburn",
        "-j", str(n_threads)
    ]
    print(" ".join(args))
    res = subprocess.check_output(args)
    print(res.decode('utf-8'))

    dbg_home = settings['dbg_home']
    dbg_apps_dir = settings['dbg_apps_dir']
    # build dbg
    print(dbg_apps_dir)
    args = [
        make_executable, "clean",
    ]
    print(" ".join(args))
    subprocess.check_output(args, cwd=dbg_apps_dir)
    print(res.decode('ascii'))

    # args = [
    #     make_executable, "cleansrc",
    # ]
    # print(" ".join(args))
    # subprocess.check_output(args, cwd=dbg_apps_dir)
    # print(res.decode('ascii'))

    args = [
        make_executable,
        "-j", str(n_threads),
    ]
    print(" ".join(args))
    subprocess.check_output(args, cwd=dbg_apps_dir)
    print(res.decode('ascii'))

if __name__ == "__main__":
    main()
