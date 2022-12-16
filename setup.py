import argparse
import logging
import os
import shutil
import sqlite3
from datetime import datetime
import subprocess
from pathlib import Path

from konect_scraper import config, column_names
import multiprocessing
import psutil


def clean_built_submodules(settings):
    submodule_paths = [
        # ParallelBatchRCM,
        settings['graph_preprocess_dir'],
        os.path.join(settings["rabbit_home"], 'demo'),
        settings['par_slashburn_dir'],
    ]
    for submodule_path in submodule_paths:
        print(
            f"rm -rf {os.path.join(submodule_path, settings['cmake_build_dir'])}")
        shutil.rmtree(
            os.path.join(submodule_path, settings['cmake_build_dir']),
            ignore_errors=True
        )

    print(f"cd {os.path.join(settings['dbg_home'], 'apps')}; make clean")
    subprocess.check_output(
        ['make', 'clean'],
        cwd=os.path.join(settings['dbg_home'], 'apps'),
    )

    args = [
        'rm', '-rf',
        f"{settings['dbg_convert_dir']}/convert",
        f"{settings['dbg_convert_dir']}/*.o",
    ]

    print(' '.join(args))
    subprocess.check_output(args)

    # build peregrine
    prgrn_dir = settings['peregrine']['prgrn_dir']
    args = ['make', 'clean',]
    print(f"cd {prgrn_dir}; make clean")
    subprocess.check_output(args, cwd=prgrn_dir)

    return


def main(args):
    config.init()
    settings = config.settings

    if args.clean:
        clean_built_submodules(settings)
        return

    cmake_build_dir = settings['cmake_build_dir']
    rabbit_cmake_build_dir = settings['rabbit_cmake_build_dir']
    cmake_build_type = settings['cmake_build_type']
    cmake_make_program = settings['cmake_make_program']  # ninja

    graph_preprocess_dir = settings["graph_preprocess_dir"]
    rabbit_home = settings["rabbit_home"]
    dbg_home = settings['dbg_home']
    par_slashburn_dir = settings['par_slashburn_dir']

    cmake_executable = settings["cmake_executable"]
    make_executable = settings["make_executable"]

    # # build and install ParallelBatchRCM
    # pbrcm_home = settings['pbrcm_home']
    # pbrcm_build_dir = os.path.join(pbrcm_home, "build")
    # Path(pbrcm_build_dir).mkdir(parents=True, exist_ok=True)
    # args = [
    #     cmake_executable,
    #     pbrcm_home,
    # ]
    # print(" ".join(args))
    # res = subprocess.check_output(args, cwd=pbrcm_build_dir)
    # print(res.decode('ascii'))
    # res = subprocess.check_output(make_executable, cwd=pbrcm_build_dir)
    # print(res.decode('ascii'))

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
    res = subprocess.check_output(args, cwd=graph_preprocess_dir)
    print(res.decode('utf-8'))

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
    print(res.decode('utf-8'))

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
    graph_preprocess_targets = [
        'graph_preprocess', 
        'slashburn', 
        'cuthill_mckee',
        'pr_experiments', 
        'convert_map_to_binary', 
        'compute_ccs', 
        'stats',
        'par_hilbert_bench'
    ]
    
    # compile
    for target in graph_preprocess_targets:
        args = [
            cmake_executable,
            "--build", cmake_build_dir,
            "--target", target,
            "-j", str(n_threads)
        ]
        print(" ".join(args))
        res = subprocess.check_output(args, cwd=graph_preprocess_dir)
        print(res.decode('utf-8'))

    # build and install abseil (only if not installed already)
    if not os.path.isdir(settings["abseil_install_include_dir"]):
        subprocess.check_output(
            ['git', 'clone', settings['abseil_repo_url']],
            cwd=settings["par_slashburn_dir"]
        )
        absl_install_dir = os.path.join(
            settings['par_slashburn_dir'], 'install')
        absl_build_dir = os.path.join(settings['abseil_repo_dir'], 'build')
        Path(absl_install_dir).mkdir(parents=True, exist_ok=True)
        Path(absl_build_dir).mkdir(parents=True, exist_ok=True)

        args = [
            cmake_executable,
            '..',
            f'-DCMAKE_INSTALL_PREFIX={absl_install_dir}'
        ]
        subprocess.check_output(args, cwd=absl_build_dir)
        args = [
            cmake_executable,
            '--build',
            '.',
            '--target',
            'install'
        ]
        subprocess.check_output(args, cwd=absl_build_dir)
    cmake_open_mp_options = '-DIPS4O_USE_OPENMP=ON -DONEDPL_PAR_BACKEND=openmp'
    # build and compile parallel slashburn
    # build graph_preprocess
    args = [
        cmake_executable,
        # tell cmake ips4o, onedpl to use
        f"-DCMAKE_BUILD_TYPE={cmake_build_type}",
        f"-DCMAKE_MAKE_PROGRAM={cmake_make_program}",
        '-DIPS4O_USE_OPENMP=ON',
        '-DONEDPL_PAR_BACKEND=openmp',
        "-G", "Ninja",
        "-S", par_slashburn_dir,
        "-B", os.path.join(par_slashburn_dir, cmake_build_dir),
        # path to local tbb install
        # (so that we don't have to use deprecated find_TBB cmake module
        # in ips4o submodule build)
        "-DCMAKE_PREFIX_PATH=/opt/intel/oneapi/tbb/latest"
    ]
    print(" ".join(args))
    res = subprocess.check_output(args, cwd=par_slashburn_dir)
    print(res.decode('utf-8'))

    for target in ['par_slashburn', 'pr']:
        args = [
            cmake_executable,
            "--build", os.path.join(par_slashburn_dir, cmake_build_dir),
            "--target", target,
            "-j", str(n_threads)
        ]
        print(" ".join(args))
        res = subprocess.check_output(args, cwd=par_slashburn_dir)
        print(res.decode('utf-8'))

    dbg_apps_dir = settings['dbg_apps_dir']
    # build dbg
    print(dbg_apps_dir)
    args = [
        make_executable, "clean",
    ]

    # compile dbg's graph-convert-utils convert
    args = [
        "g++", "-o"
        f"{settings['dbg_convert_dir']}/convert",
        '-fopenmp',
        f"{settings['dbg_convert_dir']}/convert.cpp",
    ]

    print(dbg_home)
    print(" ".join(args))
    subprocess.check_output(args, cwd=settings['dbg_convert_dir'])
    print(res.decode('utf-8'))

    print(" ".join(args))
    subprocess.check_output(args, cwd=dbg_apps_dir)
    print(res.decode('utf-8'))

    args = [
        make_executable,
        "-j", str(n_threads),
        f'graph_preprocess_dir={graph_preprocess_dir}',
        f'DBG_ROOT={dbg_home}',
    ]
    print(" ".join(args))
    subprocess.check_output(args, cwd=dbg_apps_dir)
    print(res.decode('utf-8'))

    prgrn_dir = settings['peregrine']['prgrn_dir']

    # build peregrine
    args = [
        make_executable,
        "-j", str(n_threads),
        'convert_data',
        'count',
    ]
    print(" ".join(args))
    subprocess.check_output(args, cwd=prgrn_dir)
    print(res.decode('utf-8'))

    # build webgraph
    webgraph_dir = settings['webgraph_dir']
    env = os.environ.copy()
    env['ANT_OPTS'] = '-Xmx256m'
    args = ['ant', 'ivy-setupjars', 'jar']
    print(" ".join(args))
    subprocess.check_output(args, cwd=webgraph_dir, env=env)
    print(res.decode('utf-8'))



if __name__ == "__main__":
    argparse_desc = """
    python script that sets up (i.e. builds and compiles) all the submodules 
    required by congenial enigma

    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--clean',
                        action=argparse.BooleanOptionalAction,
                        help='Clean a previously built build.')
    main(parser.parse_args())
