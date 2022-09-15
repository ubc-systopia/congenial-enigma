import subprocess

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection


def hilbert_debug(n):
    hilbert_debug_exec = "/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/sb_furhilbert"

    args = [
        hilbert_debug_exec, "-n", str(n)
    ]
    subprocess.check_output(args)


def rename_files(dir):
    import os
    path = '/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/'
    for filename in os.listdir(path):
        print(filename)
        num, name = filename[:-4].split('_')
        num = num.zfill(4)
        new_filename = num + "_" + name + ".png"
        os.rename(os.path.join(path, filename), os.path.join(path, new_filename))


def main():
    hilbert_debug_path = "/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/hilbert"

    rename_files("/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/")
    # for n in range(450, 550):
    #     # n = 120
    #
    #     hilbert_debug(n)
    #
    #     order = np.loadtxt(hilbert_debug_path).astype(np.uint32)
    #     mat = np.zeros((n, n))
    #
    #     fig, ax = plt.subplots()
    #
    #     ax.matshow(mat)
    #
    #     lines = []
    #     for previous, current in zip(order, order[1:]):
    #         lines.append([
    #             [previous[1], previous[0]],
    #             # [previous[0], previous[1]],
    #             [current[1], current[0]],
    #             # [current[0], current[1]],
    #         ])
    #     width = 1
    #     lc = LineCollection(
    #         lines,
    #         linewidths=width
    #     )
    #     ax.add_collection(lc)
    #     ax.axis('off')
    #     plt.tight_layout()
    #     fig.savefig(
    #         f'/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/{n}_hilbert.png', bbox_inches='tight', pad_inches=0,
    #         dpi=200)
    #     plt.close(fig)

    return


if __name__ == "__main__":
    main()
