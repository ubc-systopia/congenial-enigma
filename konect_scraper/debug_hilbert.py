import subprocess

import numpy as np
import psutil
from matplotlib import pyplot as plt, cm
from matplotlib.collections import LineCollection

from konect_scraper import config


def next_largest_multiple(n, d):
    multiple = np.power(2, d)
    # print(multiple)
    return int((n + multiple - 1) / multiple) * multiple
def hilbert_debug(n, n_threads):
    hilbert_debug_exec = "/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/sb_furhilbert"

    args = [
        hilbert_debug_exec, "-n", str(n),
        # "-m", str(n_threads),
    ]
    res = subprocess.check_output(args)
    print(res.decode('ascii'))


"""
uint32_t next_largest_multiple(uint32_t n, uint32_t critical_depth) {
	assert(critical_depth > 0);
	uint32_t multiple = pow(2, critical_depth);
	return ((n + multiple - 1) / multiple) * multiple;
}"""




def plot_hilbert(input_path, output_path, n):
    mat = np.zeros((n, n))

    fig, ax = plt.subplots()

    ax.matshow(mat)

    lines = []

    order = np.loadtxt(input_path).astype(np.uint32)



    for previous, current in zip(order, order[1:]):
        lines.append([
            [previous[1], previous[0]],
            [current[1], current[0]],
        ])

    color = iter(cm.rainbow(np.linspace(0, 1, len(lines))))
    colors = []
    for i in range(len(lines)):
        colors.append(next(color))
    width = 2
    lc = LineCollection(
        lines,
        linewidths=width,
        colors=colors
    )
    ax.add_collection(lc)
    ax.axis('off')
    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight', pad_inches=0, dpi=200)
    plt.close(fig)

def get_critical_depth():
    n_threads = config.settings['n_threads']
    n_quads = 1
    critical_depth = 0
    while n_quads < n_threads:
        critical_depth += 1
        n_quads = n_quads * 4
    return critical_depth
def main():
    hilbert_debug_path = "/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/hilbert"
    hilbert_plot_path = "/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/hilbert.png"
    n = 14
    critical_depth = 2
    critical_depth = get_critical_depth()
    n_threads = 2
    n = next_largest_multiple(n, critical_depth + 2)
    print(n)
    hilbert_debug(n, n_threads)
    plot_hilbert(hilbert_debug_path, hilbert_plot_path, n)

    return


if __name__ == "__main__":
    main()

"""
    for n in range(60, 70):
        hilbert_debug(n, n - 1)
        hilbert_debug(n, n)
        hilbert_debug(n, n + 1)

        i1 = f"{hilbert_debug_path}_{n}_{n - 1}"
        i2 = f"{hilbert_debug_path}_{n}_{n}"
        i3 = f"{hilbert_debug_path}_{n}_{n + 1}"

        o1 = f"{i1}.png"
        o2 = f"{i2}.png"
        o3 = f"{i3}.png"

        for i, o in zip([i1, i2, i3], [o1, o2, o3]):
            plot_hilbert(i, o, n)"""
