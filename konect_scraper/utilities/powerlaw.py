"""_summary_
use the powerlaw pip library, scipy, and igraph to compute:
- power law statistics 
- eigenvalue statistics
"""

import argparse
import igraph as ig
import powerlaw as pl
import time
from numba import njit, prange
import numpy as np


@njit(parallel=True)
def power_law_estimate(ds):
    # ds = ds[ds > 0]
    mn = np.min(ds)
    print(mn)
    n = ds.shape[0]
    print(n)
    return 1 + n * (np.reciprocal(np.sum(np.log(ds / mn))))

# @njit(parallel=True)
def degrees(g, mode='all'):
    deg = np.zeros(g.vcount()).astype(np.uint32)
    # have to sequentially populate degree np.array since numba not friends 
    # with igraph
    for i in range(deg.shape[0]):
        deg[i] = g.degree(i, mode=mode, loops=False)
    return deg
def main(args):
    
    graph_name = args.graph_path.split('/')[-2]
    t0 = time.time()
    g = ig.Graph.Read_Edgelist(args.graph_path, args.directed)
    t1 = time.time()
    print(f"Read {graph_name}: {g.vcount()} vertices, {g.ecount()} edges in {t1- t0} seconds")
    t0 = time.time()
    deg = degrees(g, mode='out')
    t1 = time.time()
    print(f"Computed degrees in {t1- t0} seconds")
    t0 = time.time()
    pl_fit = pl.Fit(deg, discrete=True)
    t1 = time.time()
    print(f"Fit a power law curve to the degree distribution in {t1 - t0} seconds")
    print(f'{pl_fit.power_law.alpha=}')
    print(f'{pl_fit.power_law.xmin=}')
    print(f'{power_law_estimate(deg)=}')
    print(f'{pl_fit.power_law}')

    return 

if __name__ == '__main__':
    argparse_desc = """
    Use the powerlaw pip library, numba, and igraph to compute missing
    konect statistics.
    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--graph-path', required=True, 
                        help="Absolute path to the data directory containing "
                        "a graph's edge list")
    parser.add_argument('-d', '--directed', action='store_true', 
    help='Is the graph directed or not.')

    
    main(parser.parse_args())
    