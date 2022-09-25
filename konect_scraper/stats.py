import numpy as np
from scipy.special import comb

from konect_scraper.sql import is_bipartite, single_val_get
from konect_scraper.util import get_n_m, get_directed


def verify_stat(col, graph_name):
    """
    If the statistic "col" needs to be recomputed:
      - either Null, 0 (where inappropriate), or is different from a manually computed value
    Recompute it.
    Otherwise, return the existing Konect stat
    """
    match col:
        case 'fill':
            return verify_float(col, graph_name)


def verify_float(col, graph_name):
    existing_val = single_val_get(col, 'statistics', graph_name)
    computed_val = float(calc_stat(col, graph_name))
    print(f"{existing_val=}, {computed_val=}, {graph_name=}, {is_bipartite(graph_name)}")
    if existing_val == 0:
        return computed_val
    # elif not np.allclose(computed_val, existing_val): trust konect
    #     return computed_val
    elif is_bipartite(graph_name):
        return existing_val
    else:
        return existing_val


def calc_stat(col, graph_name):
    match col:
        case 'fill':
            return fill(graph_name)


def loops_allowed(graph_name):
    loops = single_val_get('loops', 'metadata', graph_name)
    if loops:
        return loops == 'Contains loops'
    else:
        return loops


def fill(graph_name):
    n, m = get_n_m(graph_name)
    directed = get_directed(graph_name)
    loops = loops_allowed(graph_name)
    if directed:
        n_possible_edges = comb(n, 2, exact=True)
    else:
        n_possible_edges = comb(n, 2, exact=True)
    if loops:
        n_possible_edges += n

    return m / n_possible_edges
