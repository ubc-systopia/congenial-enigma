"""
In order to verify correctness of parallel slashburn implementation, use ground truth connected components
computed using python-igraph before and after the removal of a hub vertex
"""

import igraph as ig
import numpy as np
def main():
    graph_path = "/media/atrostan/patterson_backup/data/graphs/maayan-pdzbase/comp.net"
    # read graph
    g = ig.Graph.Read_Edgelist(graph_path, directed=False)
    g.vs['orig_id'] = [v.index for v in g.vs]
    degrees = np.array([v.degree() for v in g.vs])
    clusters = g.clusters(mode='weak')
    print(f"{degrees=}")
    hub_id = np.argmax(degrees)
    print(f"{hub_id=}")
    for c in clusters:
        mapped_c = [g.vs[v]['orig_id'] for v in c]
        print(f"{mapped_c=} {len(mapped_c)=}")
    print("\nremove hub\n")
    g.delete_vertices([hub_id])

    clusters = g.clusters(mode='weak')
    for c in clusters:
        mapped_c = [g.vs[v]['orig_id'] for v in c]
        print(f"{mapped_c=} {len(mapped_c)}")



    # compute connected components

    # remove hub

    # compute connected components

    return

if __name__ == "__main__":
    main()