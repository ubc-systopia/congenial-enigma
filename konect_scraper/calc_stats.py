
import os
from konect_scraper import config
from konect_scraper.stats import algebraic_connectivity, diameter, get_degrees, get_eigenvalues, get_laplacian, load_csr_matrix, plfit_stats, power_law_estimate, read_nx_edgelist_to_csr, save_csr_matrix, size_of_lscc, spectral_separation, tail_power_law_estimate
import numpy as np
import networkx as nx
import powerlaw as pl
import scipy.sparse as ss
import matplotlib.pyplot as plt

def save_as_scipy_csr(graph_dir, edgelist_path):
    settings = config.settings
    scipy_csr_suffix = settings['scipy_csr_suffix']
    compressed_fname = settings['compressed_el_file_name']
    mat_path = os.path.join(
        graph_dir, f'{compressed_fname}.{scipy_csr_suffix}')
    csr_mat = read_nx_edgelist_to_csr(edgelist_path)
    save_csr_matrix(mat_path, csr_mat)

    return

def main(rows):
    settings = config.settings
    graphs_dir = settings['graphs_dir']
    compressed_fname = settings['compressed_el_file_name']
    edgelist_file_suffix = settings['edgelist_file_suffix']
    scipy_csr_suffix = settings['scipy_csr_suffix']
    print(ss.coo_matrix)
    # print(ss.coo_array)
    # compute the given orders for each of the datasets
    for _, row in enumerate(rows):
        graph_name = row['graph_name']
        print(graph_name)
        graph_dir = os.path.join(graphs_dir, graph_name)
        graph_path = os.path.join(
            graph_dir, f'{compressed_fname}.{edgelist_file_suffix}')
            
        nx_graph = nx.read_edgelist(graph_path, create_using=nx.DiGraph)
        mat_path = os.path.join(
            graph_dir, f'{compressed_fname}.{scipy_csr_suffix}')
        

        # if not os.path.isfile(mat_path):            
        save_as_scipy_csr(graph_dir, graph_path)
        
        csr_mat = load_csr_matrix(mat_path)

        # symmetrize directed csr
        symm_csr_mat = csr_mat.copy()
        rs, cs = csr_mat.nonzero()
        symm_csr_mat[cs, rs] = csr_mat[rs, cs]

        print(f'{csr_mat.todense()=}')
        print(type(csr_mat))
        dense = csr_mat.todense()
        n_rows, n_cols = csr_mat.shape
        # for i in range(n_rows):
        #     str = ''
        #     for j in range(n_cols):
        #         if dense[i, j] == 1:
        #             str += '1 '
        #         else:
        #             str += '0 '
        #     print(str)
        # print(csr_mat.nnz)
        vals, vecs = get_eigenvalues(csr_mat.astype('float'))
        egv = np.linalg.eigvals(dense)

        symm_egvals, symm_eg_vecs = get_eigenvalues(symm_csr_mat.astype('float'))
        
        u, s, vh = ss.linalg.svds(csr_mat.astype('float'))
        op_2_norm = s[-1]
        print(f'{op_2_norm=}')

        print(f'{symm_egvals=}')

        print(egv[np.argmax(np.abs(egv))])
        evals = vals
        cyclic_eval = np.abs(evals[0])
        print(f'{cyclic_eval=}')
        print(evals)
        print(spectral_separation(evals[0], evals[1]))
        print(f'{size_of_lscc(csr_mat)=}')
        print(f'{diameter(nx_graph)}')

        out_degrees = get_degrees(csr_mat, mode='out')[1]
        in_degrees = get_degrees(csr_mat, mode='in')[1]
        print(f'{out_degrees=}', f'{out_degrees.shape[0]}')
        out_pl_fit = pl.Fit(
            data=out_degrees, 
            discrete=True, 
            estimate_discrete=False,
            verbose=True,
            fit_method='Likelihood' # or KS
        )
        out_xmin = out_pl_fit.power_law.xmin
        out_alpha = out_pl_fit.power_law.alpha
        print(f'{out_degrees[out_degrees > out_xmin].shape=}')
        print(f'{out_xmin=}')
        print(f'{out_alpha=}')
        print(vars(out_pl_fit.power_law))
        print(f'{out_pl_fit.n_tail=}')

        # print(f'{power_law_estimate(out_degrees)=}')
        print(f'{tail_power_law_estimate(out_degrees, xmin=8)=}')

        print(f'{power_law_estimate(in_degrees)=}')
        pl.plot_pdf(out_degrees[out_degrees > 0])
        # plt.show()

        plfit_stats(csr_mat)

    return

