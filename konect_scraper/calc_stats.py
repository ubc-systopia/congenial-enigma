
import os
from konect_scraper import config
from konect_scraper.stats import algebraic_connectivity, get_eigenvalues, get_laplacian, load_csr_matrix, read_data_file_as_coo_matrix, read_nx_edgelist_to_csr, save_csr_matrix, spectral_separation
import numpy as np

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

    # compute the given orders for each of the datasets
    for _, row in enumerate(rows):
        graph_name = row['graph_name']
        print(graph_name)
        graph_dir = os.path.join(graphs_dir, graph_name)
        graph_path = os.path.join(
            graph_dir, f'{compressed_fname}.{edgelist_file_suffix}')
            # graph_dir, f'orig.{edgelist_file_suffix}')
        mat_path = os.path.join(
            graph_dir, f'{compressed_fname}.{scipy_csr_suffix}')
        

        # if not os.path.isfile(mat_path):            
        save_as_scipy_csr(graph_dir, graph_path)
        
        csr_mat = load_csr_matrix(mat_path)
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
        print(egv[np.argmax(np.abs(egv))])
        evals = vals
        spectral_norm = np.abs(evals[0])
        print(f'{spectral_norm=}')
        print(evals)
        print(spectral_separation(evals[0], evals[1]))
        print(algebraic_connectivity(csr_mat))

    return
