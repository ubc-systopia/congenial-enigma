import numpy as np

def init():
    global meta_col_names
    global stat_col_names
    global preproc_col_names
    global konect_col_names
    global pr_expts_col_names
    global sqlite3_to_np_dtypes
    global n_m_col_names
    global sql_to_np_dtypes

    sql_to_np_dtypes = {
        'TEXT': str,
        'STRING': str,
        'INTEGER': "Int64",
        'REAL': "Float64",
        'BIGINT': object
    }

    pr_expts_col_names = {
        "graph_name": 'TEXT',
        "datetime": 'TEXT',
        "expt_num": 'INTEGER',
        "num_iters": 'INTEGER',
        "vertex_order": 'TEXT',
        "edge_order": 'TEXT',
        "runtime": 'BIGINT'
    }

    konect_col_names = {
        "graph_name": 'TEXT',
        "name": 'TEXT',
        "code": 'TEXT',
        'n': 'INTEGER',
        'm': 'INTEGER',
        'konect_url': 'TEXT',
        'data_url': 'TEXT'
    }

    preproc_col_names = {
        'graph_name': 'TEXT',
        'compress': 'INTEGER',
        'ingest': 'INTEGER',
        'sort_vs': 'INTEGER',
        'unique_vs': 'INTEGER',
        'sort_es': 'INTEGER',
        'unique_es': 'INTEGER',
        'rabbit': 'REAL',
        'cuthill_mckee': 'REAL',
        'slashburn': 'REAL',
        'random': 'REAL',
        'hubcluster': 'REAL',
        'hubsort': 'REAL',
        'sort': 'REAL',
        'dbg': 'REAL',
    }

    n_m_col_names = {
        'graph_name': 'TEXT',
        'n': 'BIGINT',
        'm': 'BIGINT'
    }

    meta_col_names = {
        'graph_name': 'TEXT',
        'code': 'TEXT',
        'internal_name': 'TEXT',
        'name': 'TEXT',
        'directory': 'TEXT',
        'data_source': 'TEXT',
        'availability': 'TEXT',
        'consistency_check': 'TEXT',
        'category': 'TEXT',
        'node_meaning': 'TEXT',
        'edge_meaning': 'TEXT',
        'network_format': 'TEXT',
        'edge_type': 'TEXT',
        'temporal_data': 'TEXT',
        'reciprocal': 'TEXT',
        'directed_cycles': 'TEXT',
        'loops': 'TEXT',
        'snapshot': 'TEXT',
        'dataset_timestamp': 'TEXT',
        'network_join': 'TEXT',
        'orientation': 'TEXT',
        'multiplicity': 'TEXT',
        'connectedness': 'TEXT',
        'zero_weights': 'TEXT',
        'skew_symmetry': 'TEXT',
        'complete': 'TEXT',
        'k_core': 'TEXT',
        'tournament': 'TEXT',
        'paths': 'TEXT',
        'txt_file_size': 'INTEGER',
        'compressed_txt_file_size': 'INTEGER',
    }

    stat_col_names = {
        'graph_name': 'STRING',
        'directed': 'INTEGER',
        'size': 'INTEGER',
        'volume': 'INTEGER',
        'n': 'INTEGER',
        'm': 'INTEGER',
        'unique_edge_count': 'INTEGER',
        'loop_count': 'INTEGER',
        'wedge_count': 'BIGINT',
        'claw_count': 'BIGINT',
        'cross_count': 'BIGINT',
        'triangle_count': 'BIGINT',
        'square_count': 'INTEGER',
        'tour_count_4': 'INTEGER',
        'maximum_degree': 'INTEGER',
        'maximum_outdegree': 'INTEGER',
        'maximum_indegree': 'INTEGER',
        #

        'average_degree': 'REAL',
        'fill': 'REAL',
        'average_edge_multiplicity': 'REAL',
        'size_of_lcc': 'INTEGER',
        'size_of_lscc': 'INTEGER',
        'relative_size_of_lscc': 'REAL',
        'diameter': 'INTEGER',
        'percentile_effective_diameter_50': 'REAL',
        'percentile_effective_diameter_90': 'REAL',
        'median_distance': 'INTEGER',
        'mean_distance': 'REAL',
        'gini_coefficient': 'REAL',
        'balanced_inequality_ratio': 'REAL',
        'outdegree_balanced_inequality_ratio': 'REAL',
        'indegree_balanced_inequality_ratio': 'REAL',
        'relative_edge_distribution_entropy': 'REAL',
        'power_law_exponent': 'REAL',
        'tail_power_law_exponent': 'REAL',
        'tail_power_law_exponent_with_p': 'REAL',
        'p_value': 'REAL',
        'degree_assortativity': 'REAL',
        'degree_assortativity_p_value': 'REAL',
        'clustering_coefficient': 'REAL',
        'directed_clustering_coefficient': 'REAL',
        'spectral_norm': 'REAL',
        'operator_2_norm': 'REAL',
        'cyclic_eigenvalue': 'REAL',
        'reciprocity': 'REAL',
        'non_bipartivity': 'REAL',
        'normalized_non_bipartivity': 'REAL',
        'algebraic_non_bipartivity': 'REAL',
        'spectral_bipartite_frustration': 'REAL',
        'relative_controllability': 'REAL',
        'controllability': 'INTEGER',
        'in_outdegree_correlation': 'REAL',
        'outdegree_tail_power_law_exponent_with_p': 'REAL',
        'indegree_p_value': 'REAL',
        'indegree_tail_power_law_exponent_with_p': 'REAL',
        'outdegree_p_value': 'REAL',
        'algebraic_connectivity': 'REAL',
        'spectral_separation': 'REAL',
        'triadic_conflict': 'REAL',
        'spectral_signed_frustration': 'REAL',
        'algebraic_conflict': 'REAL',
        'negativity': 'REAL',

        # bipartite stats - we don't consider bipartite graphs - but konect stats exist for these
        'left_balanced_inequality_ratio': 'REAL',
        'average_left_degree': 'INTEGER',
        'right_size': 'INTEGER',
        'left_size': 'INTEGER',
        'maximum_right_degree': 'INTEGER',
        'maximum_left_degree': 'INTEGER',
        'right_balanced_inequality_ratio': 'REAL',
        'right_tail_power_law_exponent_with_p': 'REAL',
        'left_p_value': 'REAL',
        'right_p_value': 'REAL',
        'left_tail_power_law_exponent_with_p': 'REAL',
        'average_right_degree': 'INTEGER',

        # computed stats - results of running slashburn, cuthill-mckee
        'orig_bandwidth': 'INTEGER',
        'cm_bandwidth': 'INTEGER',
        'sb_k': 'BIGINT',
        'sb_n_iters': 'BIGINT',
        'pr_struct_size': 'BIGINT',

    }
