import numpy as np
import scipy.io

def get_fly_connectivity_matrix_raw():
	return scipy.io.loadmat('Insect Head Direction Network/connectivity_matrix_drosophila_mine_case_5_9cols_labels1.mat')['con_matrix'].T

def get_fly_connectivity_matrix_bool():
	return get_fly_connectivity_matrix_raw() != 0

def get_fly_connectivity_matrix_string():
	connectivity_matrix = get_fly_connectivity_matrix_raw()

	PEN = slice(16)
	PEG = slice(16,16+18)
	EPG = slice(16+18,16+18+18)
	D7 = slice(16+18+18,16+18+18+8)
	connections_PEN_EPG = np.nonzero(connectivity_matrix[PEN,EPG])
	connections_PEG_EPG = np.nonzero(connectivity_matrix[PEG,EPG])
	connections_EPG_PEN = np.nonzero(connectivity_matrix[EPG,PEN])
	connections_EPG_PEG = np.nonzero(connectivity_matrix[EPG,PEG])
	connections_EPG_D7 = np.nonzero(connectivity_matrix[EPG,D7])
	connections_D7_PEN = np.nonzero(connectivity_matrix[D7,PEN])
	connections_D7_PEG = np.nonzero(connectivity_matrix[D7,PEG])
	connections_D7_D7 = np.nonzero(connectivity_matrix[D7,D7])

	conn_mat_str = np.zeros_like(connectivity_matrix, dtype=object)
	conn_mat_str[PEN,EPG][connections_PEN_EPG] = 'connections_PEN_EPG'
	conn_mat_str[PEG,EPG][connections_PEG_EPG] = 'connections_PEG_EPG'
	conn_mat_str[EPG,PEN][connections_EPG_PEN] = 'connections_EPG_PEN'
	conn_mat_str[EPG,PEG][connections_EPG_PEG] = 'connections_EPG_PEG'
	conn_mat_str[EPG,D7][connections_EPG_D7] = 'connections_EPG_D7'
	conn_mat_str[D7,PEN][connections_D7_PEN] = 'connections_D7_PEN'
	conn_mat_str[D7,PEG][connections_D7_PEG] = 'connections_D7_PEG'
	conn_mat_str[D7,D7][connections_D7_D7] = 'connections_D7_D7'

	return connectivity_matrix