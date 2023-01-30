import itertools
import os
from typing import Callable
import numpy as np
import scipy.io

# Note: Generate all the functions by passing these

class ConnectivityMatrix():
	def __init__(self, num_pen: int, num_peg: int, num_epg: int, num_d7: int, loader: Callable[[], np.ndarray]) -> None:
		self.NUM_PEN = num_pen
		self.NUM_PEG = num_peg
		self.NUM_EPG = num_epg
		self.NUM_D7 = num_d7
		self.CONNECTIVITY_MATRIX_SIZE = self.NUM_PEN + self.NUM_PEG + self.NUM_EPG + self.NUM_D7

		self.SLICE_PEN = slice(0, self.NUM_PEN)
		self.SLICE_PEG = slice(self.NUM_PEN,self.NUM_PEN+self.NUM_PEG)
		self.SLICE_EPG = slice(self.NUM_PEN+self.NUM_PEG,self.NUM_PEN+self.NUM_PEG+self.NUM_EPG)
		self.SLICE_D7 = slice(self.NUM_PEN+self.NUM_PEG+self.NUM_EPG,self.NUM_PEN+self.NUM_PEG+self.NUM_EPG+self.NUM_D7)

		self.GET_MATRIX = loader

slice_to_range = lambda slice: range(slice.start if slice.start else 0, slice.stop, slice.step if slice.step else 1)

def get_fly_connectivity_matrix_raw():
	return scipy.io.loadmat(os.path.join(os.path.dirname(__file__), 'connectivity_matrix_drosophila_mine_case_5_9cols_labels1.mat'))['con_matrix'].T

def get_locust_connectivity_matrix_raw():
	return scipy.io.loadmat(os.path.join(os.path.dirname(__file__), 'connectivity_matrix_locust_mine_case_6_with3s_labels_1.mat'))['con_matrix'].T

FLY = ConnectivityMatrix(16, 18, 18, 8, get_fly_connectivity_matrix_raw)
LOCUST = ConnectivityMatrix(16, 16, 16, 8, get_locust_connectivity_matrix_raw)

def get_connectivity_matrix_bool(CM: ConnectivityMatrix):
	return CM.GET_MATRIX() != 0

get_fly_connectivity_matrix_bool = lambda : get_connectivity_matrix_bool(FLY)
get_locust_connectivity_matrix_bool = lambda : get_connectivity_matrix_bool(LOCUST)

def get_connectivity_matrix_string(CM: ConnectivityMatrix):
	'''Get the connectivity matrix with each cell either 0 (no connection) or `'connections_{pre_name}_{post_name}'` if the neurons are connected.'''
	connectivity_matrix = CM.GET_MATRIX()
	conn_mat_str = np.zeros_like(connectivity_matrix, dtype=object)

	for (pre_slice, pre_name), (post_slice, post_name) in itertools.product(zip([CM.SLICE_PEN, CM.SLICE_PEG, CM.SLICE_EPG, CM.SLICE_D7], ['PEN', 'PEG', 'EPG', 'D7']), repeat=2):
		conn_mat_str[pre_slice,post_slice][np.nonzero(connectivity_matrix[pre_slice,post_slice])] = f'connections_{pre_name}_{post_name}'

	return conn_mat_str

get_fly_connectivity_matrix_string = lambda : get_connectivity_matrix_string(FLY)
get_locust_connectivity_matrix_string = lambda : get_connectivity_matrix_string(LOCUST)

def neuron_index_to_name(CM: ConnectivityMatrix, neuron_index: int):
	assert neuron_index >= 0 and neuron_index < CM.CONNECTIVITY_MATRIX_SIZE, f'Neuron index {neuron_index} out of range of connectivity matrix'

	for index_range, name in zip([slice_to_range(slice) for slice in [CM.SLICE_PEN, CM.SLICE_PEG, CM.SLICE_EPG, CM.SLICE_D7]], ['PEN', 'PEG', 'EPG', 'D7']):
		if neuron_index in index_range:
			return f'{name}_{index_range.index(neuron_index)}'
	
	raise IndexError(f"Can't find index {neuron_index} in connectivity matrix.")

fly_neuron_index_to_name = lambda neuron_index : neuron_index_to_name(FLY, neuron_index)
locust_neuron_index_to_name = lambda neuron_index : neuron_index_to_name(LOCUST, neuron_index)

def get_connectivity_matrix_edge_list(CM: ConnectivityMatrix):
	connectivity_matrix = CM.GET_MATRIX()
	index_to_name = lambda neuron_index : neuron_index_to_name(CM, neuron_index)

	return [(index_to_name(pre), index_to_name(post), {'weight': 5 if index_to_name(pre).split('_')[0]=='D7' else 1}) \
		for pre in range(CM.CONNECTIVITY_MATRIX_SIZE) \
		for post in range(CM.CONNECTIVITY_MATRIX_SIZE) \
		if connectivity_matrix[pre,post]]

get_fly_connectivity_matrix_edge_list = lambda : get_connectivity_matrix_edge_list(FLY)
get_locust_connectivity_matrix_edge_list = lambda : get_connectivity_matrix_edge_list(LOCUST)