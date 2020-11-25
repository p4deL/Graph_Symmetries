# Graph_Symmetries.py
#
# Created on: 2020-11-23
# Author: padelhardt
#
# From a graph list and bond files this program reads in the graphs, checks if processes on these graph are isomorphic
# to each other and counts the associated symmetry numbers. The output is a file containing a list of reduced processes
# and the associated symmetry numbers. There is also the option to account for a "reverse process" symmetry, i.e.,
# i->j = j->I.
#
#

import sys
import getopt
import os.path
import re
import time
import numpy as np
import networkx as nx
import networkx.algorithms.isomorphism as iso


def read_bondfile(name):
	"""Read in the number of vertices, number fo extra info bits and graph object of bond file.
	:param name: full path name of the bond file
	:return: (global) symmetry number, graph object
	"""

	read_line = symmetry_number = 0
	edge_list = []

	with open(name, "r") as file:
		for line in file.readlines():
			if line.startswith("# Symmetry number"):
				symmetry_number = int(line.replace("# Symmetry number = ", ""))
			elif not (line.startswith("#") or line.startswith("\n")):
				if read_line > 1:
					line = line.replace(";\n", "")
					bond = line.split(" ")
					edge = (int(bond[1]), int(bond[2]))
					edge_list.append(edge)

				read_line += 1

	graph = nx.Graph()

	graph.add_edges_from(edge_list)

	return symmetry_number, graph


def create_graph_list(dir, graph_name_list):
	"""Create a list of graph objects containing the necessary information and determine the maximal graph order.
	:param dir: path to directory containing the graph list
	:param graph_name_list: list name containing the graph names
	:return: list containing the graphs, (maximal) graph order
	"""

	graph_list = []
	order = 0

	for graph_name in graph_name_list:
		# extract graph key and order
		graph_name = graph_name.replace('\n', '')
		graph_key = re.search('_(.*)_order', graph_name).group(1)
		order_graph = int(re.search('_order(.*).cfg', graph_name).group(1))
		# find maximal order of given list
		if order_graph > order:
			order = order_graph

		# create graphs and assign global properties
		glob_sym_numb, g = read_bondfile(dir + graph_name)
		g.graph['key'] = graph_key
		g.graph['sym_numb'] = glob_sym_numb  # global symmetry number
		g.graph['checked'] = False  # we will later use this property to check if a graph was sorted out due to a sym.
		graph_list.append(g)

	return graph_list, order


def create_degree_process_dict(graph, is_observable=False):
	"""Create colored graph to incorporate start and target nodes of a process and sort the colored graphs into a
	dictionary with (degree start, degree target) keys and list values containing tuples of the form
	(node start, node target), colored graph).
	:param graph: uncolored graph which is used to generate colored graphs
	:param is_observable: If the flag is disabled the number of graphs is reduced because of the symmetry i->j = j->i
	:return: created
	"""

	degree_process_dict = {}

	# iterate over all nodes to set start and target node
	for node_start in graph.nodes:

		# if we calculate an observable we need to differentiate i->j and j->i processes
		if is_observable:
			target_nodes = list(graph.nodes)
		else:  # otherwise not, so we confine to the cas j>=i
			target_nodes = list(graph.nodes)[node_start:]

		for node_target in target_nodes:
			tmp = graph.copy()

			# we color the nodes with the attribute "hop". 1 indicates the start and 2 the target node.
			tmp.nodes[node_start]["hop"] = 1
			tmp.nodes[node_target]["hop"] = 2

			# extract the degrees of the considered nodes since we want to use such tuples as keys in the dictionary
			degrees = (tmp.degree[node_start], tmp.degree[node_target])

			# check current degrees tuple already exists as key, if yes append to the value list of the dictionary otherwise
			# create a new entry in the dictionary
			if degrees in degree_process_dict.keys():
				degree_process_dict[degrees].append(((node_start, node_target), tmp))  # for later convenience also store nodes
			else:
				process_list = [((node_start, node_target), tmp)]
				degree_process_dict[degrees] = process_list

	return degree_process_dict


# TOOOOOOO SLOW
def calc_symmetries(graph_name_list, graph_list, is_observable=False):
	# traverse all graphs
	sym_list = []

	for graph_name, graph in zip(graph_name_list, graph_list):

		# create symmetry matrix where hoppings that can be matched will have the same entry.
		# Entries are simply given in ascending order
		sz = graph.number_of_nodes()
		sym_matrix = np.zeros((sz, sz), dtype=int)

		# create 2 dimensional dictionary where first key is (degree start, degree target), second key is
		# (node start, node target) and the value if a list of accordingly colored graphs
		degree_hop_dict = create_degree_process_dict(graph, is_observable=is_observable)
		# print("graph: ", graph.graph)

		sym_counter = 0

		# iterate through outer dictionary with degree keys. As only hopping between nodes of the same degree can be matched
		for degrees, colored_graphs in degree_hop_dict.items():
			# print("	degrees: ", degrees)

			# These must be new symmetries as we have a new degree tuple

			# iterate through inner dictionary
			# choose a colored graph as reference graph
			for idx, element in enumerate(colored_graphs):
				(ref_start, ref_target), ref_graph = element
				# for each reference graph increase the symmetry counter as this might be a new symmetry
				# sym_counter += 1
				found_sym = False

				# choose a second graph to compare with the reference graph
				# print("		hops: ", ref_start, ref_target)
				for (check_start, check_target), check_graph in colored_graphs[idx:]:
					# print("			", check_start, check_target)

					# check if the two (colored) graphs are isomorphic. Therefore consider the node attributes to match the start
					# and end nodes
					# nm = iso.categorical_node_match(["start", "target"], [False, False])
					nm = iso.numerical_node_match("hop", 0)
					if nx.is_isomorphic(ref_graph, check_graph, node_match=nm):
						# print("				ok")
						if sym_matrix[check_start, check_target] == 0:
							if not found_sym:
								found_sym = True
								sym_counter += 1
							sym_matrix[check_start, check_target] = sym_counter
							if not is_observable:
								sym_matrix[check_target, check_start] = sym_counter

		print(graph.graph)
		print(sym_matrix)
		print("\n")

		for n in range(sym_matrix.size):
			indices = np.argwhere(sym_matrix == n)
			sym_numb_local = indices.shape[0]
			if sym_numb_local != 0:
				start, target = indices[0]
				sym_list.append((graph_name.replace('\n', ''), start, target, sym_numb_local))

	return sym_list


def calc_symmetries_faster(graph_name_list, graph_list, is_observable=False):
	"""Calculate the symmetry numbers of all processes for each graph in a graph list to effectively reduce the number
	of processes.
	:param graph_name_list: list containing the graph names
	:param graph_list: list containing the graph objects
	:param is_observable: If the flag is disabled the number of graphs is reduced because of the symmetry i->j = j->i
	:return: symmetry list containing tuples of the form (graph name, global symmetry, start, target, local symmetry)
	"""
	sym_list = []

	# iterate over all graphs
	for graph_name, graph in zip(graph_name_list, graph_list):

		# create dummy symmetry matrix where processes (e.g. hoppings) that can be matched will have the same entry.
		# Entries will later be given in  ascending order
		sz = graph.number_of_nodes()
		sym_matrix = np.zeros((sz, sz), dtype=int)

		# create dictionary where first key is (degree start, degree target) and the value is a list of the form
		# [(node start, node target), colored graph), (...), ...].
		degree_hop_dict = create_degree_process_dict(graph, is_observable=is_observable)

		# we use this to number the symmetries
		sym_counter = 0

		# iterate through dictionary with degree keys. As only hopping between nodes of the same degree can be matched
		for degrees, colored_graphs in degree_hop_dict.items():

			# iterate through colored graphs
			for i, ((ref_start, ref_target), ref_graph) in enumerate(colored_graphs):

				# At first, we check if a colored graph needs to be ignored because it is isomorphic to a graph that
				# was already checked. If not it is a new symmetry
				if not ref_graph.graph["checked"]:
					ref_graph.graph["checked"] = True
					sym_counter += 1
					# write symmetry number to the symmetry matrix
					sym_matrix[ref_start, ref_target] = sym_counter
					# if it is a hopping we have this symmetry: i->j = j->i
					if not is_observable and ref_start != ref_target:
						sym_matrix[ref_target, ref_start] = sym_counter

					# iterate over all other graphs
					for (check_start, check_target), check_graph in colored_graphs[i + 1:]:

						# if the current graph is unchecked, we need to check if it is isomorphic to the ref_graph
						if not check_graph.graph["checked"]:
							nm = iso.categorical_node_match("hop", 0)
							# in case we found an isomorphism we mark the isomorphic graph as checked
							if nx.is_isomorphic(ref_graph, check_graph, node_match=nm):
								check_graph.graph["checked"] = True
								# write to the symmetry matrix
								sym_matrix[check_start, check_target] = sym_counter
								if not is_observable:
									sym_matrix[check_target, check_start] = sym_counter

		# At the end we evaluate the symmetry matrix, calculate the symmetry factor and append the info to the sym_list
		for n in range(sym_matrix.size):
			indices = np.argwhere(sym_matrix == n)
			sym_numb_local = indices.shape[0]
			if sym_numb_local != 0:
				start, target = indices[0]
				sym_list.append((graph_name.replace('\n', ''), graph.graph["sym_numb"], start, target, sym_numb_local))

	return sym_list


def usage():

	return "Usage: Graph_Symmetries.py -g <path to graph list> [-o]"


def param_use(argv):

	found_g = False
	observable_flag = False
	path_to_graph_list = ""

	try:
		opts, args = getopt.getopt(argv, 'g:o', ["graphs=", "observable"])
	except getopt.GetoptError:
		usage()
		sys.exit(1)

	for opt, arg in opts:
		if opt in ('-g', '--graphs'):
			found_g = True
			path_to_graph_list = arg
		if opt in ('-o', '--observable'):
			observable_flag = True

	if not found_g:
		print("graph list not found.")
		usage()
		sys.exit(1)

	return path_to_graph_list, observable_flag


def main(argv):

	# load terminal parameters
	path_to_graph_list, observable_flag = param_use(argv)
	dir = os.path.dirname(path_to_graph_list) + "/"

	# read in all graph names
	with open(path_to_graph_list, 'r') as graph_list_file:
		graph_name_list = graph_list_file.readlines()

	# read in all graphs from the bond files
	graph_list, order = create_graph_list(dir, graph_name_list)

	# calculate the symmetries of every graph
	start_time = time.time()
	sym_list = calc_symmetries_faster(graph_name_list, graph_list, is_observable=observable_flag)
	print("--- %s seconds ---" % (time.time() - start_time))

	# different file names for hopping or observable symmetries
	if observable_flag:
		out_file_name = "output/_0_symObsGraphsList_order" + str(order) + ".cfg"
	else:
		out_file_name = "output/_0_symGraphsList_order" + str(order) + ".cfg"

	# write symmetry factors to file
	with open(out_file_name, 'w') as out_file:
		out_file.write("#graph_number  hopping_from hopping_to symmetry_number \n")
		for entry in sym_list:
			graph_name, global_sym_numb, start, target, sym_numb_local = entry
			out_file.write(graph_name + ' ' + str(global_sym_numb) + ' ' + str(start) + ' ' + str(target) + ' '
										 + str(sym_numb_local) + '\n')


if __name__ == "__main__":
	main(sys.argv[1:])
