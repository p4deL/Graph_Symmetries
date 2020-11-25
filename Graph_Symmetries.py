# Graph_Symmetries.py
#
# Created on: 2020-11-23
# Author: padelhardt
#
# Add program describtion here
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
	"""Read in the number of vertices, number fo extra info bits and graph object of bond file."""

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


def create_degree_hop_dict(graph, is_observable=False):
	degree_hop_dict = {}

	# iterate over all nodes to set start and target node
	for node_start in graph.nodes:

		if is_observable:
			target_nodes = list(graph.nodes)
		else:
			target_nodes = list(graph.nodes)[node_start:]

		for node_target in target_nodes:
			tmp = graph.copy()

			tmp.nodes[node_start]["hop"] = 1
			tmp.nodes[node_target]["hop"] = 2

			degrees = (tmp.degree[node_start], tmp.degree[node_target])

			# check current degrees tuple already exists as key
			if degrees in degree_hop_dict.keys():
				degree_hop_dict[degrees].append(((node_start, node_target), tmp))
			else:
				hop_list = [((node_start, node_target), tmp)]
				degree_hop_dict[degrees] = hop_list

	return degree_hop_dict


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
		degree_hop_dict = create_degree_hop_dict(graph, is_observable=is_observable)
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
	# traverse all graphs
	sym_list = []

	for graph_name, graph in zip(graph_name_list, graph_list):

		# create symmetry matrix where hoppings that can be matched will have the same entry.
		# Entries are simply given in ascending order
		sz = graph.number_of_nodes()
		sym_matrix = np.zeros((sz, sz), dtype=int)

		# create 2 dimensional dictionary where first key is (degree start, degree target), second key is
		# (node start, node target) and the value if a list of accordingly colored graphs
		degree_hop_dict = create_degree_hop_dict(graph, is_observable=is_observable)
		# print("graph: ", graph.graph)

		sym_counter = 0

		# iterate through outer dictionary with degree keys. As only hopping between nodes of the same degree can be matched
		for degrees, colored_graphs in degree_hop_dict.items():

			# iterate through inner dictionary
			# choose a colored graph as reference graph
			#print(graph.graph)

			for i, element in enumerate(colored_graphs):
				(ref_start, ref_target), ref_graph = element
				#print(ref_graph.nodes.data())

				if not ref_graph.graph["checked"]:
					ref_graph.graph["checked"] = True
					sym_counter += 1
					#print("check first:", (ref_start, ref_target), sym_counter)
					sym_matrix[ref_start, ref_target] = sym_counter
					if not is_observable and ref_start != ref_target:
						sym_matrix[ref_target, ref_start] = sym_counter

					for (check_start, check_target), check_graph in colored_graphs[i+1:]:

						if not check_graph.graph["checked"]:
							nm = iso.categorical_node_match("hop", 0)
							if nx.is_isomorphic(ref_graph, check_graph, node_match=nm):
								check_graph.graph["checked"] = True
								#sym_counter += 1
								#print("check other:", (check_start, check_target), sym_counter)
								sym_matrix[check_start, check_target] = sym_counter
								if not is_observable:
									sym_matrix[check_target, check_start] = sym_counter



		#print(sym_matrix)
		#print("\n")

		for n in range(sym_matrix.size):
			indices = np.argwhere(sym_matrix == n)
			sym_numb_local = indices.shape[0]
			if sym_numb_local != 0:
				start, target = indices[0]
				sym_list.append((graph_name.replace('\n', ''), start, target, sym_numb_local))

	return sym_list


def usage():
	return "Usage: Graph_Symmetries.py -g <path to graph list>"


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
	path_to_graph_list, observable_flag = param_use(argv)
	dir = os.path.dirname(path_to_graph_list) + "/"

	# read in all graph names
	with open(path_to_graph_list, 'r') as graph_list_file:
		graph_name_list = graph_list_file.readlines()

	graph_list = []

	# read in all graphs from the bond files
	for graph_name in graph_name_list:
		graph_name = graph_name.replace('\n', '')
		graph_key = re.search('_(.*)_order', graph_name).group(1)
		glob_sym_numb, g = read_bondfile(dir + graph_name)
		g.graph['key'] = graph_key
		g.graph['sym_numb'] = glob_sym_numb
		g.graph['checked'] = False
		graph_list.append(g)
	# print(g.graph, g.nodes, g.edges)

	start_time = time.time()
	sym_list = calc_symmetries_faster(graph_name_list, graph_list, is_observable=True)
	print("--- %s seconds ---" % (time.time() - start_time))

	# print(sym_list)

	if observable_flag:
		out_file_name = "output/_0_symObsGraphsList_order4.cfg"
	else:
		out_file_name = "output/_0_symGraphsList_order4.cfg"
	with open(out_file_name, 'w') as out_file:
		out_file.write("#graph_number  hopping_from hopping_to symmetry_number \n")
		for entry in sym_list:
			graph_name, start, target, sym_numb_local = entry
			out_file.write(graph_name + ' ' + str(start) + ' ' + str(target) + ' ' + str(sym_numb_local) + '\n')


if __name__ == "__main__":
	main(sys.argv[1:])
