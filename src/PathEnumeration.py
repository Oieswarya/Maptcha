# Maptcha: An efficient parallel workflow for hybrid genome scaffolding

# Oieswarya Bhowmik, Tazin Rahman, Ananth Kalyanaraman

#      (oieswarya.bhowmik@wsu.edu, tazin.rahman@wsu.edu, ananth@wsu.edu)

# Washington State University

#

# **************************************************************************************************

# Copyright (c) 2024. Washington State University ("WSU"). All Rights Reserved.
# Permission to use, copy, modify, and distribute this software and its documentation
# for educational, research, and not-for-profit purposes, without fee, is hereby
# granted, provided that the above copyright notice, this paragraph and the following
# two paragraphs appear in all copies, modifications, and distributions. For
# commercial licensing opportunities, please contact The Office of Commercialization,
# WSU, 280/286 Lighty, PB Box 641060, Pullman, WA 99164, (509) 335-5526,
# commercialization@wsu.edu<mailto:commercialization@wsu.edu>, https://commercialization.wsu.edu/

# IN NO EVENT SHALL WSU BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
# OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
# THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF WSU HAS BEEN ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# WSU SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND
# ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". WSU HAS NO
# OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
# **************************************************************************************************



import sys
import networkx as nx
import random
from tqdm import tqdm
import time
from collections import defaultdict

def load_interactions(filename):
    G = nx.Graph()
    with open(filename, 'r') as file:
        for line in file:
            node_a, node_b, score = line.strip().split()
            node_a, node_b = str(node_a), str(node_b)  # Convert to strings
            G.add_edge(node_a, node_b)
    return G

def load_attributes(filename, G):
    with open(filename, 'r') as file:
        for line in file:
            node_attributes = line.strip().split()
            node = str(node_attributes[0])  # Convert to string
            attributes = node_attributes[1:]

            # Check if the node exists in the graph before adding attributes
            if G.has_node(node):
                G.nodes[node]['attribute'] = attributes
            else:
                print(f"Node {node} not found in the graph.")

    return G

def compute_walk_progress(G):
    # Check if graph is empty
    if not G.nodes():
        print("The graph is empty.")
        return []

    # Find nodes with degree 1
    nodes_with_degree_one = [node for node, degree in G.degree() if degree == 1]

    # If there are no nodes with degree 1, choose a node with the minimum degree
    if not nodes_with_degree_one:
        min_degree = min(dict(G.degree()).values())
        min_degree_nodes = [node for node, degree in G.degree() if degree == min_degree]
        start_node = random.choice(min_degree_nodes)
    else:
        start_node = random.choice(nodes_with_degree_one)

    paths = []

    total_edges = len(G.edges())
    #pbar = tqdm(total=100)

    while G.edges():
        path = [start_node]
#         print(start_node, end='   ')
        while G.neighbors(start_node):
            neighbors = list(G.neighbors(start_node))

            # Removing the node from which we came
            if len(path) > 1 and path[-2] in neighbors:
                neighbors.remove(path[-2])
                
            # If there are no neighbors left, we've reached the end of the path.
            if not neighbors:
                break
            
            if len(neighbors) == 1:
                next_node = neighbors[0]
            else:
                # Find the neighbor with the highest intersection of attributes with both the previous and current nodes
                start_node_attributes = set(G.nodes[start_node].get('attribute', []))
                if len(path) > 1:
                    prev_node_attributes = set(G.nodes[path[-2]].get('attribute', []))
                else:
                    prev_node_attributes = set()

                intersection_counts = [max(len(start_node_attributes & set(G.nodes[neighbor].get('attribute', []))),
                                           len(prev_node_attributes & set(start_node_attributes)))
                                       for neighbor in neighbors]
                next_node = neighbors[intersection_counts.index(max(intersection_counts))]

            path.append(next_node)
#             print(next_node, end= '   ')
            # Delete the edge already used
            G.remove_edge(start_node, next_node)
            start_node = next_node

            # Update progress bar
            #pbar.update((1 - len(G.edges()) / total_edges) * 100 - pbar.n)

        paths.append(path)
#         print('got a path', path)

        # If there are still edges left, select a new start node
        if G.edges():
            min_degree = min(degree for node, degree in G.degree() if degree > 0)
            min_degree_nodes = [node for node, degree in G.degree() if degree == min_degree]
            start_node = random.choice(min_degree_nodes)

    #pbar.close()
    return paths

def time_compute_walk(G):
    start_time = time.time()
    paths = compute_walk_progress(G)
    end_time = time.time()
    #print(f"The computation took {end_time - start_time} seconds.")
    return paths

def transform_attributes_file(input_filename, output_filename):
    attributes_by_node = defaultdict(list)

    with open(input_filename, 'r') as input_file:
        next(input_file)  # skip first line
        for line in input_file:
            node, attribute = line.strip().split()
            attributes_by_node[node].append(attribute)

    with open(output_filename, 'w') as output_file:
        for node, attributes in attributes_by_node.items():
            output_file.write(node + " " + " ".join(attributes) + "\n")

def write_paths_to_file(paths, filename):
    with open(filename, 'w') as file:
        for path in paths:
            #file.write(' '.join(path) + '\n')
            file.write(",".join(path) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 PathEnumeration.py transform_attributes_file G output_file")
        sys.exit(1)

    transform_attributes_file(sys.argv[1], 'attributes.txt')

    G = load_interactions(sys.argv[2])
    G = load_attributes('attributes.txt', G)
    paths = time_compute_walk(G)

    write_paths_to_file(paths, sys.argv[3])
