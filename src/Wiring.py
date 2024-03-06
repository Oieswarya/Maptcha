import sys
import networkx as nx

def wire_contigs(file1, file2, output_file):
    # create a dictionary mapping contigs to the set of long reads that map to them
    contig_lr_map = {}
    with open(file1, 'r') as f:
        for line in f:
            contig, longread = line.strip().split()
            if contig not in contig_lr_map:
                contig_lr_map[contig] = set()
            contig_lr_map[contig].add(longread)

    # create the graph
    G = nx.Graph()

    # Read file2 and skip the first two lines
    with open(file2, 'r') as f:
        next(f)  # skip first line
        next(f)  # skip second line
        for line in f:
            contig1, contig2, weight = line.strip().split()
            if contig1 not in contig_lr_map or contig2 not in contig_lr_map:
                continue  # Skip the edge if either contig is not in contig_lr_map
            G.add_edge(contig1, contig2, weight=int(weight))

    # create a new graph for the wired contigs
    G1 = nx.Graph()

    # create a list to hold the wired nodes and the nodes to which they are wired
    wired_nodes = []

    # iterate over nodes in the original graph
    for n in G.nodes():

        # get the sorted list of neighbors based on intersection of long read IDs
        neighbors = G.neighbors(n)
        sorted_neighbors = sorted([ni for ni in neighbors if ni in contig_lr_map],
                                  key=lambda x: len(contig_lr_map[x].intersection(contig_lr_map[n])), reverse=True)

        # if the node has at least two neighbors that satisfy the minimum degree requirement
        if len(sorted_neighbors) >= 2:

            # wire the top two neighbors
            G1.add_edge(n, sorted_neighbors[0])
            G1.add_edge(n, sorted_neighbors[1])

            # append the wired nodes to the list
            wired_nodes.append((n, sorted_neighbors[0]))
            wired_nodes.append((n, sorted_neighbors[1]))

        # if the node has less than two neighbors that satisfy the minimum degree requirement
        else:
            for neighbor in G.neighbors(n):
                G1.add_edge(n, neighbor)

    # write the wired graph to a file
    with open(output_file, 'w') as f:
        for edge in sorted(G1.edges(), key=lambda x: (int(x[0]), int(x[1]))):
            contig1, contig2 = edge
            shared_lr_count = G[contig1][contig2]['weight']
            f.write(f"{contig1} {contig2} {shared_lr_count}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 Wiring.py file1 file2 output_file")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_file = sys.argv[3]

    wire_contigs(file1, file2, output_file)
