import networkx as nx
import pandas as pd
from networkx.utils import reverse_cuthill_mckee_ordering



df = pd.read_csv('~/HM7_CCPairsFromAsym.log', sep=' ', header=None)
df.head()
G = nx.from_pandas_edgelist(df, 0, 1, [2])

def reverse_cuthill_mckee(G):
    nx.set_node_attributes(G, False, "visited")
    nx.set_node_attributes(G, dict(G.degree), "degree")

    reorder = []

    for component in list(nx.connected_components(G)):

        frontier = []
        startVertex = sorted(component, key=lambda x: G.nodes[x]["degree"], reverse=False)[0]
        frontier.append(startVertex)
        G.nodes[startVertex]["visited"] = True

        while (len(frontier) > 0):
            n = frontier.pop(0)
            reorder.append(n)

            for child in sorted(G.neighbors(n), key=lambda x: G.nodes[x]["degree"], reverse=False):
                if (G.nodes[child]["visited"] == False):
                    G.nodes[child]["visited"] = True
                    frontier.append(child)

    mapping = {}

    for i in range(len(reorder)):
        mapping[reorder[i]] = i

    G = nx.relabel_nodes(G, mapping)

    return G, mapping




G = nx.from_pandas_edgelist(df, 0, 1, [2])
rcm, mapping = reverse_cuthill_mckee(G)





file_open = open('~/OrderedRCM_HM7.log', 'w')
for i, j in rcm.edges():
    file_open.write(f'{i} {j}\n')



file_open = open('~/OrderedRCM_trueID_HM7.log', 'w')
for i in mapping.keys():
    file_open.write(f'{i} {mapping[i]}\n')
    
print("RCM Creation done!!!")
