import networkx as nx
import pandas as pd
from networkx.utils import reverse_cuthill_mckee_ordering
import time

# get the start time
st = time.time()

dict_lr_set_bench = {}

file_paf_open = open('~/CE_C_LPairsFromAsym.log', 'r')
#for line in tqdm(file_paf_open):
for line in file_paf_open:
    line = line.split('\n')[0].split(' ')
    line = [int(i) for i in line]
    if line[0] not in dict_lr_set_bench:
        dict_lr_set_bench[line[0]]  = []
        dict_lr_set_bench[line[0]].append(line[1])
    else:
        if line[1] not in dict_lr_set_bench[line[0]]:
            dict_lr_set_bench[line[0]].append(line[1])

file_w = open('~/CE_CCPairsFromAsym.log', 'w')
list_contig = list(dict_lr_set_bench.keys())
#for c1 in tqdm(list_contig):
for c1 in list_contig:
    for c2 in list_contig[list_contig.index(c1):]:
        inter = set(dict_lr_set_bench[c2]).intersection(set(dict_lr_set_bench[c1]))
        if len(inter) > 0 and c1 != c2:
            # print(c1, c2)
            #for lr in inter:
            file_w.write(f'{c1} {c2} {len(inter)}\n')
                #file_w.write(f'{c1} {c2} {len(inter)}\n')

print("Graph construction done!!!")

et = time.time()

# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')

# get the start time
st = time.time()



df = pd.read_csv('~/CLPairs_from8_CE.log', sep =' ', header=None)

df.head()


# In[3]:


map_contig = {}
for i, row in df.iterrows():
    if row[0] not in map_contig.keys():
        map_contig[row[0]] = [row[1]]
    else:
        map_contig[row[0]].append(row[1])


# In[7]:


df_g = pd.read_csv('~/CCPairs_CE.log', sep =' ', header=None)

G = nx.from_pandas_edgelist(df_g, 0, 1, [2])
df_g.head()
G1 = nx.Graph()
no_edge = 2


count_lr = []
for n in G.nodes():
    if G.degree[n] > no_edge:
#         print(n)
        for x in range(no_edge):
            l = [ni for ni in G[n]]
#             print(l)
            mx = 0
            mi = 0
            mn = 0
            for i in l:
                x = len(set(map_contig[i]).intersection(set(map_contig[n])))
                count_lr.append(x)
#                 print(n, i, x)
                if mx < x:
                    mx = x
                    mi = i
                    mn = n
            G1.add_edge(mn, mi,weight=x)
    else:
        for ni in G[n]:
            x = len(set(map_contig[ni]).intersection(set(map_contig[n])))
            G1.add_edge(n, ni,weight=x)

nx.write_weighted_edgelist(G1, "LRInter_New_CE.log", delimiter=' ')

et = time.time()

# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')


df = pd.read_csv('CE_CCPairsFromAsym.log', sep=' ', header=None)
df.head()

# In[18]:


G = nx.from_pandas_edgelist(df, 0, 1, [2])


# In[19]:


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


# In[20]:


G = nx.from_pandas_edgelist(df, 0, 1, [2])
rcm, mapping = reverse_cuthill_mckee(G)

# In[21]:


file_open = open('OrderedRCM_CE.log', 'w')
for i, j in rcm.edges():
    file_open.write(f'{i} {j}\n')

# In[22]:


file_open = open('OrderedRCM_trueID_CE.log', 'w')
for i in mapping.keys():
    file_open.write(f'{i} {mapping[i]}\n')

print("RCM Creation done!!!")

et = time.time()

# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')