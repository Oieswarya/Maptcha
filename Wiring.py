import networkx as nx
import pandas as pd
import numpy as np

import collections as cs
import time

# get the start time
st = time.time()



df = pd.read_csv('~\CLPairs_from8_CE.log', sep =' ', header=None)
df.head()




map_contig = {}
for i, row in df.iterrows():
    if row[0] not in map_contig.keys():
        map_contig[row[0]] = [row[1]]
    else:
        map_contig[row[0]].append(row[1])




df_g = pd.read_csv('~\CCPairs_CE.log', sep =' ', header=None)

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
