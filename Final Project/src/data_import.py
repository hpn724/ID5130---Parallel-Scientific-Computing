import numpy as np
import networkx as nx
from networkx.algorithms.community.modularity_max import greedy_modularity_communities
import pandas as pd
import csv
import pickle

#import the graph
g = nx.karate_club_graph()

communities = greedy_modularity_communities(g)

A = nx.to_numpy_array(g)

np.savetxt('Zacharys_Karate_club.csv',A,delimiter=",",fmt='%.0f')


with open('Communities.csv','w',newline='') as csvfile:
    writer = csv.writer(csvfile)
    for community in communities:
        writer.writerow((community))
        
with open("graph_variables.pk1",'wb') as f:
    pickle.dump((g,communities,A),f)