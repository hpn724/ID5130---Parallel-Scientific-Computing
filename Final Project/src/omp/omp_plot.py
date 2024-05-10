import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from IPython.display import HTML
import pickle

def draw_kkl(nx_G, label_map, node_color, pos=None, **kwargs):
    fig, ax = plt.subplots(figsize=(10,10))
    if pos is None:
        pos = nx.spring_layout(nx_G, k=5/np.sqrt(nx_G.number_of_nodes()))

    nx.draw(
        nx_G, pos, with_labels=label_map is not None, 
        labels=label_map, 
        node_color=node_color, 
        ax=ax, **kwargs)
    

with open('graph_variables.pk1','rb') as f:
    g, communities, A = pickle.load(f)
    
colors = np.zeros(g.number_of_nodes())
for i, com in enumerate(communities):
    colors[list(com)] = i

n_classes = np.unique(colors).shape[0]

labels = np.eye(n_classes)[colors.astype(int)]

club_labels = nx.get_node_attributes(g,'club')




###########################################
"""            Initial graph            """
###########################################

fig, ax = plt.subplots(figsize=(10,10))
pos = nx.spring_layout(g, k=5/np.sqrt(g.number_of_nodes()))
kwargs = {"cmap": 'gist_rainbow', "edge_color":'gray'}
nx.draw(
    g, pos, with_labels=False, 
    node_color=colors, 
    ax=ax, **kwargs)
plt.savefig('Solution/omp/karate_club_graph_initial.png', bbox_inches='tight')



rows = 34
cols = 2


threads_used = [1,2,4,6,8,10,12]

for i in threads_used:
    file_name_1 = f"Solution/omp/GCN_embed_omp_{i}.csv"
    GCN_omp_data = np.loadtxt(open(file_name_1,'rb'),delimiter=",")

    file_name_2 = f"Solution/omp/GCN_omp_{i}_eval_params.csv"
    GCN_omp_eval = pd.read_csv(file_name_2,index_col=0,header=0)


    fig, ax = plt.subplots(figsize=(10,10))
    pos = {k: GCN_omp_data[k,:] for k in range(GCN_omp_data.shape[0])}
    kwargs = {"cmap": 'gist_rainbow', "edge_color":'gray'}
    nx.draw(
        g, pos, with_labels=False, 
        node_color=colors, 
        ax=ax, **kwargs)
    plt.savefig(f'Solution/omp/karate_club_graph_final_omp_{i}.png', bbox_inches='tight')

    omp_tt_fig = plt.figure()
    plt.plot(list(GCN_omp_eval['Train loss']), label = 'Train loss')
    plt.plot(list(GCN_omp_eval['Test loss']), label = 'Test loss')
    plt.legend()
    plt.grid()
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.savefig(f'Solution/omp/omp_{i}_tt_loss_eval.png')

plt.close('all')
    
for i in threads_used:
    file_name_2 = f"Solution/omp/GCN_omp_{i}_eval_params.csv"
    GCN_omp_eval = pd.read_csv(file_name_2,index_col=0,header=0)
    plt.plot(list(GCN_omp_eval['Accuracy']),label=f'Threads used : {i}')
    plt.grid()
    plt.ylabel('accuracy')
    plt.xlabel('Epoch')
    plt.legend()
    plt.savefig('Solution/omp/omp_accuracy_eval.png')

omp_time_taken = plt.figure()
file_name_3 = 'Solution/omp/GCN_omp_time_taken_2.csv'
time_taken_omp = pd.read_csv(file_name_3,header=0)
plt.plot(list(time_taken_omp['Number of threads']),list(time_taken_omp['Time taken']))
plt.grid()
plt.ylabel('Time taken (in seconds)')
plt.xlabel('Number of threads')
plt.savefig('Solution/omp/omp_time_taken.png')
