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
plt.savefig('Solution/Serial/karate_club_graph_initial.png', bbox_inches='tight')



depth = 7000
rows = 34
cols = 2


GCN_serial_data = np.fromfile("Solution/Serial/GCN_embed_serial.bin",dtype=np.float64)
GCN_serial_embeds = GCN_serial_data.reshape((depth,rows,cols))

GCN_serial_eval = pd.read_csv("Solution/Serial/GCN_serial_eval_params.csv",index_col=0,header=0)


fig, ax = plt.subplots(figsize=(10,10))
pos = {i: GCN_serial_embeds[-1][i,:] for i in range(GCN_serial_embeds[-1].shape[0])}
kwargs = {"cmap": 'gist_rainbow', "edge_color":'gray'}
nx.draw(
    g, pos, with_labels=False, 
    node_color=colors, 
    ax=ax, **kwargs)
plt.savefig('Solution/Serial/karate_club_graph_final.png', bbox_inches='tight')

serial_tt_fig = plt.figure()
plt.plot(list(GCN_serial_eval['Train loss']), label = 'Train loss')
plt.plot(list(GCN_serial_eval['Test loss']), label = 'Test loss')
plt.legend()
plt.grid()
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.savefig('Solution/Serial/serial_tt_loss_eval.png')

serial_acc_fig = plt.figure()
plt.plot(list(GCN_serial_eval['Accuracy']))
plt.grid()
plt.ylabel('accuracy')
plt.xlabel('Epoch')
plt.savefig('Solution/Serial/serial_accuracy_eval.png')

## To create a video of the clustering

"""
N = 100 
snapshots = np.linspace(0, len(GCN_serial_embeds)-1, N).astype(int)
fig, ax = plt.subplots(figsize=(10, 10))
kwargs = {'cmap': 'gist_rainbow', 'edge_color': 'gray'}  # Adjust visualization settings

def update(idx):
    ax.clear()
    embed = GCN_serial_embeds[snapshots[idx]]
    pos = {i: embed[i,:] for i in range(embed.shape[0])}
    nx.draw(g, pos, node_color=colors, ax=ax, **kwargs)

anim = animation.FuncAnimation(fig, update, frames=snapshots.shape[0], interval=50, repeat=False)  # Adjust interval
anim.save('Solution/Serial/Clustering.mp4', dpi=300)

# Close figures
plt.close('all')

"""



