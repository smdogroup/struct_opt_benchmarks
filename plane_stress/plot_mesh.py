import matplotlib.pyplot as plt
from matplotlib.patches import Arrow, Polygon
import numpy as np
import argparse
import pickle

# Set up parser
p = argparse.ArgumentParser('Takes in a problem object file in plk format and plot mesh')
p.add_argument('filename', metavar='cantilever.pkl', type=str)
args = p.parse_args()

# Load in pickle file
with open(args.filename, 'rb') as pklfile:
    prob_pkl = pickle.load(pklfile)

# Get data
prob_name = prob_pkl['prob_name']
nelems = prob_pkl['nelems']
nnodes = prob_pkl['nnodes']
nvars = prob_pkl['nvars']
C = prob_pkl['C']
conn = prob_pkl['conn']
X = prob_pkl['X']
vars = prob_pkl['vars']
force = prob_pkl['force']
r0 = prob_pkl['r0']

# Loop over all elememts to plot mesh edges
for i in range(nelems):
    order = [0,1,3,2,0]
    for k in range(4):
        pt1 = X[conn[i, order[k  ]]]
        pt2 = X[conn[i, order[k+1]]]

        plt.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]],'k', lw=0.5)


# Compute the size of shapes
xmax, ymax = np.amax(X, axis=0)
xmin, ymin = np.amin(X, axis=0)
domain_size = max(ymax-ymin, xmax-xmin)
shape_size = domain_size * 0.15
f_size = np.max(np.abs(force))

# Plot forces
var_indices = np.nonzero(force)[0]  # because force is 1D array, we only need 1st index
for var_index in var_indices:
    node_index, i = np.where(vars==var_index)
    node_index = node_index[0]
    i = i[0]
    x = X[node_index, 0]
    y = X[node_index, 1]
    h = force[var_index]/f_size*shape_size
    arrow = Arrow(x-h+i*h, y-i*h, h-i*h, i*h, edgecolor='red',
        width=shape_size*0.5, fill=None, lw=1.0)
    plt.gca().add_patch(arrow)

# Plot boundary condition
nodes = np.where(vars == -1)
nodes = np.array(nodes).transpose()
for node in nodes:
    node_index = node[0]
    direction = node[1]
    x = X[node_index, 0]
    y = X[node_index, 1]
    h = shape_size*0.2
    if direction == 0:
        points = [[x, y], [x-h, y+h/3**0.5], [x-h, y-h/3**0.5]]
        triangle = Polygon(points, fill=None, edgecolor='blue', lw=1.0)
    else:
        points = [[x, y], [x-h/3**0.5, y-h], [x+h/3**0.5, y-h]]
        triangle = Polygon(points, fill=None, edgecolor='green', lw=1.0)
    plt.gca().add_patch(triangle)

plt.axis('equal')
plt.title(prob_name)
plt.show()
