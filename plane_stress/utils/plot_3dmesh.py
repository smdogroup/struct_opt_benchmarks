#!/usr/bin/env python3

'''
This script can plot 3D meshes with 8-node solid elements
'''

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import Arrow, Polygon, FancyArrowPatch
import numpy as np
import os

class myArrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

def plot_3dmesh(conn, X, dof, force):

    # Get number of elements
    nelems = conn.shape[0]

    # Create 3d axes
    fig = plt.figure()
    ax = Axes3D(fig)

    # Loop over all elements
    for e in range(nelems):

        # Plot surface 1
        x = [X[conn[e, n], 0] for n in [0,1,3,2]]
        y = [X[conn[e, n], 1] for n in [0,1,3,2]]
        z = [X[conn[e, n], 2] for n in [0,1,3,2]]
        verts = [list(zip(x, y, z))]
        ax.add_collection3d(Poly3DCollection(verts, edgecolor='black',
            alpha=0, linewidth=0.5))

        # Plot surface 2
        x = [X[conn[e, n], 0] for n in [4,5,7,6]]
        y = [X[conn[e, n], 1] for n in [4,5,7,6]]
        z = [X[conn[e, n], 2] for n in [4,5,7,6]]
        verts = [list(zip(x, y, z))]
        ax.add_collection3d(Poly3DCollection(verts, edgecolor='black',
            alpha=0, linewidth=0.5))

        # Plot surface 3
        x = [X[conn[e, n], 0] for n in [0,2,6,4]]
        y = [X[conn[e, n], 1] for n in [0,2,6,4]]
        z = [X[conn[e, n], 2] for n in [0,2,6,4]]
        verts = [list(zip(x, y, z))]
        ax.add_collection3d(Poly3DCollection(verts, edgecolor='black',
            alpha=0, linewidth=0.5))

        # Plot surface 4
        x = [X[conn[e, n], 0] for n in [1,3,7,5]]
        y = [X[conn[e, n], 1] for n in [1,3,7,5]]
        z = [X[conn[e, n], 2] for n in [1,3,7,5]]
        verts = [list(zip(x, y, z))]
        ax.add_collection3d(Poly3DCollection(verts, edgecolor='black',
            alpha=0, linewidth=0.5))

    # Compute sizes of figure
    xmax, ymax, zmax = np.amax(X, axis=0)
    xmin, ymin, zmin = np.amin(X, axis=0)
    lx = xmax - xmin
    ly = ymax - ymin
    lz = zmax - zmin
    domain_size = np.max([lx, ly, lz])
    shape_size = domain_size * 0.10
    f_size = np.max(np.abs(force))

    # set limits for axes
    white_space = 0.2
    ax.set_xlim3d(xmin - white_space*lx, xmax + white_space*lx)
    ax.set_ylim3d(ymin - white_space*ly, ymax + white_space*ly)
    ax.set_zlim3d(zmin - white_space*lz, zmax + white_space*lz)

    # Plot forces
    dof_indices = np.nonzero(force)[0]  # because force is 1D array, we only need 1st index
    for dof_index in dof_indices:
        node_index, direction = np.where(dof==dof_index)
        node_index = node_index[0]
        direction = direction[0]
        x = X[node_index, 0]
        y = X[node_index, 1]
        z = X[node_index, 2]
        h = force[dof_index]/f_size*shape_size

        # x-directional force
        if direction == 0:
            arrow = myArrow3D([x-h,x], [y,y], [z,z], mutation_scale=10, lw=1.5,
                arrowstyle='simple', color='red')

        # y-directional force
        elif direction == 1:
            arrow = myArrow3D([x, x], [y-h,y], [z,z], mutation_scale=10, lw=1.5,
                arrowstyle='simple', color='red')

        # z-directional force
        elif direction == 2:
            arrow = myArrow3D([x, x], [y,y], [z-h,z], mutation_scale=10, lw=1.5,
                arrowstyle='simple', color='red')

        ax.add_artist(arrow)

    # Plot boundary condition
    nodes = np.where(dof == -1)
    nodes = np.array(nodes).transpose()
    for node in nodes:
        node_index = node[0]
        direction = node[1]
        x = X[node_index, 0]
        y = X[node_index, 1]
        z = X[node_index, 2]
        h = shape_size*0.2

        # x-directional force
        if direction == 0:
            arrow = myArrow3D([x-h,x], [y,y], [z,z], mutation_scale=10, lw=1.5,
                arrowstyle='simple', color='blue')

        # y-directional force
        elif direction == 1:
            arrow = myArrow3D([x, x], [y-h,y], [z,z], mutation_scale=10, lw=1.5,
                arrowstyle='simple', color='green')

        # z-directional force
        elif direction == 2:
            arrow = myArrow3D([x, x], [y,y], [z-h,z], mutation_scale=10, lw=1.5,
                arrowstyle='simple', color='yellow')

        ax.add_artist(arrow)


    plt.show()