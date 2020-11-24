#!/usr/bin/env python3

'''
This script can plot 2D meshes with quadrilateral elements
'''

import matplotlib.pyplot as plt
from matplotlib.patches import Arrow, Polygon
import numpy as np
import argparse
import pickle
import os

def plot_mesh(prob_pkl, savefig, outdir, paperstyle=False):

    # Get data
    prob_name = prob_pkl['prob_name']
    nelems = prob_pkl['nelems']
    conn = prob_pkl['conn']
    X = prob_pkl['X']
    dof = prob_pkl['dof']
    force = prob_pkl['force']

    ptx = X[:,0]
    pty = X[:,1]

    length = max(ptx) - min(ptx)
    height = max(pty) - min(pty)

    fig_height = 4.8
    fig_length = 4.8 * length / height

    # Set up plotting environment
    if paperstyle:
        fig, ax = plt.subplots(figsize=(fig_length, fig_height), constrained_layout=True)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.axis('off')
    else:
        fig, ax = plt.subplots()

    # Loop over all elememts to plot mesh edges
    for i in range(nelems):
        x = [X[conn[i, j], 0] for j in [0,1,3,2]]
        y = [X[conn[i, j], 1] for j in [0,1,3,2]]
        ax.fill(x, y, edgecolor='black', fill=False, lw=0.5)

    # Compute the size of shapes
    xmax, ymax = np.amax(X, axis=0)
    xmin, ymin = np.amin(X, axis=0)
    domain_size = max(ymax-ymin, xmax-xmin)
    shape_size = domain_size * 0.15
    f_size = np.max(np.abs(force))

    # Plot forces
    dof_indices = np.nonzero(force)[0]  # because force is 1D array, we only need 1st index
    for dof_index in dof_indices:
        node_index, i = np.where(dof==dof_index)
        node_index = node_index[0]
        i = i[0]
        x = X[node_index, 0]
        y = X[node_index, 1]
        h = force[dof_index]/f_size*shape_size
        arrow = Arrow(x-h+i*h, y-i*h, h-i*h, i*h, edgecolor='red',
            width=shape_size*0.4, fill=None, lw=1.0)
        ax.add_patch(arrow)

    # Plot boundary condition
    nodes = np.where(dof == -1)
    nodes = np.array(nodes).transpose()
    for node in nodes:
        node_index = node[0]
        direction = node[1]
        x = X[node_index, 0]
        y = X[node_index, 1]
        h = shape_size*0.15
        if direction == 0:
            points = [[x, y], [x-h, y+h/3**0.5], [x-h, y-h/3**0.5]]
            triangle = Polygon(points, fill=None, edgecolor='blue', lw=1.0)
        else:
            points = [[x, y], [x-h/3**0.5, y-h], [x+h/3**0.5, y-h]]
            triangle = Polygon(points, fill=None, edgecolor='green', lw=1.0)
        ax.add_patch(triangle)

    # Set the aspect ratio equal
    ax.set_aspect('equal')

    # Set title
    if not paperstyle:
        ax.set_title('{:s}\nloaded nodes:  {:d} max force: {:.2f}'.format(
            prob_name, np.count_nonzero(force), np.max(np.abs(force))))

    if savefig:
        if outdir is not None:
        # Check if specified output directory exists
            if os.path.isdir(outdir):
                # outdir exists
                pass
            # Otherwise, try to create the folder
            else:
                try:
                    os.mkdir(outdir)
                except:
                    print("\n[Warning] Cannot create directory {:s}!\n".format(outdir))
                    outdir = None

        if outdir is not None:
            if paperstyle:
                plt.savefig(outdir+'/'+prob_name+'_mesh.pdf')
            else:
                plt.savefig(outdir+'/'+prob_name+'_mesh.png')
        else:
            if paperstyle:
                plt.savefig(prob_name+'_mesh.pdf')
            else:
                plt.savefig(prob_name+'_mesh.png')

        plt.close()
    else:
        plt.show()

if __name__ == '__main__':

    # Set up parser
    p = argparse.ArgumentParser('Takes in a problem object file in plk format and plot mesh')
    p.add_argument('--pkl_folder', type=str, default='pkls')
    p.add_argument('--n', nargs='*', type=int, default=None, choices=[48, 64, 80])
    p.add_argument('--AR', nargs='*', type=float, default=None, choices=[1.0, 2.0, 3.0])
    p.add_argument('--meshtype', nargs='*', type=str, default=None, choices=[
        'structured', 'unstructured'])
    p.add_argument('--domain', nargs='*', type=str, default=None, choices=[
        'cantilever', 'michell', 'MBB', 'lbracket'])
    p.add_argument('--hole', nargs='*', type=str, default=None, choices=[
        'hole', 'nohole'])
    p.add_argument('--savefig', action='store_true')
    p.add_argument('--paperstyle', action='store_true')
    p.add_argument('--outdir', type=str)
    args = p.parse_args()

    # Value checks
    if args.n is None:
        raise ValueError('--n must be specified!')

    if args.AR is None:
        raise ValueError('--AR must be specified!')

    if args.meshtype is None:
        raise ValueError('--meshtype must be specified!')

    if args.domain is None:
        raise ValueError('--domain must be specified!')

    if args.hole is None:
        raise ValueError('--hole must be specified!')

    # Plot
    for n in args.n:
        for AR in args.AR:
            for meshtype in args.meshtype:
                for domain in args.domain:
                    for hole in args.hole:
                        if domain == 'lbracket':
                            extral = '-r1-0.4-r2-0.4'
                        else:
                            extral = ''
                        if hole == 'hole' and meshtype == 'unstructured':
                            extrah = '-hole0.2'
                        else:
                            extrah = ''

                        pklname = '{:s}-{:s}-n{:d}-AR{:.1f}{:s}-distriforce{:s}.pkl'.format(
                            meshtype, domain, n, AR, extral, extrah)
                        pkl_path = '{:s}/{:s}'.format(args.pkl_folder, pklname)

                        # Load pickle
                        with open(pkl_path, 'rb') as pklfh:
                            pkl_dict = pickle.load(pklfh)

                        try:
                            plot_mesh(pkl_dict, args.savefig, args.outdir, paperstyle=args.paperstyle)
                        except Exception as e:
                            print(e)
                            print('[Error] mesh pickle file {:s} not found!'.format(pkl_path))