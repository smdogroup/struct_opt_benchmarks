#!/usr/bin/env python3

from py_struct_opt_benchmarks.wrapper import PlaneStressAnalysis
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import argparse
import pickle
import os

def plot_stress(picklename, savefig=False):

    # Load in pickle file
    with open(picklename, 'rb') as pklfile:
        prob_pkl = pickle.load(pklfile)

    # Get data
    nnodes = prob_pkl['nnodes']
    C = prob_pkl['C']
    conn = prob_pkl['conn']
    X = prob_pkl['X']
    dof = prob_pkl['dof']
    force = prob_pkl['force']
    r0 = prob_pkl['r0']
    density = prob_pkl['density']

    try:
        x = prob_pkl['x']
    except:
        print("\n[Warning] pkl file doesn't contain design!\n")
        x = None

    try:
        qval = prob_pkl['qval']
    except:
        print("\n[Warning] pkl file doesn't contain qval, set qval = 3.0.\n")
        qval = 3.0

    try:
        epsilon = prob_pkl['epsilon']
    except:
        print("\n[Warning] pkl file doesn't contain epsilon, set epsilon = 0.1.\n")
        epsilon = 0.1

    try:
        ks_parameter = prob_pkl['ks_parameter']
    except:
        print("\n[Warning] pkl file doesn't contain ks_parameter, set ks_parameter = 50.0.\n")
        ks_parameter = 50.0

    try:
        design_freq = prob_pkl['design_freq']
    except:
        design_freq = None

    try:
        design_stress = prob_pkl['design_stress']
    except:
        design_stress = None

    # Check if we have a solution in pickle or not
    if x is None:
        x = np.zeors(nnodes)

    # Instantiate analysis because we need to compute filtered design
    analysis = PlaneStressAnalysis(conn, dof, X, force, r0, C, density,
        qval, epsilon, ks_parameter, design_stress=design_stress,
        design_freq=design_freq, compute_comp=True, compute_mass=True,
        compute_freq=True, compute_stress=True)

    # Plot or save fig
    figname = os.path.splitext(picklename)[0]
    figname += '_stress.png'
    analysis.plot_stress(x, savefig=savefig, filename=figname)

if __name__ == '__main__':

    # Set up parser
    p = argparse.ArgumentParser()
    p.add_argument('picklename', type=str)
    p.add_argument('--savefig', action='store_true')
    args = p.parse_args()

    plot_stress(args.picklename, savefig=args.savefig)