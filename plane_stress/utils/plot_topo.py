#!/usr/bin/env python3

from plane_stress.wrapper import PlaneStressAnalysis
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import argparse
import pickle
import os

def plot_topo(picklename, savefig=False, filename=None, paperstyle=False):

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
    analysis.plot_topology(x, savefig=savefig, filename=filename, paperstyle=paperstyle)

if __name__ == '__main__':

    # Set up parser
    p = argparse.ArgumentParser()
    p.add_argument('--result_folder', type=str, default='results')
    p.add_argument('--opt_problem', nargs='*', type=str, default=None, choices=[
        'comp_min_mass_constr', 'comp_min_massfreq_constr',
        'comp_min_massstress_constr', 'comp_min_massfreqstress_constr',
        'stress_min_mass_constr', 'mass_min_stress_constr'])
    p.add_argument('--optimizer', nargs='*', type=str, default=None, choices=[
        'ParOpt', 'ParOptAdapt', 'ParOptFilter', 'ParOptFilterSoc', 'IPOPT', 'SNOPT'])
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
    args = p.parse_args()

    # Value checks
    if args.opt_problem is None:
        raise ValueError('--opt_problem must be specified!')

    if args.optimizer is None:
        raise ValueError('--optimizer must be specified!')

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
    for opt_problem in args.opt_problem:
        for optimizer in args.optimizer:
            for n in args.n:
                for AR in args.AR:
                    for meshtype in args.meshtype:
                        for domain in args.domain:
                            for hole in args.hole:
                                if hole == 'hole' and meshtype == 'unstructured':
                                    extra = '-r1-0.4-r2-0.4'
                                else:
                                    extra = ''

                                foldername = '{:s}-{:s}-{:s}-n{:d}-AR{:.1f}{:s}-distriforce'.format(
                                    opt_problem, meshtype, domain, n, AR, extra)

                                try:
                                    files = os.listdir('{:s}/{:s}'.format(args.result_folder, foldername))
                                except:
                                    print('[Warning] problem folder {}/{} not found!'.format(
                                        args.result_folder, foldername))
                                else:
                                    for f in files:
                                        try:
                                            optimizer = f.split('-')[1]
                                        except:
                                            continue
                                        else:
                                            if '.pkl' in f and optimizer in args.optimizer:
                                                if args.paperstyle:
                                                    filename = '{:s}-{:s}_topo.pdf'.format(optimizer, foldername)
                                                else:
                                                    filename = '{:s}-{:s}_topo.png'.format(optimizer, foldername)

                                                pklfullpath = '{:s}/{:s}/{:s}'.format(args.result_folder, foldername, f)
                                                plot_topo(pklfullpath, savefig=args.savefig, filename=filename, paperstyle=args.paperstyle)