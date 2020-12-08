#!/usr/bin/env python3

"""
This script can be used to generate a batch of pkl files
"""

import argparse
from py_struct_opt_benchmarks.preprocessors import preproc

def generate_pkls(ns, ARs, plot_mesh):

    ratio1 = 0.4
    ratio2 = 0.4
    hole_radius = 0.25
    nr0 = 32
    outdir = None
    use_concentrated_force = False

    force_magnitude = 25.0  # total force applied to structure
    forced_portion = 0.2  # portion of the edge of which load is applied to
    MBB_bc_portion = 0.1
    ly = 1.0  # for this project we always set height of domain to be 1.0
    density = 2700.0  # this is only used for frequency analysis
    E = 70e3  # Young's modulus
    nu = 0.3  # Poisson's ratio

    for n in ns:
        for AR in ARs:
            for prob in ['cantilever', 'michell', 'MBB', 'lbracket']:
                for meshtype in ['structured', 'unstructured']:
                    for use_hole in [False, True]:
                        if meshtype == 'structured' and use_hole == True:
                            continue
                        else:
                            preproc(n, AR, prob, meshtype, ratio1, ratio2,
                            use_concentrated_force, use_hole, hole_radius,
                            nr0, outdir, plot_mesh, force_magnitude,
                            forced_portion, MBB_bc_portion, ly, density, E, nu)

if __name__ == '__main__':

    p = argparse.ArgumentParser()
    p.add_argument('--n', nargs='*', type=int, default=None,
        help='number of elements on short edge')
    p.add_argument('--AR', nargs='*', type=float, default=None,
        help='domain aspect ratio, AR = length / height')
    p.add_argument('--plot_mesh', action='store_true')
    args = p.parse_args()

    generate_pkls(args.n, args.AR, args.plot_mesh)
