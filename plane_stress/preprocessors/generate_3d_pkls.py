#!/usr/bin/env python3

"""
This script can be used to generate a batch of pkl files
"""

import argparse
from plane_stress.preprocessors import preproc3d

def generate_pkls(ns, ARs, loaded_thicks, nr0, plot_mesh):

    outdir = None
    force_magnitude = 25.0  # total force applied to structure
    forced_portion = 0.3  # portion of the edge of which load is applied to
    MBB_bc_portion = 0.3
    ly = 1.0  # for this project we always set height of domain to be 1.0
    density = 2700.0  # this is only used for frequency analysis
    E = 70e3  # Young's modulus
    nu = 0.3  # Poisson's ratio
    ratio1 = 0.4
    ratio2 = 0.4
    hole_r = 0.25


    for n in ns:
        for AR in ARs:
            for loaded_thick in loaded_thicks:
                for prob in ['cantilever', 'michell', 'MBB', 'lbracket']:
                    for meshtype in ['structured', 'unstructured']:
                            preproc3d(n, AR, prob, ratio1, ratio2, hole_r, meshtype, nr0,
                                    outdir, plot_mesh, forced_portion, force_magnitude,
                                    MBB_bc_portion, loaded_thick, ly, density, E, nu)

if __name__ == '__main__':

    p = argparse.ArgumentParser()
    p.add_argument('--n', nargs='*', type=int, default=None,
        help='number of elements on short edge')
    p.add_argument('--AR', nargs='*', type=float, default=None,
        help='domain aspect ratio, AR = length / height')
    p.add_argument('--loaded_thickness', nargs='*', type=float, default=None)
    p.add_argument('--nr0', type=int, default=6)
    p.add_argument('--plot_mesh', action='store_true')
    args = p.parse_args()

    generate_pkls(args.n, args.AR, args.loaded_thickness, args.nr0, args.plot_mesh)
