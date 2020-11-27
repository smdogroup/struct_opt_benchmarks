#!/usr/bin/env python3

import pickle
import os
import argparse

def pkltovtk(pklfile):

    # Get name
    base = os.path.basename(pklfile)
    name = os.path.splitext(base)[0]
    vtkfile = '{:s}.vtk'.format(name)

    # Load pickle file
    with open(pklfile, 'rb') as fh:
        pkl_dict = pickle.load(fh)

    # Get pkl data
    nnodes = pkl_dict['nnodes']
    nelems = pkl_dict['nelems']
    X = pkl_dict['X']
    conn = pkl_dict['conn']

    try:
        x = pkl_dict['x']
    except:
        has_design = False
        print('\n[Warning] No design contained in pickle file!')
    else:
        has_design = True

    # Create a empty vtk file and write headers
    with open(vtkfile, 'w') as fh:
        fh.write('# vtk DataFile Version 3.0\n')
        fh.write('my example\n')
        fh.write('ASCII\n')
        fh.write('DATASET UNSTRUCTURED_GRID\n')

        # Write nodal points
        fh.write('POINTS {:d} float\n'.format(nnodes))
        for i in range(nnodes):
            fh.write('{:f} {:f} {:f}\n'.format(X[i, 0], X[i, 1], X[i, 2]))

        # Write connectivity
        fh.write('CELLS {:d} {:d}\n'.format(nelems, 9*nelems))
        for i in range(nelems):
            fh.write('8 {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d}\n'.format(
                conn[i, 0], conn[i, 1], conn[i, 3], conn[i, 2],
                conn[i, 4], conn[i, 5], conn[i, 7], conn[i, 6]))

        # Write cell type
        fh.write('CELL_TYPES {:d}\n'.format(nelems))
        for i in range(nelems):
            fh.write('12\n')

        # Write nodal density
        if has_design:
            fh.write('POINT_DATA {:d}\n'.format(nnodes))
            fh.write('SCALARS density float 1\n')
            fh.write('LOOKUP_TABLE default\n')
            for i in range(nnodes):
                fh.write('{:f}\n'.format(x[i]))


if __name__ == '__main__':

    # Set up argparser
    p = argparse.ArgumentParser()
    p.add_argument('pklfile', type=str)
    args = p.parse_args()

    # Call function
    pkltovtk(args.pklfile)



