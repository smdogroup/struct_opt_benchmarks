#!/usr/bin/env python3

"""
This script generates a 2D plane stress finite element analysis
problem, including geometry, mesh, boundary condition and load,
and stores the problem object to a python pickle file.

The structure of the problem object is as follows:

prob_pkl: dict
    |---prob_name:    problem name,            str
    |---nelems:       number of elements,      int
    |---nnodes:       number of mesh nodes,    int
    |---ndof:         number of nodal dof,     int
    |---C:            constitutive matrix,     3*3 ndarray
    |---conn:         connectivity,            nelems*4 ndarray
    |---X:            nodal position,          nnodes*2 ndarray
    |---dof:          nodal degree of freedom, nnodes*2 ndarray
    |---force:        nodal forces,            ndof*1 ndarray
    |---r0:           filter radius,           float
    |---density:      material density,        float

domain:
rectangle
         ______________________________
        |                              |
 ny = n |                              |
 ly=1.0 |                              |
        |______________________________|
            nx = AR*n, lx = AR*ly

L-bracket

        nx1 = nx*ratio1
        lx1 = lx*ratio1
        fixed boundary
        --------------
        |  ratio1    |
        |            |ny2
        |            |
        |            |ly2           ^ force
 ny = n |            |              |
 ly=1.0 |            |______________|
        |                           |
        |                           |ratio2
        |                           |ny1 = ny*ratio2
        |                           |ly1 = ly*ratio2
        |___________________________|
            nx = AR*n, lx = AR*ly
"""

import numpy as np
import pickle
import argparse
import os

def preproc(n, AR, prob, meshtype, ratio1, ratio2,
    use_concentrated_force, use_hole, hole_radius,
    nr0, outdir, plot_mesh, force_magnitude,
    forced_portion, MBB_bc_portion, ly, density, E, nu):

    # Check input values
    if hole_radius >= 0.5:
        raise ValueError("Hole radius must be smaller than 0.5!")
    if n <= nr0:
        print('\n[Warning] Element size is larger than filter radius!\n')

    # nnodes, nelems and prob_name
    if prob == 'lbracket':

        ny = n
        nx = round(n*AR)
        lx = ly*AR
        ny1 = round(ny*ratio2)
        nx1 = round(nx*ratio1)
        ly1 = ly*ratio2
        lx1 = lx*ratio1
        ny2 = ny - ny1
        nx2 = nx - nx1
        ly2 = ly - ly1
        lx2 = lx - lx1

        prob_name = '{:s}-{:s}-n{:d}-AR{:.1f}-r1-{:.1f}-r2-{:.1f}'.format(
            meshtype, prob, n, AR, ratio1, ratio2)

        if meshtype == 'structured':
            nelems = nx1*ny1 + nx1*ny2 + nx2*ny1
            nnodes = (nx1+1)*(ny1+1)+(nx1+1)*(ny2+1)+(nx2+1)*(ny1+1)-(nx1+1)-(ny1+1)

    else:

        lx = ly*AR
        ny = n
        nx = round(ny*AR)
        prob_name = '{:s}-{:s}-n{:d}-AR{:.1f}'.format(
            meshtype, prob, n, AR)

        if meshtype == 'structured':
            nelems = nx*ny
            nnodes = (nx+1)*(ny+1)

    # Update prob_name
    if use_concentrated_force:
        prob_name += '-ccforce'
    else:
        prob_name += '-distriforce'

    if use_hole:
        if meshtype == 'unstructured':
            prob_name += '-hole{:.1f}'.format(hole_radius)

    # r0
    r0 = ly / nr0

    # C
    C = np.zeros((3, 3))
    C[0, 0] = E/(1.0 - nu**2)
    C[0, 1] = nu*E/(1.0 - nu**2)
    C[1, 0] = C[0, 1]
    C[1, 1] = C[0, 0]
    C[2, 2] = 0.5*E/(1.0 + nu)

    # We use tmr to generate unstructured mesh
    if meshtype == 'unstructured':

        from plane_stress.preprocessors.unstructured_utils import create_mesh

        nelems, nnodes, ndof, conn, X, dof, force = create_mesh(
            n, AR, prob, ratio1, ratio2,
            forced_portion, MBB_bc_portion, force_magnitude,
            use_concentrated_force, use_hole, hole_radius)

    else:
        # ndof, conn, X, dof
        conn = np.zeros((nelems, 4), dtype=np.intc)
        dof = -np.ones((nnodes, 2), dtype=np.intc)
        X = np.zeros((nnodes, 2))
        ndof = 0

        if prob == 'lbracket':
            nodei = 0
            for j in range(ny+1):
                for i in range(nx+1):
                    if i <= nx1 or j <= ny1:
                        if i <= nx1:
                            X[nodei, 0] = lx1*i/nx1
                        else:
                            X[nodei, 0] = lx1 + (i-nx1)*lx2/nx2

                        if j <= ny1:
                            X[nodei, 1] = ly1*j/ny1
                        else:
                            X[nodei, 1] = ly1 + (j-ny1)*ly2/ny2

                        if j < ny:
                            dof[nodei, 0] = ndof
                            ndof += 1
                            dof[nodei, 1] = ndof
                            ndof += 1

                        nodei += 1

            elemi = 0
            for j in range(ny):
                for i in range(nx):
                    if i < nx1 or j < ny1:
                        if j <= ny1:
                            conn[elemi, 0] = i + (nx+1)*j
                            conn[elemi, 1] = i+1 + (nx+1)*j
                            conn[elemi, 2] = i + (nx+1)*(j+1)
                            conn[elemi, 3] = i+1 + (nx+1)*(j+1)

                        else:
                            conn[elemi, 0] = i + (nx1+1)*(j-ny1-1) + (nx+1)*(ny1+1)
                            conn[elemi, 1] = i+1 + (nx1+1)*(j-ny1-1) + (nx+1)*(ny1+1)
                            conn[elemi, 2] = i + (nx1+1)*(j-ny1) + (nx+1)*(ny1+1)
                            conn[elemi, 3] = i+1 + (nx1+1)*(j-ny1) + (nx+1)*(ny1+1)

                        elemi += 1
        else:
            for j in range(ny):
                for i in range(nx):
                    conn[i + j*nx, 0] = i + (nx+1)*j
                    conn[i + j*nx, 1] = i+1 + (nx+1)*j
                    conn[i + j*nx, 2] = i + (nx+1)*(j+1)
                    conn[i + j*nx, 3] = i+1 + (nx+1)*(j+1)

            if prob == 'MBB':
                for j in range(ny+1):
                    for i in range(nx+1):
                        X[i + j*(nx+1), 0] = lx*i/nx
                        X[i + j*(nx+1), 1] = ly*j/ny
                        if i > round(nx*(1-MBB_bc_portion)) and j == 0:
                            dof[i + j*(nx+1), 0] = ndof
                            ndof += 1
                        elif i == 0:
                            dof[i + j*(nx+1), 1] = ndof
                            ndof += 1
                        else:
                            dof[i + j*(nx+1), 0] = ndof
                            ndof += 1
                            dof[i + j*(nx+1), 1] = ndof
                            ndof += 1
            else:
                for j in range(ny+1):
                    for i in range(nx+1):
                        X[i + j*(nx+1), 0] = lx*i/nx
                        X[i + j*(nx+1), 1] = ly*j/ny
                        if i > 0:
                            dof[i + j*(nx+1), 0] = ndof
                            ndof += 1
                            dof[i + j*(nx+1), 1] = ndof
                            ndof += 1

        # force
        force = np.zeros(ndof)

        if prob == 'lbracket':
            if use_concentrated_force:
                force[dof[nx + ny1*(nx+1), 1]] = -force_magnitude
            else:
                nforce = round(ny1*forced_portion)
                for j in range(ny1+1-nforce, ny1+1):
                    force[dof[nx + j*(nx+1), 1]] = -force_magnitude
                force /= nforce

        elif prob == 'MBB':
            if use_concentrated_force:
                force[dof[ny*(nx+1), 1]] = -force_magnitude

            else:
                nforce = round(nx*forced_portion)
                for i in range(nforce):
                    force[dof[i + ny*(nx+1), 1]] = -force_magnitude
                force /= nforce
        elif prob == 'cantilever':

            if use_concentrated_force:
                force[dof[nx, 1]] = -force_magnitude

            else:
                nforce = round(ny*forced_portion)
                for j in range(nforce):
                    force[dof[nx + j*(nx+1), 1]] = -force_magnitude
                force /= nforce

        elif prob == 'michell':
            if use_concentrated_force:
                force[dof[nx+int(ny/2)*(nx+1), 1]] = -force_magnitude

            else:
                nforce = round(ny*forced_portion)
                if nforce % 2 == 0:
                    nforce += 1
                nhalf = int((nforce-1)/2)
                midj = int(ny/2)
                for j in range(midj-nhalf, midj+nhalf+1):
                    force[dof[nx + j*(nx+1), 1]] = -force_magnitude
                force /= nforce

    # Generate pickle file
    prob_pkl = dict()
    prob_pkl['prob_name'] = prob_name
    prob_pkl['nelems'] = nelems
    prob_pkl['nnodes'] = nnodes
    prob_pkl['ndof'] = ndof
    prob_pkl['C'] = C
    prob_pkl['conn'] = conn
    prob_pkl['X'] = X
    prob_pkl['dof'] = dof
    prob_pkl['force'] = force
    prob_pkl['r0'] = r0
    prob_pkl['density'] = density

    # Save pickle file
    outname = prob_pkl['prob_name']+'.pkl'
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
        outname = outdir + '/' + outname
    with open(outname, 'wb') as pklfile:
        pickle.dump(prob_pkl, pklfile)

    # Save mesh plot
    if plot_mesh:
        from plane_stress.utils import plot_mesh
        plot_mesh(prob_pkl, savefig=True, outdir=outdir)

    return

if __name__ == '__main__':

    # Set up parser
    p = argparse.ArgumentParser()
    p.add_argument('n', type=int,
        help='number of elements on short edge')
    p.add_argument('AR', type=float,
        help='domain aspect ratio, AR = length / height')
    p.add_argument('prob', type=str,
        choices=['cantilever', 'michell', 'MBB', 'lbracket'])
    p.add_argument('meshtype', type=str,
        choices=['structured', 'unstructured'])
    p.add_argument('--ratio1', type=float, default=0.4)
    p.add_argument('--ratio2', type=float, default=0.4)
    p.add_argument('--use_concentrated_force', action='store_true')
    p.add_argument('--use_hole', action='store_true')
    p.add_argument('--hole_radius', type=float, default=0.25,
        help='absolute radius = ly*hole_radius or ly1*hole_radius for lbracket')
    p.add_argument('--nr0', type=int, default=32,
        help='nr0 controls filter radius, r0 = height / nr0')
    p.add_argument('--outdir', type=str, default=None,
        help='directory for pkl output')
    p.add_argument('--no_plot_mesh', action='store_false')
    args = p.parse_args()


    # Set up constants
    force_magnitude = 25.0  # total force applied to structure
    forced_portion = 0.2  # portion of the edge of which load is applied to
    MBB_bc_portion = 0.1  # portion of the right-bottom corner where bc is applied to
    ly = 1.0  # for this project we always set height of domain to be 1.0
    density = 2700.0  # material density
    E = 70e3  # Young's modulus
    nu = 0.3  # Poisson's ratio

    # call preproc
    preproc(args.n, args.AR, args.prob, args.meshtype,
            args.ratio1, args.ratio2,args.use_concentrated_force,
            args.use_hole, args.hole_radius, args.nr0,
            args.outdir, args.no_plot_mesh, force_magnitude,
            forced_portion, MBB_bc_portion, ly, density, E, nu)