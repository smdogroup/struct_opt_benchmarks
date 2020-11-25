#!/usr/bin/env python3

import numpy as np
import pickle
import argparse
import os

def preproc3d(n, AR, prob, meshtype, nr0, outdir, plot_mesh,
              force_magnitude, forced_portion, MBB_bc_portion,
              ly, density, E, nu):

    # prob_name
    prob_name = '3D-{:s}-{:s}-n{:d}-AR{:.1f}'.format(meshtype, prob, n, AR)

    # dimensions
    ny = n
    nx = round(ny*AR)
    nz = n
    lx = ly*AR
    lz = ly

    # nnodes, nelems and prob_name
    if prob == 'lbracket':

        ratio1 = ratio2 = 0.4
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

        if meshtype == 'structured':
            nelems_layer = nx1*ny1 + nx1*ny2 + nx2*ny1
            nnodes_layer = (nx1+1)*(ny1+1)+(nx1+1)*(ny2+1)+(nx2+1)*(ny1+1)-(nx1+1)-(ny1+1)
            nelems = nelems_layer * nz
            nnodes = nnodes_layer * (nz+1)

    else:

        if meshtype == 'structured':
            nelems_layer = nx*ny
            nnodes_layer = (nx+1)*(ny+1)
            nelems = nelems_layer * nz
            nnodes = nnodes_layer * (nz+1)

    # r0
    r0 = ly / nr0

    # C
    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = C[2, 2] = 1 - nu
    C[0, 1] = C[0, 2] = C[1, 2] = nu
    C[1, 0] = C[2, 0] = C[2, 1] = nu
    C[3, 3] = C[4, 4] = C[5, 5] = 0.5 - nu
    C *= E/(1+nu)/(1-2*nu)

    # We use tmr to generate unstructured mesh
    if meshtype == 'unstructured':
        pass

    else:
        # conn, dof, X
        conn = np.zeros((nelems, 8), dtype=np.intc)
        dof = -np.ones((nnodes, 3), dtype=np.intc)
        X = np.zeros((nnodes, 3))
        ndof = 0

        if prob == 'lbracket':
            nodei = 0
            for k in range(nz+1):
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

                            X[nodei, 2] = lz*k/nz

                            if j < ny:
                                dof[nodei, 0] = ndof
                                ndof += 1
                                dof[nodei, 1] = ndof
                                ndof += 1
                                dof[nodei, 2] = ndof
                                ndof += 1

                            nodei += 1

            elemi = 0
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        if i < nx1 or j < ny1:
                            if j <= ny1:
                                conn[elemi, 0] = i + (nx+1)*j + nnodes_layer*k
                                conn[elemi, 1] = i+1 + (nx+1)*j + nnodes_layer*k
                                conn[elemi, 2] = i + (nx+1)*(j+1) + nnodes_layer*k
                                conn[elemi, 3] = i+1 + (nx+1)*(j+1) + nnodes_layer*k
                                conn[elemi, 4] = i + (nx+1)*j + nnodes_layer*(k+1)
                                conn[elemi, 5] = i+1 + (nx+1)*j + nnodes_layer*(k+1)
                                conn[elemi, 6] = i + (nx+1)*(j+1) + nnodes_layer*(k+1)
                                conn[elemi, 7] = i+1 + (nx+1)*(j+1) + nnodes_layer*(k+1)

                            else:
                                conn[elemi, 0] = i + (nx1+1)*(j-ny1-1) + (nx+1)*(ny1+1) + nnodes_layer*k
                                conn[elemi, 1] = i+1 + (nx1+1)*(j-ny1-1) + (nx+1)*(ny1+1) + nnodes_layer*k
                                conn[elemi, 2] = i + (nx1+1)*(j-ny1) + (nx+1)*(ny1+1) + nnodes_layer*k
                                conn[elemi, 3] = i+1 + (nx1+1)*(j-ny1) + (nx+1)*(ny1+1) + nnodes_layer*k
                                conn[elemi, 4] = i + (nx1+1)*(j-ny1-1) + (nx+1)*(ny1+1) + nnodes_layer*(k+1)
                                conn[elemi, 5] = i+1 + (nx1+1)*(j-ny1-1) + (nx+1)*(ny1+1) + nnodes_layer*(k+1)
                                conn[elemi, 6] = i + (nx1+1)*(j-ny1) + (nx+1)*(ny1+1) + nnodes_layer*(k+1)
                                conn[elemi, 7] = i+1 + (nx1+1)*(j-ny1) + (nx+1)*(ny1+1) + nnodes_layer*(k+1)

                            elemi += 1

        else:
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        conn[i + j*nx + k*nelems_layer, 0] = i + (nx+1)*j + nnodes_layer*k
                        conn[i + j*nx + k*nelems_layer, 1] = i+1 + (nx+1)*j + nnodes_layer*k
                        conn[i + j*nx + k*nelems_layer, 2] = i + (nx+1)*(j+1) + nnodes_layer*k
                        conn[i + j*nx + k*nelems_layer, 3] = i+1 + (nx+1)*(j+1) + nnodes_layer*k
                        conn[i + j*nx + k*nelems_layer, 4] = i + (nx+1)*j + nnodes_layer*(k+1)
                        conn[i + j*nx + k*nelems_layer, 5] = i+1 + (nx+1)*j + nnodes_layer*(k+1)
                        conn[i + j*nx + k*nelems_layer, 6] = i + (nx+1)*(j+1) + nnodes_layer*(k+1)
                        conn[i + j*nx + k*nelems_layer, 7] = i+1 + (nx+1)*(j+1) + nnodes_layer*(k+1)

            if prob == 'MBB':
                for k in range(nz+1):
                    for j in range(ny+1):
                        for i in range(nx+1):
                            X[i + j*(nx+1) + nnodes_layer*k, 0] = lx*i/nx
                            X[i + j*(nx+1) + nnodes_layer*k, 1] = ly*j/ny
                            X[i + j*(nx+1) + nnodes_layer*k, 2] = lz*k/nz
                            nbc = round(nx*MBB_bc_portion)
                            if nbc < 2:
                                nbc = 2
                            if i > nx - nbc and j == 0:
                                dof[i + j*(nx+1) + nnodes_layer*k, 0] = ndof
                                ndof += 1
                            elif i == 0:
                                dof[i + j*(nx+1) + nnodes_layer*k, 1] = ndof
                                ndof += 1
                            else:
                                dof[i + j*(nx+1) + nnodes_layer*k, 0] = ndof
                                ndof += 1
                                dof[i + j*(nx+1) + nnodes_layer*k, 1] = ndof
                                ndof += 1
                                dof[i + j*(nx+1) + nnodes_layer*k, 2] = ndof
                                ndof += 1

            else:
                for k in range(nz+1):
                    for j in range(ny+1):
                        for i in range(nx+1):
                            X[i + j*(nx+1) + nnodes_layer*k, 0] = lx*i/nx
                            X[i + j*(nx+1) + nnodes_layer*k, 1] = ly*j/ny
                            X[i + j*(nx+1) + nnodes_layer*k, 2] = lz*k/nz
                            if i > 0:
                                dof[i + j*(nx+1) + nnodes_layer*k, 0] = ndof
                                ndof += 1
                                dof[i + j*(nx+1) + nnodes_layer*k, 1] = ndof
                                ndof += 1
                                dof[i + j*(nx+1) + nnodes_layer*k, 2] = ndof
                                ndof += 1


        # force
        force = np.zeros(ndof)
        if prob == 'lbracket':
            nforce = round(ny1*forced_portion)
            if nforce < 2:
                nforce = 2
            for k in range(nz+1):
                for j in range(ny1+1-nforce, ny1+1):
                    force[dof[nx + j*(nx+1) + nnodes_layer*k, 1]] = -force_magnitude
            force /= (nforce*(nz+1))

        elif prob == 'MBB':
            nforce = round(nx*forced_portion)
            if nforce < 2:
                nforce = 2
            for k in range(nz+1):
                for i in range(nforce):
                    force[dof[i + ny*(nx+1) + nnodes_layer*k, 1]] = -force_magnitude
            force /= (nforce*(nz+1))

        elif prob == 'cantilever':
            nforce = round(ny*forced_portion)
            if nforce < 2:
                nforce = 2
            for k in range(nz+1):
                for j in range(nforce):
                    force[dof[nx + j*(nx+1) + nnodes_layer*k, 1]] = -force_magnitude
            force /= (nforce*(nz+1))


        elif prob == 'michell':
            nforce = round(ny*forced_portion)
            midj = int(ny/2)
            if nforce <= 2:
                nforce = 2
                for k in range(nz+1):
                    for j in range(midj, midj+2):
                        force[dof[nx + j*(nx+1) + nnodes_layer*k, 1]] = -force_magnitude
            elif nforce > 2:
                nhalf = int((nforce-1)/2)
                for k in range(nz+1):
                    for j in range(midj-nhalf, midj+nhalf+1):
                        force[dof[nx + j*(nx+1) + nnodes_layer*k, 1]] = -force_magnitude
            force /= (nforce*(nz+1))


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
        from plane_stress.utils import plot_3dmesh
        plot_3dmesh(prob_pkl, savefig=True)

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
        choices=['structured'])
    p.add_argument('--nr0', type=int, default=6,
        help='nr0 controls filter radius, r0 = height / nr0')
    p.add_argument('--outdir', type=str, default=None,
        help='directory for pkl output')
    p.add_argument('--no_plot_mesh', action='store_false')
    args = p.parse_args()

    # Set up constants
    force_magnitude = 25.0  # total force applied to structure
    forced_portion = 0.3  # portion of the edge of which load is applied to
    MBB_bc_portion = 0.3  # portion of the right-bottom corner where bc is applied to
    ly = 1.0  # for this project we always set height of domain to be 1.0
    density = 2700.0  # material density
    E = 70e3  # Young's modulus
    nu = 0.3  # Poisson's ratio

    # call
    preproc3d(args.n, args.AR, args.prob, args.meshtype, args.nr0,
              args.outdir, args.no_plot_mesh, force_magnitude, forced_portion,
              MBB_bc_portion, ly, density, E, nu)