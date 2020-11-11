"""
This script generates a 2D plane stress finite element analysis
problem, and stores the problem object to a python pickle file.

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
    |---qval:         SIMP penalty,            float
    |---x:            nodal design variable,   None or nnodes*1 ndarray
    |---opt_settings: optimization consts,     None or dict

domain:
L-bracket

    fixed boundary
    --------------
    |  nx1, lx1  |
    |            |ny2
    |            |
    |            |ly2           ^ force
    |            |              |
    |            |___nx2, lx2___|
    |                           |
    |                           |ny1
    |                           |
    |                           |ly1
    |___________________________|



"""

import numpy as np
import pickle
import argparse
import os



# Set up parser
p = argparse.ArgumentParser()
p.add_argument('--nx1', type=int, default=20)
p.add_argument('--nx2', type=int, default=30)
p.add_argument('--ny1', type=int, default=20)
p.add_argument('--ny2', type=int, default=30)
p.add_argument('--lx1', type=float, default=0.4)
p.add_argument('--lx2', type=float, default=0.6)
p.add_argument('--ly1', type=float, default=0.4)
p.add_argument('--ly2', type=float, default=0.6)
p.add_argument('--qval', type=float, default=5.0)
p.add_argument('--nr0', type=int, default=32,
        help='r0 = 1 divided by nr0')
p.add_argument('--outdir', type=str, default='',
    help='directory for pkl output')
args = p.parse_args()

# nelems and nnodes
nx1 = args.nx1
nx2 = args.nx2
ny1 = args.ny1
ny2 = args.ny2
lx1 = args.lx1
lx2 = args.lx2
ly1 = args.ly1
ly2 = args.ly2
nelems = nx1*ny1 + nx1*ny2 + nx2*ny1
nnodes = (nx1+1)*(ny1+1) + (nx1+1)*(ny2+1) + (nx2+1)*(ny1+1) - (nx1 + 1) - (ny1 + 1)

# prob_name
prob_type = 'lbracket'
prob_name = '{:s}-nx{:d}-{:d}-ny{:d}-{:d}-lx{:.1f}-{:.1f}-ly{:.1f}-{:.1f}'.format(
    prob_type, nx1, nx2, ny1, ny2, lx1, lx2, ly1, ly2)

# r0
r0 = 1.0 / args.nr0

# density
density = 2700.0

# qval
qval = args.qval

# C
C = np.zeros((3, 3))
E = 70e3
nu = 0.3
C[0, 0] = E/(1.0 - nu**2)
C[0, 1] = nu*E/(1.0 - nu**2)
C[1, 0] = C[0, 1]
C[1, 1] = C[0, 0]
C[2, 2] = 0.5*E/(1.0 + nu)

# ndof, conn, X, dof
conn = np.zeros((nelems, 4), dtype=np.intc)
dof = -np.ones((nnodes, 2), dtype=np.intc)
X = np.zeros((nnodes, 2))
ndof = 0
nx = nx1 + nx2
ny = ny1 + ny2
lx = lx1 + lx2
ly = ly1 + ly2

hx1 = lx1/nx1
hx2 = lx2/nx2
hy1 = ly1/ny1
hy2 = ly2/ny2
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

# force
forceval = 25.0
force = np.zeros(ndof)
i = nx
for j in range(int(0.75*(ny1+1)), ny1+1):
    force[dof[i + j*(nx+1), 1]] = -forceval

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
prob_pkl['qval'] = qval
prob_pkl['x'] = None
prob_pkl['opt_settings'] = None

outname = prob_pkl['prob_name']+'.pkl'
if args.outdir != '':
    try:
        os.mkdir(args.outdir)
    except:
        pass
    outname = args.outdir + '/' + outname
with open(outname, 'wb') as pklfile:
    pickle.dump(prob_pkl, pklfile)
