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
"""

import numpy as np
import pickle
import argparse

# Set up parser
p = argparse.ArgumentParser()
p.add_argument('--nx', type=int, default=20)
p.add_argument('--ny', type=int, default=20)
p.add_argument('--lx', type=float, default=1.0)
p.add_argument('--ly', type=float, default=1.0)
p.add_argument('--outdir', type=str, default='',
    help='directory for pkl output')
args = p.parse_args()

# nelems and nnodes
nx = args.nx
ny = args.ny
lx = args.lx
ly = args.ly
nelems = nx*ny
nnodes = (nx+1)*(ny+1)

# prob_name
prob_name = 'cantilever-nx{:d}-ny{:d}-lx{:.1f}-ly{:.1f}'.format(nx, ny, lx, ly)

# r0
r0 = 1.0/40.0

# density
density = 2700.0

# qval
qval = 5.0

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
for j in range(ny):
    for i in range(nx):
        conn[i + j*nx, 0] = i + (nx+1)*j
        conn[i + j*nx, 1] = i+1 + (nx+1)*j
        conn[i + j*nx, 2] = i + (nx+1)*(j+1)
        conn[i + j*nx, 3] = i+1 + (nx+1)*(j+1)
ndof = 0
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
forceval = 25.0
force = np.zeros(ndof)
i = nx
j = 0
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
    outname = args.outdir + '/' + outname
with open(outname, 'wb') as pklfile:
    pickle.dump(prob_pkl, pklfile)
