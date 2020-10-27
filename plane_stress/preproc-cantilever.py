"""
This script generates a 2D plane stress finite element analysis
problem, and stores the problem object to a python pickle file.

The structure of the problem object is as follows:

prob_pkl: dict
    |---prob_name: str
    |---nelems: int
    |---nnodes: int
    |---nnodes: int
    |---nvars: int
    |---C: 3*3 ndarray
    |---conn: nelems*4 ndarray
    |---X: nnodes*2 ndarray
    |---vars: nnodes*2 ndarray
    |---force: nvars*1 ndarray
    |---r0: float
    |---x: None or nnodes*1 ndarray
"""

import numpy as np
import pickle
import argparse

# Set up parser
p = argparse.ArgumentParser()
p.add_argument('--nx', type=int, default=64)
p.add_argument('--ny', type=int, default=64)
p.add_argument('--lx', type=float, default=1.0)
p.add_argument('--ly', type=float, default=1.0)
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
r0 = 2*np.max([lx/nx, ly/ny])

# C
C = np.zeros((3, 3))
E = 70e3
nu = 0.3
C[0, 0] = E/(1.0 - nu**2)
C[0, 1] = nu*E/(1.0 - nu**2)
C[1, 0] = C[0, 1]
C[1, 1] = C[0, 0]
C[2, 2] = 0.5*E/(1.0 + nu)

# nvars, conn, X, vars
conn = np.zeros((nelems, 4), dtype=np.intc)
vars = -np.ones((nnodes, 2), dtype=np.intc)
X = np.zeros((nnodes, 2))
for j in range(ny):
    for i in range(nx):
        conn[i + j*nx, 0] = i + (nx+1)*j
        conn[i + j*nx, 1] = i+1 + (nx+1)*j
        conn[i + j*nx, 2] = i + (nx+1)*(j+1)
        conn[i + j*nx, 3] = i+1 + (nx+1)*(j+1)
nvars = 0
for j in range(ny+1):
    for i in range(nx+1):
        X[i + j*(nx+1), 0] = lx*i/nx
        X[i + j*(nx+1), 1] = ly*j/ny
        if i > 0:
            vars[i + j*(nx+1), 0] = nvars
            nvars += 1
            vars[i + j*(nx+1), 1] = nvars
            nvars += 1

# force
forceval = 25.0
force = np.zeros(nvars)
i = nx
j = 0
force[vars[i + j*(nx+1), 1]] = -forceval

# Generate pickle file
prob_pkl = dict()
prob_pkl['prob_name'] = prob_name
prob_pkl['nelems'] = nelems
prob_pkl['nnodes'] = nnodes
prob_pkl['nvars'] = nvars
prob_pkl['C'] = C
prob_pkl['conn'] = conn
prob_pkl['X'] = X
prob_pkl['vars'] = vars
prob_pkl['force'] = force
prob_pkl['r0'] = r0
prob_pkl['x'] = None

outname = prob_pkl['prob_name']+'.pkl'
with open(outname, 'wb') as pklfile:
    pickle.dump(prob_pkl, pklfile)
