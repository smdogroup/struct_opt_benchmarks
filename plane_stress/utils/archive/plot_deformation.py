from mpi4py import MPI
import numpy as np
import matplotlib.pylab as plt
import plane_stress
from scipy import sparse
from scipy.sparse import linalg
from scipy.spatial import KDTree
import matplotlib.pylab as plt
import matplotlib.tri as tri
from paropt import ParOpt
import timeit
import argparse

def plotShape(X, conn, u=None, exaggerate=1.0, title=None,
              savefig=False):
    for i in range(nelems):
        order = [0,1,3,2,0]
        for k in range(4):
            pt1 = X[conn[i, order[k  ]]]
            pt2 = X[conn[i, order[k+1]]]
            d1 = [0., 0.]
            d2 = [0., 0.]
            if u is not None:
                if (vars[conn[i, order[k]], 0] >= 0):
                    d1[0] = u[vars[conn[i, order[k]], 0]]*exaggerate
                if (vars[conn[i, order[k]], 1] >= 0):
                    d1[1] = u[vars[conn[i, order[k]], 1]]*exaggerate
                if (vars[conn[i, order[k+1]], 0] >= 0):
                    d2[0] = u[vars[conn[i, order[k+1]], 0]]*exaggerate
                if (vars[conn[i, order[k+1]], 1] >= 0):
                    d2[1] = u[vars[conn[i, order[k+1]], 1]]*exaggerate

            plt.plot([pt1[0]+d1[0], pt2[0]+d2[0]],
                     [pt1[1]+d1[1], pt2[1]+d2[1]],'k')

    plt.axis('equal')
    if title is not None:
        plt.title(title)
    if savefig is True:
        plt.savefig(title+'.png')
        plt.close()
    else:
        plt.show()


# Start main timer
tstart = timeit.default_timer()

# Take cmd line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--nx', type=int, default=64)
parser.add_argument('--ny', type=int, default=64)
parser.add_argument('--lx', type=float, default=1.0)
parser.add_argument('--ly', type=float, default=1.0)
args = parser.parse_args()

# Construct cantilever mesh
nxelems = args.nx
nyelems = args.ny
lx = args.lx
ly = args.ly
nxnodes = nxelems + 1
nynodes = nyelems + 1
nelems = nxelems * nyelems
nnodes = nxnodes * nynodes

# Construct connectivity
conn = np.zeros((nelems, 4), dtype=np.intc)
for j in range(nyelems):
    for i in range(nxelems):
        conn[i + j*nxelems, 0] = i + j*nxnodes
        conn[i + j*nxelems, 1] = i+1 + j*nxnodes
        conn[i + j*nxelems, 2] = i + (j+1)*nxnodes
        conn[i + j*nxelems, 3] = i+1 + (j+1)*nxnodes

# Construct nodal positions and vars list
X = np.zeros((nnodes, 2))
vars = -np.ones((nnodes, 2), dtype=np.intc)
nvars = 0
for j in range(nynodes):
    for i in range(nxnodes):
        X[i + j*nxnodes, 0] = lx * i/nxelems
        X[i + j*nxnodes, 1] = ly * j/nyelems
        if (i > 0):
            vars[i + j*nxnodes, 0] = nvars
            nvars += 1
            vars[i + j*nxnodes, 1] = nvars
            nvars += 1

# Form constitutive matrix
E = 70e3
nu = 0.3
density = 2700.0
C = np.zeros((3, 3))
C[0, 0] = E/(1.0 - nu**2)
C[0, 1] = nu*E/(1.0 - nu**2)
C[1, 0] = C[0, 1]
C[1, 1] = C[0, 0]
C[2, 2] = 0.5*E/(1.0 + nu)

# Specify the force and location
force = np.zeros(nvars)
j = 0
i = nxelems
force[vars[i + j*nxnodes, 1]] = -1e2

# Switch off filter
qval = 0.0
rho = np.ones(nnodes)

# Compute the non-zero pattern for the sparse matrix
rowp = np.zeros(nvars+1, dtype=np.intc)
cols = np.zeros(1, dtype=np.intc)

# Compute the dimension of the cols array required
ncols = plane_stress.computenzpattern(conn.T, vars.T, rowp, cols)

# Allocate the required dimension of the cols array
cols_temp = np.zeros(ncols, dtype=np.intc)
plane_stress.computenzpattern(conn.T, vars.T, rowp, cols_temp)

# Truncate the cols array to only include
cols = np.zeros(rowp[-1], dtype=np.intc)
cols[:] = cols_temp[:rowp[-1]]
rowp = rowp

# Construct global stiffness matrix K
Kvals = np.zeros(cols.shape)
plane_stress.computekmat(conn.T, vars.T, X.T, qval,
                         C.T, rho, rowp, cols, Kvals)
Kmat = sparse.csr_matrix((Kvals, cols, rowp), shape=(nvars, nvars))

# Solve for displacements
tstart_u = timeit.default_timer()
LU = linalg.dsolve.factorized(Kmat.tocsc())
u = LU(force)
tend_u = timeit.default_timer()

# Construct mass matrix M and matrix A
Mvals = np.zeros(cols.shape)
lambda0 = 20.0
Avals = Kvals - lambda0*Mvals
plane_stress.computemmat(conn.T, vars.T, X.T, density,
                         rho, rowp, cols, Mvals)
Mmat = sparse.csr_matrix((Mvals, cols, rowp), shape=(nvars, nvars))
Amat = sparse.csr_matrix((Avals, cols, rowp), shape=(nvars, nvars))
num_eigs = 8

# solve eigenvalue problem K*eigvec = eigval*M*eigvec
tstart_eig = timeit.default_timer()
eigvals, eigvecs = linalg.eigsh(Kmat, k=num_eigs, M=Mmat, sigma=0.0, which='LM')
tend_eig = timeit.default_timer()
freqs = np.sqrt(eigvals)

# solve eigenvalue problem A*eigvec = eigval*eigvec
tstart_eigA = timeit.default_timer()
eigvalsA, eigvecsA = linalg.eigsh(Amat, k=num_eigs, sigma=0.0, which='LA')
tend_eigA = timeit.default_timer()

# End main timer
tend = timeit.default_timer()

# Output results
print('degree of freedom:', nvars)
print('static solution wall time:         {:.3e} s'.format(tend_u-tstart_u))
print('eigenproblem solution wall time:   {:.3e} s'.format(tend_eig-tstart_eig))
print('eigenproblem A solution wall time: {:.3e} s'.format(tend_eigA-tstart_eigA))
print('total wall time:                   {:.3e} s'.format(tend-tstart))
print('eigenvalues :\n', eigvals)
print('A problem eigenvalues: \n', eigvalsA)
# for i in range(num_eigs):
#     plotShape(X, conn, eigvecs[:, i], exaggerate=100,
#               title='modal-shape-No-'+str(i+1))