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

from compliance_minimization import ComplianceMinimization
from compliance_frequency import ComplianceFrequency

def plot_solution(X, u, conn):
    triangles = []
    for i in range(conn.shape[0]):
        triangles.append([conn[i, 0], conn[i, 1], conn[i, 2]])
        triangles.append([conn[i, 1], conn[i, 3], conn[i, 2]])

    # Create the triangles
    x = X[:,0]
    y = X[:,1]
    tri_obj = tri.Triangulation(x, y, triangles)

    # Plot the result as a figure
    fig, ax = plt.subplots()

    # Set the aspect ratio equal
    ax.set_aspect('equal')

    # Make sure that there are no ticks on either axis (these affect the bounding-box
    # and make extra white-space at the corners). Finally, turn off the axis so its
    # not visible.
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.axis('off')

    # Create the contour plot
    ax.tricontourf(tri_obj, u, cmap='coolwarm')

    plt.show()

    # # Save the figure. The bounding box is made tight to the figure, and the pading is
    # # determined via the pad_inches argument. Transparent sets the background transparent.
    # plt.savefig(outfile, dpi=500, transparent=True,
    #             bbox_inches='tight', pad_inches=0.01)

    # # Close the figure
    # plt.close()

nx = 64
ny = 64

r0 = 1.0/32.0

nelems = nx*ny
nnodes = (nx+1)*(ny+1)

conn = np.zeros((nelems, 4), dtype=np.intc)
vars = -np.ones((nnodes, 2), dtype=np.intc)
X = np.zeros((nnodes, 2))
rho = np.ones(nnodes)

C = np.zeros((3, 3))
qval = 5.0

E = 70e3
nu = 0.3
C[0, 0] = E/(1.0 - nu**2)
C[0, 1] = nu*E/(1.0 - nu**2)
C[1, 0] = C[0, 1]
C[1, 1] = C[0, 0]
C[2, 2] = 0.5*E/(1.0 + nu)

for j in range(ny):
    for i in range(nx):
        conn[i + j*nx, 0] = i + (nx+1)*j
        conn[i + j*nx, 1] = i+1 + (nx+1)*j
        conn[i + j*nx, 2] = i + (nx+1)*(j+1)
        conn[i + j*nx, 3] = i+1 + (nx+1)*(j+1)

nvars = 0
for j in range(ny+1):
    for i in range(nx+1):
        X[i + j*(nx+1), 0] = 1.0*i/nx
        X[i + j*(nx+1), 1] = 1.0*j/ny
        if i > 0:
            vars[i + j*(nx+1), 0] = nvars
            nvars += 1
            vars[i + j*(nx+1), 1] = nvars
            nvars += 1

# Specify the force and location
force = np.zeros(nvars)
i = nx
j = 0
force[vars[i + j*(nx+1), 1]] = -1.0e2

# problem = ComplianceMinimization(conn, vars, X, force, r0, qval, C)

density = 1.0
lambda0 = 1.0
ks_parameter = 100.0

problem = ComplianceFrequency(conn, vars, X, force, r0, qval, C,
    density, lambda0, ks_parameter)

problem.checkGradients()

options = {
    'algorithm': 'tr',
    'tr_init_size': 0.05,
    'tr_min_size': 1e-6,
    'tr_max_size': 10.0,
    'tr_eta': 0.25,
    'tr_infeas_tol': 1e-6,
    'tr_l1_tol': 1e-3,
    'tr_linfty_tol': 0.0,
    'tr_adaptive_gamma_update': True,
    'tr_max_iterations': 1000,
    'penalty_gamma': 10.0,
    'qn_subspace_size': 10,
    'qn_type': 'bfgs',
    'qn_diag_type': 'yts_over_sts',
    'abs_res_tol': 1e-8,
    'starting_point_strategy': 'affine_step',
    'barrier_strategy': 'mehrotra_predictor_corrector',
    'tr_steering_barrier_strategy':
        'mehrotra_predictor_corrector',
    'tr_steering_starting_point_strategy': 'affine_step',
    'use_line_search': False}

# Set up the optimizer
opt = ParOpt.Optimizer(problem, options)

# Set a new starting point
opt.optimize()
x, z, zw, zl, zu = opt.getOptimizedPoint()

rho = problem.F.dot(x[:])

plot_solution(problem.X, rho, problem.conn)
