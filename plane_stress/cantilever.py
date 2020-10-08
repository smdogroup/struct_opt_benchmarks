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


class ComplianceMinimization(ParOpt.Problem):
    def __init__(self, conn, vars, X, force, qval, C, r0):
        """
        The constructor for the topology optimization class.

        This function sets up the data that is requried to perform a
        plane stress analysis of a square, plane stress structure.
        This is probably only useful for topology optimization.
        """

        # Save the data
        self.conn = conn
        self.vars = vars
        self.X = X
        self.force = force
        self.qval = qval
        self.C = C

        # Set the number of variables and the number of nodes
        self.nvars = np.max(self.vars) + 1
        self.nnodes = np.max(self.conn) + 1

        super(ComplianceMinimization, self).__init__(MPI.COMM_SELF, self.nnodes, 1)

        # Compute the non-zero pattern for the sparse matrix
        rowp = np.zeros(self.nvars+1, dtype=np.intc)
        cols = np.zeros(1, dtype=np.intc)

        # Compute the dimension of the cols array required
        ncols = plane_stress.computenzpattern(conn.T, vars.T, rowp, cols)

        # Allocate the required dimension of the cols array
        cols_temp = np.zeros(ncols, dtype=np.intc)
        plane_stress.computenzpattern(self.conn.T, self.vars.T, rowp, cols_temp)

        # Truncate the cols array to only include
        self.cols = np.zeros(rowp[-1], dtype=np.intc)
        self.cols[:] = cols_temp[:rowp[-1]]
        self.rowp = rowp

        # Allocate space for the entries of the matrix
        self.Kvals = np.zeros(self.cols.shape)

        # Compute the mass (area) of the structure with a full density
        rho = np.ones(self.nnodes)
        self.total_mass = plane_stress.computemass(self.conn.T, self.X.T, rho)

        # Now, compute the filter weights and store them as a sparse
        # matrix
        F = sparse.lil_matrix((self.nnodes, self.nnodes))

        # Form a KDTree
        tree = KDTree(X)
        result = tree.query_ball_tree(tree, r0)

        for i, rlist in enumerate(result):
            w = []
            wvars = []
            for j in rlist:
                r = np.sqrt(np.dot(X[i,:] - X[j,:], X[i,:] - X[j,:]))
                if r < r0:
                    w.append((r0 - r)/r0)
                    wvars.append(j)

            # Normalize the weights
            w = np.array(w)
            w /= np.sum(w)

            # Set the weights into the filter matrix W
            F[i, wvars] = w

        # Covert the matrix to a CSR data format
        self.F = F.tocsr()

        return

    def mass(self, x):
        """
        Compute the mass of the structure
        """

        mass = plane_stress.computemass(self.conn.T, self.X.T, x)

        return mass

    def mass_grad(self, x):
        """
        Compute the derivative of the mass
        """
        dmdx = np.zeros(x.shape)
        plane_stress.computemassderiv(self.conn.T, self.X.T, dmdx)

        return dmdx

    def compliance(self, x):
        """
        Compute the structural compliance
        """

        # Compute the filtered compliance. Note that 'dot' is scipy
        # matrix-vector multiplicataion
        rho = self.F.dot(x)

        # Compute the stiffness matrix
        plane_stress.computekmat(self.conn.T, self.vars.T, self.X.T,
            self.qval, self.C.T, rho, self.rowp, self.cols, self.Kvals)

        # Form the matrix
        Kmat = sparse.csr_matrix((self.Kvals, self.cols, self.rowp), shape=(self.nvars, self.nvars))
        self.Kmat = Kmat.tocsc()
        self.LU = linalg.dsolve.factorized(self.Kmat)

        # Compute the solution to the linear system K*u = f
        self.u = self.LU(self.force)

        # Return the compliance
        return np.dot(self.force, self.u)

    def compliance_grad(self, x):
        """
        Compute the gradient of the compliance using the adjoint
        method.

        Since the governing equations are self-adjoint, and the
        function itself takes a special form:

        K*psi = 0.5*f => psi = 0.5*u

        So we can skip the adjoint computation itself since we have
        the displacement vector u from the solution.

        d(compliance)/dx = - 0.5*u^{T}*d(K*u - f)/dx = - 0.5*u^{T}*dK/dx*u
        """

        # Compute the filtered variables
        rho = self.F.dot(x)

        # First compute the derivative with respect to the filtered
        # variables
        dKdrho = np.zeros(x.shape)

        plane_stress.computekmatderiv(self.conn.T, self.vars.T, self.X.T,
            self.qval, self.C.T, rho, self.u, self.u, dKdrho)

        # Now evaluate the effect of the filter
        dcdx = -(self.F.transpose()).dot(dKdrho)

        return dcdx

    def getVarsAndBounds(self, x, lb, ub):
        """Get the variable values and bounds"""
        lb[:] = 1e-3
        ub[:] = 1.0
        x[:] = 0.95
        return

    def evalObjCon(self, x):
        """
        Return the objective, constraint and fail flag
        """

        fail = 0
        obj = self.compliance(x[:])
        con = np.array([0.4*self.total_mass - self.mass(x[:])])

        return fail, obj, con

    def evalObjConGradient(self, x, g, A):
        """
        Return the objective, constraint and fail flag
        """

        fail = 0
        g[:] = self.compliance_grad(x[:])
        A[0][:] = -self.mass_grad(x[:])

        # self.write_output(x[:])

        return fail

    # def write_output(self, x):
    #     """
    #     Write out something to the screen
    #     """

    #     if self.draw_figure:
    #         if not hasattr(self, 'fig'):
    #             plt.ion()
    #             self.fig, self.ax = plt.subplots()
    #             plt.draw()

    #         xfilter = self.F.dot(x)

    #         # Prepare a pixel visualization of the design vars
    #         image = np.zeros((self.nyelems, self.nxelems))
    #         for j in range(self.nyelems):
    #             for i in range(self.nxelems):
    #                 image[j, i] = xfilter[i + j*self.nxelems]

    #         x = np.linspace(0, self.Lx, self.nxelems)
    #         y = np.linspace(0, self.Ly, self.nyelems)

    #         self.ax.contourf(x, y, image)
    #         self.ax.set_aspect('equal', 'box')
    #         plt.draw()
    #         plt.pause(0.001)

    #     return

nx = 128
ny = 128

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

problem = ComplianceMinimization(conn, vars, X, force, qval, C, r0)

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


# rowp = np.zeros(nvars+1, dtype=np.intc)
# cols = np.zeros(1, dtype=np.intc)

# # Compute the dimension of the cols array required
# ncols = plane_stress.computenzpattern(conn.T, vars.T, rowp, cols)

# # Allocate the required dimension of the cols array
# cols_temp = np.zeros(ncols, dtype=np.intc)
# plane_stress.computenzpattern(conn.T, vars.T, rowp, cols_temp)

# # Truncate the cols array to only include
# cols = np.zeros(rowp[-1], dtype=np.intc)
# cols[:] = cols_temp[:rowp[-1]]

# # Allocate space for the entries of the matrix
# K = np.zeros(cols.shape)
# plane_stress.computekmat(conn.T, vars.T, X.T, qval, C.T, rho, rowp, cols, K)

# Kmat = sparse.csr_matrix((K, cols, rowp), shape=(nvars, nvars))
# Kmat = Kmat.tocsc()
# LU = linalg.dsolve.factorized(Kmat)

# f = np.ones(nvars)

# # Compute the solution to the linear system K*u = f
# u = LU(f)

# dfdx = np.zeros(nnodes)

# plane_stress.computekmatderiv(conn.T, vars.T, X.T, qval, C.T, rho, u, u, dfdx)



# uf = np.zeros(nnodes)
# for i in range(nnodes):
#     if vars[i,0] >= 0:
#         uf[i] = u[vars[i,0]]

# plot_solution(X, uf, conn)

# # plt.figure()
# for i in range(len(rowp)-1):
#     for jp in range(rowp[i], rowp[i+1]):
#         j = cols[jp]
#         plt.plot(i, j, 'o')
# plt.show()