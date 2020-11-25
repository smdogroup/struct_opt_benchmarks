import numpy as np
import openmdao.api as om
import matplotlib.pylab as plt
from plane_stress import solid_lib
from scipy import sparse
from scipy.sparse import linalg
from scipy.spatial import KDTree
import matplotlib.pylab as plt
import matplotlib.tri as tri
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pprint import pprint
from plane_stress.utils import plot_3dmesh
import pyamg
import timeit

class SolidAnalysis(om.ExplicitComponent):

    def __init__(self, conn, dof, X, force, r0, C, density,
                 qval, epsilon, ks_parameter,
                 design_stress=None, design_freq=None,
                 compute_comp=False, compute_mass=False,
                 compute_inertia=False, compute_freq=False,
                 compute_stress=False, use_pyamg=True):

        super().__init__()

        # Value check
        if compute_freq is True or design_freq is not None:
            print('\n[Warning] frequency computation is not supported for 3D problems yet!')

        if compute_inertia is True:
            print('\n[Warning] inertia computation is not supported for 3D problems yet!')

        # Save the data
        self.conn = conn
        self.dof = dof
        self.X = X
        self.force = force
        self.qval = qval
        self.C = C
        self.density = density
        self.qval = qval
        self.epsilon = epsilon
        self.ks_parameter = ks_parameter
        self.design_stress = design_stress
        self.compute_comp = compute_comp
        self.compute_mass = compute_mass
        self.compute_stress = compute_stress
        self.use_pyamg = use_pyamg

        # Set the number of variables and the number of nodes
        self.ndof = np.max(self.dof) + 1
        self.nnodes = np.max(self.conn) + 1
        self.nelems = conn.shape[0]

        # Compute the non-zero pattern for the sparse matrix
        rowp = np.zeros(self.ndof+1, dtype=np.intc)
        cols = np.zeros(1, dtype=np.intc)

        # Compute the dimension of the cols array required
        ncols = solid_lib.computenzpattern(conn.T, dof.T, rowp, cols)

        # Allocate the required dimension of the cols array
        cols_temp = np.zeros(ncols, dtype=np.intc)
        solid_lib.computenzpattern(self.conn.T, self.dof.T, rowp, cols_temp)

        # Truncate the cols array to only include
        self.cols = np.zeros(rowp[-1], dtype=np.intc)
        self.cols[:] = cols_temp[:rowp[-1]]
        self.rowp = rowp

        # Allocate space for the entries of the matrix
        self.Kvals = np.zeros(self.cols.shape)
        self.Mvals = np.zeros(self.cols.shape)

        # Compute the mass (area) of the structure with a full density
        rho = np.ones(self.nnodes)

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

        # Get adjacent element numbers and number of adjacent elements
        # of each node
        nsurelems = np.zeros(self.nnodes, dtype=np.intc)
        _surelems = [[] for i in range(self.nnodes)]
        for i in range(self.nelems):
            for j in range(8):
                _surelems[self.conn[i, j]].append(i)
                nsurelems[self.conn[i, j]] += 1

        # Convert surelems into 1D array
        surelems = []
        for i in _surelems:
            surelems.extend(i)
        surelems = np.array(surelems, dtype=np.intc)

        self.nsurelems = nsurelems
        self.surelems = surelems

    def setup(self):
        """
        The setup function required by openMDAO
        """

        # Inputs
        self.add_input('x', shape=(self.nnodes,), desc='topology design variable')

        # Outputs
        if self.compute_comp:
            self.add_output('c', shape=(1,), desc='compliance')

        if self.compute_mass:
            self.add_output('m', shape=(1,), desc='mass')

        if self.compute_stress:
            self.add_output('ks_nodal_stress', shape=(1,), desc='ks aggregation of nodal stress')

        # Declare partials
        if self.compute_comp:
            self.declare_partials(of='c', wrt='x')

        if self.compute_mass:
            self.declare_partials(of='m', wrt='x')

        if self.compute_stress:
            self.declare_partials(of='ks_nodal_stress', wrt='x')
        return

    def mass(self, x):
        """
        Compute the mass of the structure, since mass is linear,
        we always use density = 1.0 here, i.e. mass = volume
        """

        mass = solid_lib.computemass(self.conn.T, self.X.T, x)

        return mass

    def mass_grad(self, x):
        """
        Compute the derivative of the mass
        """
        dmdx = np.zeros(x.shape)
        solid_lib.computemassderiv(self.conn.T, self.X.T, dmdx)

        return dmdx

    def compliance(self, x):
        """
        Compute the structural compliance
        """

        # Compute the filtered compliance. Note that 'dot' is scipy
        # matrix-vector multiplicataion
        rho = self.F.dot(x)

        # Compute the stiffness matrix
        solid_lib.computekmat(self.conn.T, self.dof.T, self.X.T,
            self.qval, self.C.T, rho, self.rowp, self.cols, self.Kvals)

        # Form the matrix
        Kmat = sparse.csr_matrix((self.Kvals, self.cols, self.rowp),
                                 shape=(self.ndof, self.ndof))

        if self.use_pyamg:
            self.ml = pyamg.smoothed_aggregation_solver(Kmat, max_coarse=10)
            self.u = self.ml.solve(self.force, accel='cg', cycle='V', tol=1e-10)

        else:
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

        K*psi = f => psi = u

        So we can skip the adjoint computation itself since we have
        the displacement vector u from the solution.

        d(compliance)/dx = - u^{T}*d(K*u - f)/dx = - u^{T}*dK/dx*u
        """

        # Compute the filtered variables
        rho = self.F.dot(x)

        # First compute the derivative with respect to the filtered
        # variables
        dKdrho = np.zeros(x.shape)

        solid_lib.computekmatderiv(self.conn.T, self.dof.T, self.X.T,
            self.qval, self.C.T, rho, self.u, self.u, dKdrho)

        # Now evaluate the effect of the filter
        dcdx = -(self.F.transpose()).dot(dKdrho)

        return dcdx

    def nodal_stress(self, x):
        """
        Compute reconstructed nodal stress using either superconvergent
        patch recovery method
        Note that this nodal stress is dimensional
        """

        # Compute the filtered compliance. Note that 'dot' is scipy
        # matrix-vector multiplicataion
        rho = self.F.dot(x)

        # Compute the stiffness matrix
        solid_lib.computekmat(self.conn.T, self.dof.T, self.X.T,
            self.qval, self.C.T, rho, self.rowp, self.cols, self.Kvals)

        # Form the matrix
        Kmat = sparse.csr_matrix((self.Kvals, self.cols, self.rowp),
                                 shape=(self.ndof, self.ndof))
        self.Kmat = Kmat.tocsc()

        # Compute the solution to the linear system K*u = f
        if self.use_pyamg:
            self.ml = pyamg.smoothed_aggregation_solver(Kmat, max_coarse=10)
            self.u = self.ml.solve(self.force, accel='cg', cycle='V', tol=1e-10)
        else:
            self.LU = linalg.dsolve.factorized(self.Kmat)
            self.u = self.LU(self.force)

        # Compute nodal stress
        stress = np.zeros(self.nnodes)
        solid_lib.computenodalstress(self.surelems, self.nsurelems,
            self.conn.T, self.dof.T, self.X.T, self.epsilon,
            self.C.T, self.u, rho, stress)

        return stress

    def ks_nodal_stress(self, x):
        """
        Compute the KS approximation of the maximum nodal stress
        Note that this is the dimensionless stress
        """

        # design_stress must be given
        if self.design_stress is None:
            raise ValueError("\ndesign_stress must be specified to call this function!\n")

        # Compute nodal stress
        nodal_stress = self.nodal_stress(x)

        # Normalize by the material yield stress
        nodal_stress = nodal_stress / self.design_stress

        # Compute KS aggregation
        max_stress = np.max(nodal_stress)
        self.nodaletas = np.exp(self.ks_parameter*(nodal_stress - max_stress))
        etas_sum = np.sum(self.nodaletas)
        ks = max_stress + np.log(etas_sum) / self.ks_parameter
        self.nodaletas = self.nodaletas / etas_sum

        return ks

    def ks_nodal_stress_grad(self, x):
        """
        Compute the gradient of KS nodal stress
        Note that this is the gradient of dimensionless stress
        """

        # design_stress must be given
        if self.design_stress is None:
            raise ValueError("\ndesign_stress must be specified to call this function!\n")

        nelems = self.conn.shape[0]
        nnodes = self.nnodes
        # nsurelems = np.zeros(self.nnodes)

        # Compute the filtered variables
        rho = self.F.dot(x)

        # Compute derivative of ks w.r.t. quadrature stress
        dksds = np.zeros(8*nelems)
        dksdns = self.nodaletas
        solid_lib.computenodalstressderiv(self.surelems, self.nsurelems,
            self.conn.T, self.X.T, dksdns, dksds)

        # Compute the derivative of nodal stress w.r.t. quadrature stress
        dksds = dksds.reshape((-1, 8))

        # Compute dfdu
        dfdu = np.zeros(self.ndof)
        solid_lib.computestressstatederiv(self.conn.T, self.dof.T, self.X.T,
            self.epsilon, self.C.T, self.u, rho, dksds.T, dfdu)

        # Compute the adjoint variables
        if self.use_pyamg:
            psi = self.ml.solve(dfdu, accel='cg', cycle='V', tol=1e-10)
        else:
            psi = self.LU(dfdu)

        # Compute derivative of ks w.r.t. design variable
        dfdrho = np.zeros(nnodes)
        solid_lib.computestressderiv(self.conn.T, self.dof.T, self.X.T,
            self.epsilon, self.C.T, self.u, rho, dksds.T, dfdrho)

        # Compute the total derivative
        temp = np.zeros(self.nnodes)
        solid_lib.computekmatderiv(self.conn.T, self.dof.T, self.X.T,
            self.qval, self.C.T, rho, psi, self.u, temp)

        # Compute the remainder of the total derivative
        dfdrho -= temp

        dksdx = (self.F.transpose()).dot(dfdrho)

        dksdx /= self.design_stress

        return dksdx


    def compute(self, inputs, outputs):
        """
        The compute function required by openMDAO
        """

        # Inputs
        x = inputs['x']

        # Outputs
        if self.compute_comp:
            outputs['c'] = self.compliance(x[:])

        if self.compute_mass:
            outputs['m'] = self.mass(x[:])

        if self.compute_stress:
            outputs['ks_nodal_stress'] =self.ks_nodal_stress(x[:])

        return

    def compute_partials(self, inputs, partials):
        """
        The compute_partials function required by openMDAO
        """

        # Inputs
        x = inputs['x']

        # Derivatives
        if self.compute_comp:
            partials['c', 'x'] = self.compliance_grad(x[:])

        if self.compute_mass:
            partials['m', 'x'] = self.mass_grad(x[:])

        if self.compute_stress:
            partials['ks_nodal_stress', 'x'] = self.ks_nodal_stress_grad(x[:])

        return


# Test gradient implmentation using finite difference
# With a simple cantilever example
if __name__ == '__main__':

    print("\nRunning finite difference gradient checks...\n")

    nx = 2
    ny = 2
    nz = 4
    r0 = 1.0/6.0
    lx = 1.0
    ly = 1.0
    lz = 4.0
    force_magnitude = 25.0
    density = 2700.0

    nelems = nx*ny*nz
    nnodes = (nx+1)*(ny+1)*(nz+1)
    conn = np.zeros((nelems, 8), dtype=np.intc)
    dof = -np.ones((nnodes, 3), dtype=np.intc)
    X = np.zeros((nnodes, 3))
    rho = np.ones(nnodes)

    E = 70e3
    nu = 0.3
    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = C[2, 2] = 1 - nu
    C[0, 1] = C[0, 2] = C[1, 2] = nu
    C[1, 0] = C[2, 0] = C[2, 1] = nu
    C[3, 3] = C[4, 4] = C[5, 5] = 0.5 - nu
    C *= E/(1+nu)/(1-2*nu)

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                conn[i + j*nx + k*nx*ny, 0] = i + (nx+1)*j + (nx+1)*(ny+1)*k
                conn[i + j*nx + k*nx*ny, 1] = i+1 + (nx+1)*j + (nx+1)*(ny+1)*k
                conn[i + j*nx + k*nx*ny, 2] = i + (nx+1)*(j+1) + (nx+1)*(ny+1)*k
                conn[i + j*nx + k*nx*ny, 3] = i+1 + (nx+1)*(j+1) + (nx+1)*(ny+1)*k
                conn[i + j*nx + k*nx*ny, 4] = i + (nx+1)*j + (nx+1)*(ny+1)*(k+1)
                conn[i + j*nx + k*nx*ny, 5] = i+1 + (nx+1)*j + (nx+1)*(ny+1)*(k+1)
                conn[i + j*nx + k*nx*ny, 6] = i + (nx+1)*(j+1) + (nx+1)*(ny+1)*(k+1)
                conn[i + j*nx + k*nx*ny, 7] = i+1 + (nx+1)*(j+1) + (nx+1)*(ny+1)*(k+1)

    ndof = 0
    for k in range(nz+1):
        for j in range(ny+1):
            for i in range(nx+1):
                X[i + j*(nx+1)+k*(nx+1)*(ny+1), 0] = lx*i/nx
                X[i + j*(nx+1)+k*(nx+1)*(ny+1), 1] = ly*j/ny
                X[i + j*(nx+1)+k*(nx+1)*(ny+1), 2] = lz*k/nz
                if k > 0:
                    dof[i + j*(nx+1)+k*(nx+1)*(ny+1), 0] = ndof
                    ndof += 1
                    dof[i + j*(nx+1)+k*(nx+1)*(ny+1), 1] = ndof
                    ndof += 1
                    dof[i + j*(nx+1)+k*(nx+1)*(ny+1), 2] = ndof
                    ndof += 1

    force = np.zeros(ndof)
    k = nz
    nforce = 0
    for j in range(ny+1):
        for i in range(nx+1):
            force[dof[i + j*(nx+1)+k*(nx+1)*(ny+1), 1]] = -force_magnitude
            nforce += 1
    force /= nforce

    # We set random material for each element
    np.random.seed(0)
    x = np.random.rand(nnodes)

    # Plot mesh
    # plot_3dmesh(conn, X, dof, force)

    # Find out which one is better:
    N = 1

    analysis = SolidAnalysis(conn, dof, X, force, r0, C, density, qval=5.0,
                                epsilon=0.1, ks_parameter=100.0, design_stress=1.0,
                                use_pyamg=True, compute_comp=True, compute_mass=True, compute_stress=True)

    t1 = timeit.default_timer()
    for i in range(N):
        c2 = analysis.compliance(x)
    t2 = timeit.default_timer()
    print('compliance: {:e}'.format(c2))
    print('Pyamg runtime: {:e}'.format(t2-t1))

    # analysis = SolidAnalysis(conn, dof, X, force, r0, C, density, qval=5.0,
    #                             epsilon=0.1, ks_parameter=100.0, design_stress=1.0,
    #                             use_pyamg=False, compute_comp=True, compute_mass=True, compute_stress=True)


    # t1 = timeit.default_timer()
    # for i in range(N):
    #     c1 = analysis.compliance(x)
    # t2 = timeit.default_timer()
    # print('Scipy runtime: {:e}'.format(t2-t1))

    # print('relerr       : {:e}'.format((c2-c1)/c1))

    # Create openMDAO problem instance
    prob = om.Problem()
    indeps = prob.model.add_subsystem('indeps', om.IndepVarComp())
    indeps.add_output('x', x)
    prob.model.add_subsystem('topo', analysis)
    prob.model.connect('indeps.x', 'topo.x')

    # Execute
    prob.setup()
    prob.run_model()

    # Check partials using finite difference
    prob.check_partials()


