import numpy as np
import openmdao.api as om
import matplotlib.pylab as plt
import plane_stress
from scipy import sparse
from scipy.sparse import linalg
from scipy.spatial import KDTree
import matplotlib.pylab as plt
import matplotlib.tri as tri
from mpl_toolkits.axes_grid1 import make_axes_locatable

class PlaneStressAnalysis(om.ExplicitComponent):

    def __init__(self, conn, vars, X, force, r0, qval, C, density=1.0,
                 epsilon=1.0, ys=1.0, ks_parameter=100.0,
                 freqconstr=False, lambda0=0.0, num_eigs=8, eigshsigma=-100.0):
        super().__init__()

        # Save the data
        self.conn = conn
        self.vars = vars
        self.X = X
        self.force = force
        self.qval = qval
        self.C = C
        self.density = density
        self.epsilon = epsilon
        self.ys = ys
        self.freqconstr = freqconstr
        self.lambda0 = lambda0
        self.ks_parameter = ks_parameter
        self.num_eigs = num_eigs
        self.eigshsigma = eigshsigma

        # Set the number of variables and the number of nodes
        self.nvars = np.max(self.vars) + 1
        self.nnodes = np.max(self.conn) + 1
        self.nelems = conn.shape[0]

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
            for j in range(4):
                _surelems[self.conn[i, j]].append(i)
                nsurelems[self.conn[i, j]] += 1

        # Convert surelems into 1D array
        surelems = []
        for i in _surelems:
            surelems.extend(i)
        surelems = np.array(surelems, dtype=np.intc)

        self.nsurelems = nsurelems
        self.surelems = surelems

        return

    def setup(self):
        self.add_input('x', shape=(self.nnodes,),
                       desc='topology design variable')
        self.add_output('c', shape=(1,), desc='compliance')
        self.declare_partials(of='c', wrt='x')
        self.add_output('m', shape=(1,), desc='mass constraint')
        self.declare_partials(of='m', wrt='x')
        self.add_output('inertia', shape=(1,), desc='moment of inertia')
        self.declare_partials(of='inertia', wrt='x')
        # self.add_output('ks_stress', shape=(1,), desc='ks aggregation of stress')
        # self.declare_partials(of='ks_stress', wrt='x')
        self.add_output('ks_nodal_stress', shape=(1,),
                        desc='ks aggregation of nodal stress')
        self.declare_partials(of='ks_nodal_stress', wrt='x')
        if self.freqconstr is True:
            self.add_output('freq', shape=(1,),
                            desc='frequency constraint')
            self.declare_partials(of='freq', wrt='x')

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

    def inertia(self, x):
        """
        Compute the moment of inertia of the structure
        around center of mass
        """

        inertia = plane_stress.computemomofinertia(
            self.conn.T, self.X.T, x)

        return inertia

    def inertia_grad(self, x):
        """
        Compute the derivative of moment of inertia of
        the structure around center of mass
        """

        dinertiadx = np.zeros(x.shape)
        plane_stress.computemomofinertiaderiv(
            self.conn.T, self.X.T, x, dinertiadx)

        return dinertiadx

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
        Kmat = sparse.csr_matrix((self.Kvals, self.cols, self.rowp),
                                 shape=(self.nvars, self.nvars))
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

        plane_stress.computekmatderiv(self.conn.T, self.vars.T, self.X.T,
            self.qval, self.C.T, rho, self.u, self.u, dKdrho)

        # Now evaluate the effect of the filter
        dcdx = -(self.F.transpose()).dot(dKdrho)

        return dcdx

    def _ks_stress(self, x):
        """
        Compute the KS approximation of the maximum stress
        """

        # Compute the filtered compliance. Note that 'dot' is scipy
        # matrix-vector multiplicataion
        rho = self.F.dot(x)

        # Compute the stiffness matrix
        plane_stress.computekmat(self.conn.T, self.vars.T, self.X.T,
            self.qval, self.C.T, rho, self.rowp, self.cols, self.Kvals)

        # Form the matrix
        Kmat = sparse.csr_matrix((self.Kvals, self.cols, self.rowp),
                                 shape=(self.nvars, self.nvars))
        self.Kmat = Kmat.tocsc()
        self.LU = linalg.dsolve.factorized(self.Kmat)

        # Compute the solution to the linear system K*u = f
        self.u = self.LU(self.force)

        # Compute the stress values
        self.stress = np.zeros((self.conn.shape[0], 4))
        plane_stress.computestress(self.conn.T, self.vars.T, self.X.T,
            self.epsilon, self.C.T, self.u, rho, self.stress.T)

        # Normalize by the yield stress
        self.stress = self.stress/self.ys

        max_stress = np.max(self.stress)
        self.etas = np.exp(self.ks_parameter*(self.stress - max_stress))
        etas_sum = np.sum(self.etas)

        ks = max_stress + np.log(etas_sum)/self.ks_parameter
        self.etas = self.etas/etas_sum

        # Return the compliance
        return ks

    def _ks_stress_grad(self, x):
        """
        Compute the gradient of the approximate maximum stress
        """

        # Compute the filtered variables
        rho = self.F.dot(x)

        # Compute the derivative of the ks function with respect to rho
        dfdrho = np.zeros(self.nnodes)
        temp = np.zeros(self.nnodes)

        # Compute dfdu
        dfdu = np.zeros(self.nvars)
        plane_stress.computestressstatederiv(self.conn.T, self.vars.T, self.X.T,
            self.epsilon, self.C.T, self.u, rho, self.etas.T, dfdu)

        # Compute the adjoint variables
        psi = self.LU(dfdu)

        # Compute the derivative w.r.t.
        plane_stress.computestressderiv(self.conn.T, self.vars.T, self.X.T,
            self.epsilon, self.C.T, self.u, rho, self.etas.T, dfdrho)

        # Compute the total derivative
        plane_stress.computekmatderiv(self.conn.T, self.vars.T, self.X.T,
            self.qval, self.C.T, rho, psi, self.u, temp)

        # Compute the remainder of the total derivative
        dfdrho -= temp

        dksdx = (self.F.transpose()).dot(dfdrho)

        dksdx /= self.ys

        return dksdx

    def nodal_stress(self, x):
        """
        Compute reconstructed nodal stress using either superconvergent
        patch recovery method
        """

        # Compute the filtered compliance. Note that 'dot' is scipy
        # matrix-vector multiplicataion
        rho = self.F.dot(x)

        # Compute the stiffness matrix
        plane_stress.computekmat(self.conn.T, self.vars.T, self.X.T,
            self.qval, self.C.T, rho, self.rowp, self.cols, self.Kvals)

        # Form the matrix
        Kmat = sparse.csr_matrix((self.Kvals, self.cols, self.rowp),
                                 shape=(self.nvars, self.nvars))
        self.Kmat = Kmat.tocsc()
        self.LU = linalg.dsolve.factorized(self.Kmat)

        # Compute the solution to the linear system K*u = f
        self.u = self.LU(self.force)

        # Compute nodal stress
        stress = np.zeros(self.nnodes)
        plane_stress.computenodalstress(self.surelems, self.nsurelems,
            self.conn.T, self.vars.T, self.X.T, self.epsilon,
            self.C.T, self.u, rho, stress)

        return stress

    def ks_nodal_stress(self, x):
        """
        Compute the KS approximation of the maximum nodal stress
        """

        # Compute nodal stress
        nodal_stress = self.nodal_stress(x)

        # Normalize by the material yield stress
        nodal_stress = nodal_stress / self.ys

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
        """

        nelems = self.conn.shape[0]
        nnodes = self.nnodes
        # nsurelems = np.zeros(self.nnodes)

        # Compute the filtered variables
        rho = self.F.dot(x)

        # Compute derivative of ks w.r.t. quadrature stress
        dksds = np.zeros(4*nelems)
        dksdns = self.nodaletas
        plane_stress.computenodalstressderiv(self.surelems, self.nsurelems,
            self.conn.T, self.X.T, dksdns, dksds)

        # Compute the derivative of nodal stress w.r.t. quadrature stress
        dksds = dksds.reshape((-1, 4))

        # Compute dfdu
        dfdu = np.zeros(self.nvars)
        plane_stress.computestressstatederiv(self.conn.T, self.vars.T, self.X.T,
            self.epsilon, self.C.T, self.u, rho, dksds.T, dfdu)

        # Compute the adjoint variables
        psi = self.LU(dfdu)

        # Compute derivative of ks w.r.t. design variable
        dfdrho = np.zeros(nnodes)
        plane_stress.computestressderiv(self.conn.T, self.vars.T, self.X.T,
            self.epsilon, self.C.T, self.u, rho, dksds.T, dfdrho)

        # Compute the total derivative
        temp = np.zeros(self.nnodes)
        plane_stress.computekmatderiv(self.conn.T, self.vars.T, self.X.T,
            self.qval, self.C.T, rho, psi, self.u, temp)

        # Compute the remainder of the total derivative
        dfdrho -= temp

        dksdx = (self.F.transpose()).dot(dfdrho)

        dksdx /= self.ys

        return dksdx

    def frequency(self, x):

        # Apply the filter to obtain the filtered values
        rho = self.F.dot(x)

        # Compute the stiffness matrix
        plane_stress.computekmat(self.conn.T, self.vars.T, self.X.T,
            self.qval, self.C.T, rho, self.rowp, self.cols, self.Kvals)
        plane_stress.computemmat(self.conn.T, self.vars.T, self.X.T,
            self.density, rho, self.rowp, self.cols, self.Mvals)

        # Compute the A-matrix
        Avals = self.Kvals - self.lambda0*self.Mvals

        # Form the matrix
        Amat = sparse.csr_matrix((Avals, self.cols, self.rowp),
                                 shape=(self.nvars, self.nvars))

        # Find the smallest eigenvalues close to zero that have the smallest real part
        self.eigvals, self.eigvecs = linalg.eigsh(Amat, k=self.num_eigs, sigma=self.eigshsigma,
                                                  which='LM', tol=1e-8)

        # Compute the smallest eigenvalue
        etaf = np.exp(-self.ks_parameter*(self.eigvals - np.min(self.eigvals)))

        ksvalue = np.min(self.eigvals) - np.log(np.sum(etaf))/self.ks_parameter

        self.etaf = etaf/np.sum(etaf)

        # print("ks approx. min eigval: {:.2e}, real min eigval: {:.2e}".format(
        #     ksvalue, np.min(self.eigvals)))

        return ksvalue

    def frequency_grad(self, x):

        # Apply the filter to obtain the filtered values
        rho = self.F.dot(x)

        dfdx = np.zeros(x.shape)
        temp = np.zeros(x.shape)

        for i in range(self.num_eigs):
            # Compute the derivative of d(eigvec^{T}*K(x)*eigvec)
            psi = self.eigvecs[:,i]
            plane_stress.computekmatderiv(self.conn.T, self.vars.T, self.X.T,
                self.qval, self.C.T, rho, psi.T, psi.T, temp)
            dfdx += self.etaf[i]*temp

            # Compute the derivative of d(eigvec^{T}*M(x)*eigvec)
            plane_stress.computemmatderiv(self.conn.T, self.vars.T, self.X.T,
                self.density, psi.T, psi.T, temp)
            dfdx -= self.lambda0*self.etaf[i]*temp

        dfdx = (self.F.transpose()).dot(dfdx)

        return dfdx

    def base_frequencies(self, x):
        """
        Compute the first [self.num_eigs] frequencies of the design
        """
        # Apply the filter to obtain the filtered values
        rho = self.F.dot(x)

        # Compute the stiffness and mass matrix
        plane_stress.computekmat(self.conn.T, self.vars.T, self.X.T,
            self.qval, self.C.T, rho, self.rowp, self.cols, self.Kvals)
        plane_stress.computemmat(self.conn.T, self.vars.T, self.X.T,
            self.density, rho, self.rowp, self.cols, self.Mvals)
        Kmat = sparse.csr_matrix((self.Kvals, self.cols, self.rowp),
                                 shape=(self.nvars, self.nvars))
        Mmat = sparse.csr_matrix((self.Mvals, self.cols, self.rowp),
                                 shape=(self.nvars, self.nvars))

        # Find the smallest eigenvalues close to zero that have the smallest real part
        # Note that we set sigma = 0.0 here because for a physically-realistic system,
        # eigenvalues for generalized system (K, M) should be always positive
        self.eigvals2, self.eigvecs2 = linalg.eigsh(Kmat, k=self.num_eigs, M=Mmat,
                                                  sigma=0.0, which='LM', tol=1e-8)
        frequencies = np.sqrt(self.eigvals2)/(np.pi*2)
        # print("base frequency: {:.2e}".format(frequencies[0]))

        return frequencies

    def compute(self, inputs, outputs):
        x = inputs['x']
        outputs['c'] = self.compliance(x[:])
        outputs['m'] = self.mass(x[:])
        outputs['inertia'] = self.inertia(x[:])
        # outputs['ks_stress'] = self.ks_nodal_stress(x[:])
        outputs['ks_nodal_stress'] =self.ks_nodal_stress(x[:])
        if self.freqconstr is True:
            outputs['freq'] = self.frequency(x[:])

    def compute_partials(self, inputs, partials):
        x = inputs['x']
        partials['c', 'x'] = self.compliance_grad(x[:])
        partials['m', 'x'] = self.mass_grad(x[:])
        partials['inertia', 'x'] = self.inertia_grad(x[:])
        # partials['ks_stress', 'x'] = self.ks_nodal_stress_grad(x[:])
        partials['ks_nodal_stress', 'x'] = self.ks_nodal_stress_grad(x[:])
        if self.freqconstr is True:
            partials['freq', 'x'] = self.frequency_grad(x[:])

    def plot_solution(self, x, savefig=False, name='fig'):
        triangles = []
        for i in range(self.conn.shape[0]):
            triangles.append([self.conn[i, 0], self.conn[i, 1], self.conn[i, 2]])
            triangles.append([self.conn[i, 1], self.conn[i, 3], self.conn[i, 2]])

        # Create the triangles
        ptx = self.X[:,0]
        pty = self.X[:,1]
        tri_obj = tri.Triangulation(ptx, pty, triangles)

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
        rho = self.F.dot(x[:])
        cntr = ax.tricontourf(tri_obj, rho, cmap='coolwarm')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(cntr, cax=cax)
        fig.set_size_inches(6.4,5.0)
        fig.subplots_adjust(top=0.85)
        fig.subplots_adjust(bottom=0.05)

        # Compute compliance
        compliance = self.compliance(x)

        # Compute mass
        mass = self.mass(x)
        xfull = np.ones(self.nnodes)
        mass /= self.mass(xfull)

        # Compute base frequency
        base_freq = self.base_frequencies(x)[0]

        ax.set_title('compliance: {:.2e}\nnormalized mass: {:.2e}\nminimal frequency: {:.2e}'.
                      format(compliance, mass, base_freq))

        # Compute original eigenvalues
        eig = self.frequency(x)

        print("sigma: {:.2e}, ks_eigs(A): {:.2e}, minimal freq: {:.2e}".format(
            self.eigshsigma, eig, base_freq))

        # Plot the figure
        if savefig is False:
            plt.show()
        else:
            plt.savefig(name+'.png')
            plt.close()

    def plot_stress(self, x):
        """
        Generate stress contour plot
        """

        # Compute stress
        stress = self.nodal_stress(x)
        triangles = []
        for i in range(self.conn.shape[0]):
            triangles.append([self.conn[i, 0], self.conn[i, 1], self.conn[i, 2]])
            triangles.append([self.conn[i, 1], self.conn[i, 3], self.conn[i, 2]])

        # Create the triangles
        ptx = self.X[:,0]
        pty = self.X[:,1]
        tri_obj = tri.Triangulation(ptx, pty, triangles)

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
        cntr = ax.tricontourf(tri_obj, stress, cmap='coolwarm')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(cntr, cax=cax)
        fig.set_size_inches(6.4,5.0)
        fig.subplots_adjust(top=0.85)
        fig.subplots_adjust(bottom=0.05)

        # Compute mass
        mass = self.mass(x)
        xfull = np.ones(self.nnodes)
        mass /= self.mass(xfull)

        # Compute normalized maximum stress
        ks_stress = self.ks_nodal_stress(x)

        # Set title
        ax.set_title('mass: {:.2e}\nks stress: {:.2e}'.format(
            mass, ks_stress))


        plt.show()


    def rho(self, x):
        """
        Compute filtered variable rho out of x
        """

        # Apply the filter to obtain the filtered values
        rho = self.F.dot(x)

        return rho


# Test gradient implmentation using finite difference
# With a simple cantilever example
if __name__ == '__main__':

    print("\nRunning finite difference gradient checks...\n")

    # Define geometry, mesh and boundary conditions
    nx = 8
    ny = 8
    r0 = 1.0/32.0
    lx = 4.0
    ly = 1.0
    forceval = 25.0

    nelems = nx*ny
    nnodes = (nx+1)*(ny+1)

    conn = np.zeros((nelems, 4), dtype=np.intc)
    vars = -np.ones((nnodes, 2), dtype=np.intc)
    X = np.zeros((nnodes, 2))
    rho = np.ones(nnodes)

    C = np.zeros((3, 3))
    qval = 5.0

    density = 2700.0 # kg/m^3

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
            X[i + j*(nx+1), 0] = lx*i/nx
            X[i + j*(nx+1), 1] = ly*j/ny
            if i > 0:
                vars[i + j*(nx+1), 0] = nvars
                nvars += 1
                vars[i + j*(nx+1), 1] = nvars
                nvars += 1

    force = np.zeros(nvars)
    i = nx
    j = 0
    force[vars[i + j*(nx+1), 1]] = -forceval

    epsilon = 0.3
    ys = 200.0

    # Create analysis object instance
    analysis = PlaneStressAnalysis(conn, vars, X, force,
        r0, qval, C, density, epsilon, ys, freqconstr=True)

    # We set random material for each element
    np.random.seed(0)
    x = np.random.rand(nnodes)

    # compute mass and inertia
    mass = analysis.mass(x)
    inertia = analysis.inertia(x)

    # Create openMDAO problem instance
    prob = om.Problem()
    indeps = prob.model.add_subsystem('indeps', om.IndepVarComp())
    indeps.add_output('x', x)
    prob.model.add_subsystem('topo', analysis)
    prob.model.connect('indeps.x', 'topo.x')
    prob.model.add_design_var('indeps.x', lower=1e-3, upper=1.0)
    prob.model.add_constraint('topo.m', lower=0.0)
    prob.model.add_constraint('topo.c', lower=0.0)
    prob.model.add_constraint('topo.inertia', lower=0.0)
    prob.model.add_constraint('topo.freq', lower=0.0)

    # Execute
    prob.setup()
    prob.run_model()

    # Check partials using finite difference
    prob.check_partials()