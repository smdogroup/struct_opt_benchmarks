import numpy as np
import openmdao.api as om
import matplotlib.pylab as plt
from plane_stress import plane_lib
from scipy import sparse
from scipy.sparse import linalg
from scipy.spatial import KDTree
import matplotlib.pylab as plt
import matplotlib.tri as tri
from mpl_toolkits.axes_grid1 import make_axes_locatable

class PlaneStressAnalysis(om.ExplicitComponent):

    def __init__(self, conn, dof, X, force, r0, C, density,
                 qval, epsilon, ks_parameter,
                 design_stress=None, design_freq=None,
                 compute_comp=False, compute_mass=False,
                 compute_inertia=False, compute_freq=False,
                 compute_stress=False, num_eigs=8,eigsh_sigma=-100.0):
        """
        Initialize the analysis class object
        """

        super().__init__()

        # Save the data
        self.conn = conn
        self.dof = dof
        self.X = X
        self.force = force
        self.qval = qval
        self.C = C
        self.density = density
        self.epsilon = epsilon
        self.design_stress = design_stress
        self.design_freq = design_freq
        self.ks_parameter = ks_parameter
        self.num_eigs = num_eigs
        self.eigsh_sigma = eigsh_sigma
        self.compute_comp = compute_comp
        self.compute_mass = compute_mass
        self.compute_inertia = compute_inertia
        self.compute_freq = compute_freq
        self.compute_stress = compute_stress

        # Set the number of variables and the number of nodes
        self.ndof = np.max(self.dof) + 1
        self.nnodes = np.max(self.conn) + 1
        self.nelems = conn.shape[0]

        # Compute the non-zero pattern for the sparse matrix
        rowp = np.zeros(self.ndof+1, dtype=np.intc)
        cols = np.zeros(1, dtype=np.intc)

        # Compute the dimension of the cols array required
        ncols = plane_lib.computenzpattern(conn.T, dof.T, rowp, cols)

        # Allocate the required dimension of the cols array
        cols_temp = np.zeros(ncols, dtype=np.intc)
        plane_lib.computenzpattern(self.conn.T, self.dof.T, rowp, cols_temp)

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

        if self.compute_inertia:
            self.add_output('inertia', shape=(1,), desc='moment of inertia')

        if self.compute_stress:
            self.add_output('ks_nodal_stress', shape=(1,), desc='ks aggregation of nodal stress')

        if self.compute_freq:
            self.add_output('freqc', shape=(1,), desc='frequency constraint')

        # Declare partials
        if self.compute_comp:
            self.declare_partials(of='c', wrt='x')

        if self.compute_mass:
            self.declare_partials(of='m', wrt='x')

        if self.compute_inertia:
            self.declare_partials(of='inertia', wrt='x')

        if self.compute_stress:
            self.declare_partials(of='ks_nodal_stress', wrt='x')

        if self.compute_freq:
            self.declare_partials(of='freqc', wrt='x')

    def mass(self, x):
        """
        Compute the mass of the structure, since mass is linear,
        we always use density = 1.0 here, i.e. mass = area
        """

        mass = plane_lib.computemass(self.conn.T, self.X.T, x)

        return mass

    def mass_grad(self, x):
        """
        Compute the derivative of the mass
        """
        dmdx = np.zeros(x.shape)
        plane_lib.computemassderiv(self.conn.T, self.X.T, dmdx)

        return dmdx

    def inertia(self, x):
        """
        Compute the moment of inertia of the structure around
        center of mass, since inertia is linear, we always use
        density = 1.0 here
        """

        inertia = plane_lib.computemomofinertia(
            self.conn.T, self.X.T, x)

        return inertia

    def inertia_grad(self, x):
        """
        Compute the derivative of moment of inertia of
        the structure around center of mass
        """

        dinertiadx = np.zeros(x.shape)
        plane_lib.computemomofinertiaderiv(
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
        plane_lib.computekmat(self.conn.T, self.dof.T, self.X.T,
            self.qval, self.C.T, rho, self.rowp, self.cols, self.Kvals)

        # Form the matrix
        Kmat = sparse.csr_matrix((self.Kvals, self.cols, self.rowp),
                                 shape=(self.ndof, self.ndof))
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

        plane_lib.computekmatderiv(self.conn.T, self.dof.T, self.X.T,
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
        plane_lib.computekmat(self.conn.T, self.dof.T, self.X.T,
            self.qval, self.C.T, rho, self.rowp, self.cols, self.Kvals)

        # Form the matrix
        Kmat = sparse.csr_matrix((self.Kvals, self.cols, self.rowp),
                                 shape=(self.ndof, self.ndof))
        self.Kmat = Kmat.tocsc()
        self.LU = linalg.dsolve.factorized(self.Kmat)

        # Compute the solution to the linear system K*u = f
        self.u = self.LU(self.force)

        # Compute nodal stress
        stress = np.zeros(self.nnodes)
        plane_lib.computenodalstress(self.surelems, self.nsurelems,
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
        dksds = np.zeros(4*nelems)
        dksdns = self.nodaletas
        plane_lib.computenodalstressderiv(self.surelems, self.nsurelems,
            self.conn.T, self.X.T, dksdns, dksds)

        # Compute the derivative of nodal stress w.r.t. quadrature stress
        dksds = dksds.reshape((-1, 4))

        # Compute dfdu
        dfdu = np.zeros(self.ndof)
        plane_lib.computestressstatederiv(self.conn.T, self.dof.T, self.X.T,
            self.epsilon, self.C.T, self.u, rho, dksds.T, dfdu)

        # Compute the adjoint variables
        psi = self.LU(dfdu)

        # Compute derivative of ks w.r.t. design variable
        dfdrho = np.zeros(nnodes)
        plane_lib.computestressderiv(self.conn.T, self.dof.T, self.X.T,
            self.epsilon, self.C.T, self.u, rho, dksds.T, dfdrho)

        # Compute the total derivative
        temp = np.zeros(self.nnodes)
        plane_lib.computekmatderiv(self.conn.T, self.dof.T, self.X.T,
            self.qval, self.C.T, rho, psi, self.u, temp)

        # Compute the remainder of the total derivative
        dfdrho -= temp

        dksdx = (self.F.transpose()).dot(dfdrho)

        dksdx /= self.design_stress

        return dksdx

    def freq_constr(self, x):
        """
        This computes the frequency constraint value freqc,
        note that this is not the actual base frequency of structure,
        instead, freqc >= 0.0 equals base frequency >= design frequency
        """

        # design_freq must be given
        if self.design_freq is None:
            raise ValueError("\design_freq must be specified to call this function!\n")

        # Apply the filter to obtain the filtered values
        rho = self.F.dot(x)

        # Compute the stiffness matrix
        plane_lib.computekmat(self.conn.T, self.dof.T, self.X.T,
            self.qval, self.C.T, rho, self.rowp, self.cols, self.Kvals)
        plane_lib.computemmat(self.conn.T, self.dof.T, self.X.T,
            self.density, rho, self.rowp, self.cols, self.Mvals)

        # Compute the A-matrix
        lambda0 = (2.0*np.pi*self.design_freq)**2
        Avals = self.Kvals - lambda0*self.Mvals

        # Form the matrix
        Amat = sparse.csr_matrix((Avals, self.cols, self.rowp),
                                 shape=(self.ndof, self.ndof))

        # Find the smallest eigenvalues close to zero that have the smallest real part
        self.eigvals, self.eigvecs = linalg.eigsh(Amat, k=self.num_eigs, sigma=self.eigsh_sigma,
                                                  which='LM', tol=1e-8)

        # Compute the smallest eigenvalue
        etaf = np.exp(-self.ks_parameter*(self.eigvals - np.min(self.eigvals)))

        freqc = np.min(self.eigvals) - np.log(np.sum(etaf))/self.ks_parameter

        self.etaf = etaf/np.sum(etaf)

        return freqc

    def freq_constr_grad(self, x):
        """
        Computes the gradient of frequency constraint value freqc
        """

        # design_freq must be given
        if self.design_freq is None:
            raise ValueError("\design_freq must be specified to call this function!\n")

        # Apply the filter to obtain the filtered values
        rho = self.F.dot(x)

        dfreqcdx = np.zeros(x.shape)
        temp = np.zeros(x.shape)

        for i in range(self.num_eigs):
            # Compute the derivative of d(eigvec^{T}*K(x)*eigvec)
            psi = self.eigvecs[:,i]
            plane_lib.computekmatderiv(self.conn.T, self.dof.T, self.X.T,
                self.qval, self.C.T, rho, psi.T, psi.T, temp)
            dfreqcdx += self.etaf[i]*temp

            # Compute the derivative of d(eigvec^{T}*M(x)*eigvec)
            plane_lib.computemmatderiv(self.conn.T, self.dof.T, self.X.T,
                self.density, psi.T, psi.T, temp)
            lambda0 = (2.0*np.pi*self.design_freq)**2
            dfreqcdx -= lambda0*self.etaf[i]*temp

        dfreqcdx = (self.F.transpose()).dot(dfreqcdx)

        return dfreqcdx

    def base_frequency(self, x):
        """
        Compute the base frequency of the design
        """
        # Apply the filter to obtain the filtered values
        rho = self.F.dot(x)

        # Compute the stiffness and mass matrix
        plane_lib.computekmat(self.conn.T, self.dof.T, self.X.T,
            self.qval, self.C.T, rho, self.rowp, self.cols, self.Kvals)
        plane_lib.computemmat(self.conn.T, self.dof.T, self.X.T,
            self.density, rho, self.rowp, self.cols, self.Mvals)
        Kmat = sparse.csr_matrix((self.Kvals, self.cols, self.rowp),
                                 shape=(self.ndof, self.ndof))
        Mmat = sparse.csr_matrix((self.Mvals, self.cols, self.rowp),
                                 shape=(self.ndof, self.ndof))

        # Find the smallest eigenvalues close to zero that have the smallest real part
        # Note that we set sigma = 0.0 here because for a physically-realistic system,
        # eigenvalues for generalized system (K, M) should be always positive
        self.eigvals2, self.eigvecs2 = linalg.eigsh(Kmat, k=self.num_eigs, M=Mmat,
                                                  sigma=0.0, which='LM', tol=1e-8)

        frequencies = np.sqrt(self.eigvals2)/(np.pi*2)
        base_freq = frequencies[0]

        return base_freq

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

        if self.compute_inertia:
            outputs['inertia'] = self.inertia(x[:])

        if self.compute_stress:
            outputs['ks_nodal_stress'] =self.ks_nodal_stress(x[:])

        if self.compute_freq:
            outputs['freqc'] = self.freq_constr(x[:])

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

        if self.compute_inertia:
            partials['inertia', 'x'] = self.inertia_grad(x[:])

        if self.compute_stress:
            partials['ks_nodal_stress', 'x'] = self.ks_nodal_stress_grad(x[:])

        if self.compute_freq:
            partials['freqc', 'x'] = self.freq_constr_grad(x[:])

        return

    def plot_topology(self, x, savefig=False, filename='topology.png', paperstyle=False):
        """
        Plot the topology of design x
        """

        # Create the triangles
        triangles = []
        for i in range(self.conn.shape[0]):
            triangles.append([self.conn[i, 0], self.conn[i, 1], self.conn[i, 2]])
            triangles.append([self.conn[i, 1], self.conn[i, 3], self.conn[i, 2]])

        ptx = self.X[:,0]
        pty = self.X[:,1]
        tri_obj = tri.Triangulation(ptx, pty, triangles)

        length = max(ptx) - min(ptx)
        height = max(pty) - min(pty)

        fig_height = 4.8
        fig_length = 4.8 * length / height

        # Plot the result as a figure
        if paperstyle:
            fig, ax = plt.subplots(figsize=(fig_length, fig_height), constrained_layout=True)
        else:
            fig, ax = plt.subplots(figsize=(6.4, 5.0))

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

        if not paperstyle:
            fig.subplots_adjust(top=0.80)
            fig.subplots_adjust(bottom=0.02)

        # Set colorbar
        if not paperstyle:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(cntr, cax=cax)

        # Compute compliance
        compliance = self.compliance(x)

        # Compute mass
        xfull = np.ones(self.nnodes)
        mass = self.mass(x)
        normalized_mass = mass / self.mass(xfull)

        # Compute base frequency
        base_freq = self.base_frequency(x)

        # Compute maximum stress
        stress = self.nodal_stress(x)

        # Compute problem specific values
        if self.design_stress is not None:
            ks_stress = '{:.4e}'.format(self.ks_nodal_stress(x))
        else:
            ks_stress = '-'

        if self.design_freq is not None:
            freqc = '{:.4e}'.format(self.freq_constr(x))
        else:
            freqc = '-'

        # Set title
        if not paperstyle:
            title = '{:<15s}{:<12.6f}{:<15s}{:<12.6f}\n' \
                    '{:<15s}{:<12.6f}{:<15s}{:<12.6f}\n' \
                    '{:<15s}{:<12s}{:<15s}{:<12s}\n' \
                    '{:<10s}{:<8.2f}{:<10s}{:<8.2f}{:<10s}{:<8.2f}'.format(
                    'compliance:', compliance, 'norm-ed mass:', normalized_mass,
                    'base freq:', base_freq, 'max stress:', np.max(stress),
                    'freq constr:', freqc, 'ks stress:', ks_stress,
                    'SIMP qval:', self.qval, 'epsilon:', self.epsilon, 'ks param:', self.ks_parameter)

            ax.set_title(title, loc='left')

        # Plot or save the figure
        if savefig is False:
            plt.show()

        else:
            plt.savefig(filename)
            plt.close()

        return

    def plot_stress(self, x, savefig=False, filename='stress.png'):
        """
        Plot stress contour plot
        """

        # Compute stress
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
        stress = self.nodal_stress(x)
        cntr = ax.tricontourf(tri_obj, stress, cmap='coolwarm')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(cntr, cax=cax)
        fig.set_size_inches(6.4,5.0)
        fig.subplots_adjust(top=0.80)
        fig.subplots_adjust(bottom=0.02)

        # Compute compliance
        compliance = self.compliance(x)

        # Compute mass
        xfull = np.ones(self.nnodes)
        mass = self.mass(x)
        normalized_mass = mass / self.mass(xfull)

        # Compute base frequency
        base_freq = self.base_frequency(x)

        # Compute maximum stress
        stress = self.nodal_stress(x)

        # Compute problem specific values
        if self.design_stress is not None:
            ks_stress = '{:.4e}'.format(self.ks_nodal_stress(x))
        else:
            ks_stress = '-'

        if self.design_freq is not None:
            freqc = '{:.4e}'.format(self.freq_constr(x))
        else:
            freqc = '-'

        # Set title
        title = '{:<15s}{:<12.6f}{:<15s}{:<12.6f}\n' \
                '{:<15s}{:<12.6f}{:<15s}{:<12.6f}\n' \
                '{:<15s}{:<12s}{:<15s}{:<12s}\n' \
                '{:<10s}{:<8.2f}{:<10s}{:<8.2f}{:<10s}{:<8.2f}'.format(
                'compliance:', compliance, 'norm-ed mass:', normalized_mass,
                'base freq:', base_freq, 'max stress:', np.max(stress),
                'freq constr:', freqc, 'ks stress:', ks_stress,
                'SIMP qval:', self.qval, 'epsilon:', self.epsilon, 'ks param:', self.ks_parameter)

        ax.set_title(title, loc='left')

        # Plot or save the figure
        if savefig is False:
            plt.show()

        else:
            plt.savefig(filename)
            plt.close()

        return


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
    dof = -np.ones((nnodes, 2), dtype=np.intc)
    X = np.zeros((nnodes, 2))
    rho = np.ones(nnodes)

    C = np.zeros((3, 3))

    density = 2700.0

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

    force = np.zeros(ndof)
    i = nx
    j = 0
    force[dof[i + j*(nx+1), 1]] = -forceval


    # Create analysis object instance
    analysis = PlaneStressAnalysis(conn, dof, X, force,
        r0, C, density, qval=3.0, epsilon=0.1, ks_parameter=50.0,
        design_stress=200.0, design_freq=0.0, compute_comp=True,
        compute_mass=True, compute_inertia=True, compute_freq=True,
        compute_stress=True)

    # We set random material for each element
    np.random.seed(0)
    x = np.random.rand(nnodes)

    analysis.mass(x)

    # # Create openMDAO problem instance
    # prob = om.Problem()
    # indeps = prob.model.add_subsystem('indeps', om.IndepVarComp())
    # indeps.add_output('x', x)
    # prob.model.add_subsystem('topo', analysis)
    # prob.model.connect('indeps.x', 'topo.x')

    # # Execute
    # prob.setup()
    # prob.run_model()

    # # Check partials using finite difference
    # prob.check_partials()
