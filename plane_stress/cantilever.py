
import numpy as np
import matplotlib.pylab as plt
import plane_stress
from scipy import sparse
from scipy.sparse import linalg
import matplotlib.pylab as plt
import matplotlib.tri as tri

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


nx = 128
ny = 128

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

# Allocate space for the entries of the matrix
K = np.zeros(cols.shape)
plane_stress.computekmat(conn.T, vars.T, X.T, qval, C.T, rho, rowp, cols, K)

Kmat = sparse.csr_matrix((K, cols, rowp), shape=(nvars, nvars))
Kmat = Kmat.tocsc()
LU = linalg.dsolve.factorized(Kmat)

f = np.ones(nvars)

# Compute the solution to the linear system K*u = f
u = LU(f)

uf = np.zeros(nnodes)
for i in range(nnodes):
    if vars[i,0] >= 0:
        uf[i] = u[vars[i,0]]

plot_solution(X, uf, conn)

# plt.figure()
# for i in range(len(rowp)-1):
#     for jp in range(rowp[i], rowp[i+1]):
#         j = cols[jp]
#         plt.plot(i, j, 'o')
# plt.show()