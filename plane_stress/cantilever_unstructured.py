from mpi4py import MPI
from tmr import TMR
import numpy as np
import matplotlib.pylab as plt
import plane_stress
import matplotlib.tri as tri
from paropt import ParOpt

from compliance_minimization import ComplianceMinimization
from compliance_frequency import ComplianceFrequency
from stress_minimization import StressMinimization

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

def create_cantilever(Lx=2.0, Ly=1.0, use_hole=False):

    # Create the surface in the x-y plane
    nu = 2
    nv = 2
    x = np.linspace(0.0, Lx, nu)
    y = np.linspace(0.0, Ly, nv)
    pts = np.zeros((nu, nv, 3))
    for j in range(nv):
        for i in range(nu):
            pts[i,j,0] = x[i]
            pts[i,j,1] = y[j]

    tu = np.array([0.0, 0.0, Lx, Lx])
    tv = np.array([0.0, 0.0, Ly, Ly])

    # Create the b-spline surface
    surf = TMR.BsplineSurface(pts, tu=tu, tv=tv)
    face = TMR.FaceFromSurface(surf)

    # Create the vertices on the surface
    v1 = TMR.VertexFromFace(face, 0.0, 0.0)
    v2 = TMR.VertexFromFace(face, Lx, 0.0)
    v3 = TMR.VertexFromFace(face, Lx, (0.5 - 0.125)*Ly)
    v4 = TMR.VertexFromFace(face, Lx, (0.5 + 0.125)*Ly)
    v5 = TMR.VertexFromFace(face, Lx, Ly)
    v6 = TMR.VertexFromFace(face, 0.0, Ly)
    verts = [v1, v2, v3, v4, v5, v6]

    # Set up the edges
    pcurve1 = TMR.BsplinePcurve(np.array([[0.0, 0.0], [Lx, 0.0]]))
    edge1 = TMR.EdgeFromFace(face, pcurve1)
    edge1.setVertices(v1, v2)
    edge1.setName('1')

    pcurve2 = TMR.BsplinePcurve(np.array([[Lx, 0.0], [Lx, (0.5 - 0.125)*Ly]]))
    edge2 = TMR.EdgeFromFace(face, pcurve2)
    edge2.setVertices(v2, v3)
    edge2.setName('2')

    pcurve3 = TMR.BsplinePcurve(np.array([[Lx, (0.5 - 0.125)*Ly], [Lx, (0.5 + 0.125)*Ly]]))
    edge3 = TMR.EdgeFromFace(face, pcurve3)
    edge3.setVertices(v3, v4)
    edge3.setName('3')

    pcurve4 = TMR.BsplinePcurve(np.array([[Lx, (0.5 + 0.125)*Ly], [Lx, Ly]]))
    edge4 = TMR.EdgeFromFace(face, pcurve4)
    edge4.setVertices(v4, v5)
    edge4.setName('4')

    pcurve5 = TMR.BsplinePcurve(np.array([[Lx, Ly], [0.0, Ly]]))
    edge5 = TMR.EdgeFromFace(face, pcurve5)
    edge5.setVertices(v5, v6)
    edge5.setName('5')

    pcurve6 = TMR.BsplinePcurve(np.array([[0.0, Ly], [0.0, 0.0]]))
    edge6 = TMR.EdgeFromFace(face, pcurve6)
    edge6.setVertices(v6, v1)
    edge6.setName('6')
    edges = [edge1, edge2, edge3, edge4, edge5, edge6]

    dirs = [1, 1, 1, 1, 1, 1]
    loop = TMR.EdgeLoop([edge1, edge2, edge3, edge4, edge5, edge6], dirs)
    face.addEdgeLoop(1, loop)

    # Set the points for the hole location
    if use_hole:
        r = 0.25*Ly
        v7 = TMR.VertexFromFace(face, 0.5*Lx - r, 0.5*Ly)

        pts = [[-r, 0.0], [-r, r], [0.0, r], [r, r],
            [r, 0.0], [r, -r], [0.0, -r], [-r, -r], [-r, 0.0]]
        pts = np.array(pts)
        for i in range(pts.shape[0]):
            pts[i,0] += 0.5*Lx
            pts[i,1] += 0.5*Ly

        wts = [1.0, 1.0/np.sqrt(2), 1.0, 1.0/np.sqrt(2),
            1.0, 1.0/np.sqrt(2), 1.0, 1.0/np.sqrt(2), 1.0]
        Tu = [0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0]
        pcurve7 = TMR.BsplinePcurve(np.array(pts),
                                    tu=np.array(Tu), wts=np.array(wts), k=3)
        edge7 = TMR.EdgeFromFace(face, pcurve7)
        edge7.setVertices(v7, v7)

        loop = TMR.EdgeLoop([edge7], [1])
        face.addEdgeLoop(-1, loop)

        verts.append(v7)
        edges.append(edge7)

    faces = [face]

    # Create the TMRModel
    geo = TMR.Model(verts, edges, faces)
    return geo

geo = create_cantilever(use_hole=True)

# Create the mesh
comm = MPI.COMM_WORLD
mesh = TMR.Mesh(comm, geo)

# Mesh the part
opts = TMR.MeshOptions()
opts.frontal_quality_factor = 1.25
opts.num_smoothing_steps = 20
opts.write_mesh_quality_histogram = 1

# Mesh the geometry with the given target size
n = 48
htarget = 1.0/n
mesh.mesh(htarget, opts=opts)
mesh.writeToVTK('surface-mesh.vtk')

# Create a model from the mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
forest = TMR.QuadForest(comm)
forest.setTopology(topo)

# Create random trees and balance the mesh. Print the output file
forest.createTrees()
forest.balance(1)

# Create the nodes
forest.createNodes()

# Get the mesh connectivity
conn = forest.getMeshConn()

# Get the node locations
X = forest.getPoints()

# Get the nodes with the specified name
bcs = forest.getNodesWithName('6')

# Create the vars array with the number of nodes
nnodes = np.max(conn)+1
var = np.zeros((nnodes, 2), dtype=int)

# Set the vars to a negative index where we have a constraint
var[bcs, :] = -1

# Assign variables to the nodes
nvars = 0
for i in range(nnodes):
    if var[i,0] >= 0:
        var[i,0] = nvars
        nvars += 1
    if var[i,1] >= 0:
        var[i,1] = nvars
        nvars += 1

# Set up the force vector
force = np.zeros(nvars)

# Assign the force vector
quads_on_edge = forest.getQuadsWithName('3')

for q in quads_on_edge:
    if q.info < 2:
        n1 = q.info
        n2 = q.info + 2
    else:
        n1 = q.info//2
        n2 = n1 + 1

    v1 = var[conn[q.face, n1], 1]
    v2 = var[conn[q.face, n2], 1]

    force[v1] += 0.5
    force[v2] += 0.5


r0 = 1.0/32.0
qval = 5.0

E = 70e3
nu = 0.3
C = np.zeros((3, 3))
C[0, 0] = E/(1.0 - nu**2)
C[0, 1] = nu*E/(1.0 - nu**2)
C[1, 0] = C[0, 1]
C[1, 1] = C[0, 0]
C[2, 2] = 0.5*E/(1.0 + nu)

# Set the nodal positions using only the x and y locations
X = np.array(X[:,:2])

problem = ComplianceMinimization(conn, var, X, force, r0, qval, C)

# problem = ComplianceFrequency(conn, var, X, force, r0, qval, C,
#     density, lambda0, ks_parameter)

problem.checkGradients()

options = {
    'algorithm': 'tr',
    'output_level': 2,
    'norm_type': 'l1',
    'tr_init_size': 0.05,
    'tr_min_size': 0.01,
    'tr_max_size': 10.0,
    'tr_eta': 0.25,
    'tr_infeas_tol': 1e-6,
    'tr_l1_tol': 0.0, # 1e-3,
    'tr_linfty_tol': 0.0,
    'tr_adaptive_gamma_update': False,
    'tr_max_iterations': 500,
    'penalty_gamma': 50.0,
    'qn_subspace_size': 2,
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

problem.checkGradients(x=x)

rho = problem.F.dot(x[:])

plot_solution(problem.X, rho, problem.conn)