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

from tmr import TMR
import numpy as np
from mpi4py import MPI
import pickle
import argparse

# Set up parser
p = argparse.ArgumentParser()
p.add_argument('--n', type=int, default=64,
        help='meshsize = ly / n')
p.add_argument('--lx', type=float, default=1.0)
p.add_argument('--ly', type=float, default=1.0)
p.add_argument('--qval', type=float, default=5.0)
p.add_argument('--nr0', type=int, default=32,
        help='r0 = ly divided by nr0')
p.add_argument('--use_hole', action='store_true',
        help='put a circular hole in mesh')
p.add_argument('--type', type=str, default='cantilever',
        choices=['cantilever', 'michell', 'michelldistri', 'MBB'])
p.add_argument('--savevtk', action='store_true')
p.add_argument('--outdir', type=str, default='',
    help='directory for pkl output')
args = p.parse_args()

# prob_name
prob_type = 'unstruct-' + args.type
prob_name = '{:s}-n{:d}-lx{:.1f}-ly{:.1f}'.format(prob_type, args.n, args.lx, args.ly)

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

geo = create_cantilever(Lx=args.lx, Ly=args.ly, use_hole=args.use_hole)

# Create the mesh
comm = MPI.COMM_WORLD
mesh = TMR.Mesh(comm, geo)

# Mesh the part
opts = TMR.MeshOptions()
opts.frontal_quality_factor = 1.25
opts.num_smoothing_steps = 20
opts.write_mesh_quality_histogram = 1

# Mesh the geometry with the given target size
htarget = args.ly/args.n
mesh.mesh(htarget, opts=opts)
if args.savevtk:
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

# Set the nodal positions using only the x and y locations
X = np.array(X[:,:2])

# Get the nodes with the specified name
bcs = forest.getNodesWithName('6')

# Create the vars array with the number of nodes
nnodes = np.max(conn)+1
var = np.zeros((nnodes, 2), dtype=int)

# Set the vars to a negative index where we have a constraint
if args.type == 'MBB':
    # If it is MBB problem, we only fix x degree of freedom at left edge
    var[bcs, 0] = -1

    # We also want to fix y degree of freedom of lower-right corner node
    south_east_corner_node = -1
    xpos = X[0, 0]
    ypos = X[0, 0]
    for i in range(nnodes):
        if X[i, 0] >= xpos and X[i, 1] <= ypos:
            south_east_corner_node = i
            xpos = X[i, 0]
            ypos = X[i, 1]
    var[south_east_corner_node, 1] =-1

else:
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
forceval = 25.0
if args.type == 'cantilever':
    south_east_corner_node = -1
    xpos = X[0, 0]
    ypos = X[0, 0]
    for i in range(nnodes):
        if X[i, 0] >= xpos and X[i, 1] <= ypos:
            south_east_corner_node = i
            xpos = X[i, 0]
            ypos = X[i, 1]

    force[var[south_east_corner_node, 1]] = -forceval

elif args.type == 'michell-distri':
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

        force[v1] = -forceval
        force[v2] = -forceval

    force[:] /= np.count_nonzero(force)

elif args.type == 'MBB':
    north_west_corner_node = -1
    xpos = X[0, 0]
    ypos = X[0, 0]
    for i in range(nnodes):
        if X[i, 0] <= xpos and X[i, 1] >= ypos:
            north_west_corner_node = i
            xpos = X[i, 0]
            ypos = X[i, 1]

    force[var[north_west_corner_node, 1]] = -forceval

else:
    distance = args.lx**2 + args.ly**2
    xtarget = args.lx
    ytarget = args.ly / 2
    middle_node = -1
    for i in range(nnodes):
        xpos = X[i, 0]
        ypos = X[i, 1]
        if (xpos-xtarget)**2 + (ypos-ytarget)**2 <= distance:
            middle_node = i
            distance = (xpos-xtarget)**2 + (ypos-ytarget)**2

    force[var[middle_node, 1]] = -forceval

# r0
r0 = args.ly / args.nr0

# density
density = 2700.0

# qval
qval = args.qval

E = 70e3
nu = 0.3
C = np.zeros((3, 3))
C[0, 0] = E/(1.0 - nu**2)
C[0, 1] = nu*E/(1.0 - nu**2)
C[1, 0] = C[0, 1]
C[1, 1] = C[0, 0]
C[2, 2] = 0.5*E/(1.0 + nu)

prob_pkl = dict()
prob_pkl['prob_name'] = prob_name
prob_pkl['nelems'] = len(conn)
prob_pkl['nnodes'] = nnodes
prob_pkl['ndof'] = nvars
prob_pkl['C'] = C
prob_pkl['conn'] = conn
prob_pkl['X'] = X
prob_pkl['dof'] = var
prob_pkl['force'] = force
prob_pkl['r0'] = r0
prob_pkl['density'] = density
prob_pkl['qval'] = qval
prob_pkl['x'] = None
prob_pkl['opt_settings'] = None

outname = prob_pkl['prob_name'] + '.pkl'
with open(outname, 'wb') as pklfile:
    pickle.dump(prob_pkl, pklfile)


outname = prob_pkl['prob_name']+'.pkl'
if args.outdir != '':
    try:
        os.mkdir(args.outdir)
    except:
        pass
    outname = args.outdir + '/' + outname
with open(outname, 'wb') as pklfile:
    pickle.dump(prob_pkl, pklfile)