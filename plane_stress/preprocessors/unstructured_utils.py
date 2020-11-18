import numpy as np
from mpi4py import MPI
try:
    from tmr import TMR
except:
    print('\n[Warning] Cannot import tmr, unstructured mesh generator is unavailable!')

def create_geo(AR, prob, forced_portion, MBB_bc_portion,
    ratio1, ratio2, use_hole, hole_radius):

    # Create the surface in the x-y plane
    Ly = 1.0
    Lx = Ly*AR
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

    if prob == 'cantilever':

        # Create the vertices on the surface
        v1 = TMR.VertexFromFace(face, 0.0, 0.0)
        v2 = TMR.VertexFromFace(face, Lx, 0.0)
        v3 = TMR.VertexFromFace(face, Lx, forced_portion*Ly)
        v4 = TMR.VertexFromFace(face, Lx, Ly)
        v5 = TMR.VertexFromFace(face, 0.0, Ly)
        verts = [v1, v2, v3, v4, v5]

        # Set up the edges
        pcurve1 = TMR.BsplinePcurve(np.array([[0.0, 0.0], [Lx, 0.0]]))
        pcurve2 = TMR.BsplinePcurve(np.array([[Lx, 0.0], [Lx, forced_portion*Ly]]))
        pcurve3 = TMR.BsplinePcurve(np.array([[Lx, forced_portion*Ly], [Lx, Ly]]))
        pcurve4 = TMR.BsplinePcurve(np.array([[Lx, Ly], [0.0, Ly]]))
        pcurve5 = TMR.BsplinePcurve(np.array([[0.0, Ly], [0.0, 0.0]]))

        edge1 = TMR.EdgeFromFace(face, pcurve1)
        edge2 = TMR.EdgeFromFace(face, pcurve2)
        edge3 = TMR.EdgeFromFace(face, pcurve3)
        edge4 = TMR.EdgeFromFace(face, pcurve4)
        edge5 = TMR.EdgeFromFace(face, pcurve5)

        edge1.setVertices(v1, v2)
        edge2.setVertices(v2, v3)
        edge3.setVertices(v3, v4)
        edge4.setVertices(v4, v5)
        edge5.setVertices(v5, v1)

        edge1.setName('1')
        edge2.setName('2')
        edge3.setName('3')
        edge4.setName('4')
        edge5.setName('5')

        edges = [edge1, edge2, edge3, edge4, edge5]
        dirs = [1, 1, 1, 1, 1]
        loop = TMR.EdgeLoop([edge1, edge2, edge3, edge4, edge5], dirs)
        face.addEdgeLoop(1, loop)

    elif prob == 'michell':

        # Create the vertices on the surface
        v1 = TMR.VertexFromFace(face, 0.0, 0.0)
        v2 = TMR.VertexFromFace(face, Lx, 0.0)
        v3 = TMR.VertexFromFace(face, Lx, 0.5*(1-forced_portion)*Ly)
        v4 = TMR.VertexFromFace(face, Lx, 0.5*(1+forced_portion)*Ly)
        v5 = TMR.VertexFromFace(face, Lx, Ly)
        v6 = TMR.VertexFromFace(face, 0.0, Ly)
        verts = [v1, v2, v3, v4, v5, v6]

        # Set up the edges
        pcurve1 = TMR.BsplinePcurve(np.array([[0.0, 0.0], [Lx, 0.0]]))
        pcurve2 = TMR.BsplinePcurve(np.array([[Lx, 0.0], [Lx, 0.5*(1-forced_portion)*Ly]]))
        pcurve3 = TMR.BsplinePcurve(np.array([[Lx, 0.5*(1-forced_portion)*Ly], [Lx, 0.5*(1+forced_portion)*Ly]]))
        pcurve4 = TMR.BsplinePcurve(np.array([[Lx, 0.5*(1+forced_portion)*Ly], [Lx, Ly]]))
        pcurve5 = TMR.BsplinePcurve(np.array([[Lx, Ly], [0.0, Ly]]))
        pcurve6 = TMR.BsplinePcurve(np.array([[0.0, Ly], [0.0, 0.0]]))

        edge1 = TMR.EdgeFromFace(face, pcurve1)
        edge2 = TMR.EdgeFromFace(face, pcurve2)
        edge3 = TMR.EdgeFromFace(face, pcurve3)
        edge4 = TMR.EdgeFromFace(face, pcurve4)
        edge5 = TMR.EdgeFromFace(face, pcurve5)
        edge6 = TMR.EdgeFromFace(face, pcurve6)

        edge1.setVertices(v1, v2)
        edge2.setVertices(v2, v3)
        edge3.setVertices(v3, v4)
        edge4.setVertices(v4, v5)
        edge5.setVertices(v5, v6)
        edge6.setVertices(v6, v1)

        edge1.setName('1')
        edge2.setName('2')
        edge3.setName('3')
        edge4.setName('4')
        edge5.setName('5')
        edge6.setName('6')

        edges = [edge1, edge2, edge3, edge4, edge5, edge6]
        dirs = [1, 1, 1, 1, 1, 1]
        loop = TMR.EdgeLoop([edge1, edge2, edge3, edge4, edge5, edge6], dirs)
        face.addEdgeLoop(1, loop)

    elif prob == 'MBB':

        # Create the vertices on the surface
        v1 = TMR.VertexFromFace(face, 0.0, 0.0)
        v2 = TMR.VertexFromFace(face, Lx*(1-MBB_bc_portion), 0.0)
        v3 = TMR.VertexFromFace(face, Lx, 0.0)
        v4 = TMR.VertexFromFace(face, Lx, Ly)
        v5 = TMR.VertexFromFace(face, Lx*forced_portion, Ly)
        v6 = TMR.VertexFromFace(face, 0.0, Ly)
        verts = [v1, v2, v3, v4, v5, v6]

        # Set up the edges
        pcurve1 = TMR.BsplinePcurve(np.array([[0.0, 0.0], [Lx*(1-MBB_bc_portion), 0.0]]))
        pcurve2 = TMR.BsplinePcurve(np.array([[Lx*(1-MBB_bc_portion), 0.0], [Lx, 0.0]]))
        pcurve3 = TMR.BsplinePcurve(np.array([[Lx, 0.0], [Lx, Ly]]))
        pcurve4 = TMR.BsplinePcurve(np.array([[Lx, Ly], [Lx*forced_portion, Ly]]))
        pcurve5 = TMR.BsplinePcurve(np.array([[Lx*forced_portion, Ly], [0.0, Ly]]))
        pcurve6 = TMR.BsplinePcurve(np.array([[0.0, Ly], [0.0, 0.0]]))

        edge1 = TMR.EdgeFromFace(face, pcurve1)
        edge2 = TMR.EdgeFromFace(face, pcurve2)
        edge3 = TMR.EdgeFromFace(face, pcurve3)
        edge4 = TMR.EdgeFromFace(face, pcurve4)
        edge5 = TMR.EdgeFromFace(face, pcurve5)
        edge6 = TMR.EdgeFromFace(face, pcurve6)

        edge1.setVertices(v1, v2)
        edge2.setVertices(v2, v3)
        edge3.setVertices(v3, v4)
        edge4.setVertices(v4, v5)
        edge5.setVertices(v5, v6)
        edge6.setVertices(v6, v1)

        edge1.setName('1')
        edge2.setName('2')
        edge3.setName('3')
        edge4.setName('4')
        edge5.setName('5')
        edge6.setName('6')

        edges = [edge1, edge2, edge3, edge4, edge5, edge6]
        dirs = [1, 1, 1, 1, 1, 1]
        loop = TMR.EdgeLoop([edge1, edge2, edge3, edge4, edge5, edge6], dirs)
        face.addEdgeLoop(1, loop)

    elif prob == 'lbracket':

        # Create the vertices on the surface
        v1 = TMR.VertexFromFace(face, 0.0, 0.0)
        v2 = TMR.VertexFromFace(face, Lx, 0.0)
        v3 = TMR.VertexFromFace(face, Lx, Ly*ratio2*(1-forced_portion))
        v4 = TMR.VertexFromFace(face, Lx, Ly*ratio2)
        v5 = TMR.VertexFromFace(face, Lx*ratio1, Ly*ratio2)
        v6 = TMR.VertexFromFace(face, Lx*ratio1, Ly)
        v7 = TMR.VertexFromFace(face, 0.0, Ly)
        verts = [v1, v2, v3, v4, v5, v6, v7]

        # Set up the edges
        pcurve1 = TMR.BsplinePcurve(np.array([[0.0, 0.0], [Lx, 0.0]]))
        pcurve2 = TMR.BsplinePcurve(np.array([[Lx, 0.0], [Lx, Ly*ratio2*(1-forced_portion)]]))
        pcurve3 = TMR.BsplinePcurve(np.array([[Lx, Ly*ratio2*(1-forced_portion)], [Lx, Ly*ratio2]]))
        pcurve4 = TMR.BsplinePcurve(np.array([[Lx, Ly*ratio2], [Lx*ratio1, Ly*ratio2]]))
        pcurve5 = TMR.BsplinePcurve(np.array([[Lx*ratio1, Ly*ratio2], [Lx*ratio1, Ly]]))
        pcurve6 = TMR.BsplinePcurve(np.array([[Lx*ratio1, Ly], [0.0, Ly]]))
        pcurve7 = TMR.BsplinePcurve(np.array([[0.0, Ly], [0.0, 0.0]]))

        edge1 = TMR.EdgeFromFace(face, pcurve1)
        edge2 = TMR.EdgeFromFace(face, pcurve2)
        edge3 = TMR.EdgeFromFace(face, pcurve3)
        edge4 = TMR.EdgeFromFace(face, pcurve4)
        edge5 = TMR.EdgeFromFace(face, pcurve5)
        edge6 = TMR.EdgeFromFace(face, pcurve6)
        edge7 = TMR.EdgeFromFace(face, pcurve7)

        edge1.setVertices(v1, v2)
        edge2.setVertices(v2, v3)
        edge3.setVertices(v3, v4)
        edge4.setVertices(v4, v5)
        edge5.setVertices(v5, v6)
        edge6.setVertices(v6, v7)
        edge7.setVertices(v7, v1)

        edge1.setName('1')
        edge2.setName('2')
        edge3.setName('3')
        edge4.setName('4')
        edge5.setName('5')
        edge6.setName('6')
        edge7.setName('7')

        edges = [edge1, edge2, edge3, edge4, edge5, edge6, edge7]
        dirs = [1, 1, 1, 1, 1, 1, 1]
        loop = TMR.EdgeLoop([edge1, edge2, edge3, edge4, edge5, edge6, edge7], dirs)
        face.addEdgeLoop(1, loop)

    # Set up the hole
    if use_hole:

        if prob == 'lbracket':
            r = hole_radius*Ly*ratio2
            xc = 0.5*(1+ratio1)*Lx
            yc = 0.5*Ly*ratio2

        else:
            r = hole_radius *Ly
            xc = 0.5*Lx
            yc = 0.5*Ly

        vc = TMR.VertexFromFace(face, xc - r, yc)

        pts = [[-r, 0.0], [-r, r], [0.0, r], [r, r],
            [r, 0.0], [r, -r], [0.0, -r], [-r, -r], [-r, 0.0]]
        pts = np.array(pts)
        for i in range(pts.shape[0]):
            pts[i,0] += xc
            pts[i,1] += yc

        wts = [1.0, 1.0/np.sqrt(2), 1.0, 1.0/np.sqrt(2),
            1.0, 1.0/np.sqrt(2), 1.0, 1.0/np.sqrt(2), 1.0]
        Tu = [0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0]

        pcurvec = TMR.BsplinePcurve(np.array(pts),
            tu=np.array(Tu), wts=np.array(wts), k=3)

        edgec = TMR.EdgeFromFace(face, pcurvec)
        edgec.setVertices(vc, vc)

        loop = TMR.EdgeLoop([edgec], [1])
        face.addEdgeLoop(-1, loop)

        verts.append(vc)
        edges.append(edgec)


    # Create the TMRModel
    faces = [face]
    geo = TMR.Model(verts, edges, faces)

    return geo

def create_mesh(n, AR, prob, ratio1, ratio2, forced_portion, MBB_bc_portion,
    force_magnitude, use_concentrated_force, use_hole, hole_radius):

    Ly = 1.0
    Lx = Ly*AR

    # Create tmr geometry object
    geo = create_geo(AR, prob, forced_portion, MBB_bc_portion,
        ratio1, ratio2, use_hole, hole_radius)

    # Create the mesh
    comm = MPI.COMM_WORLD
    mesh = TMR.Mesh(comm, geo)

    # Mesh the part
    opts = TMR.MeshOptions()
    opts.frontal_quality_factor = 1.25
    opts.num_smoothing_steps = 20
    opts.write_mesh_quality_histogram = 1

    # Mesh the geometry with the given target size
    htarget = Ly / n
    mesh.mesh(htarget, opts=opts)

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
    if prob == 'cantilever':
        bcs = forest.getNodesWithName('5')
    elif prob == 'MBB':
        bc1 = forest.getNodesWithName('6')
        bc2 = forest.getNodesWithName('2')
    elif prob == 'michell':
        bcs = forest.getNodesWithName('6')
    elif prob == 'lbracket':
        bcs = forest.getNodesWithName('6')

    # Create the vars array with the number of nodes
    nnodes = np.max(conn)+1
    dof = np.zeros((nnodes, 2), dtype=int)

    # Set the vars to a negative index where we have a constraint
    if prob == 'MBB':

        # If it is MBB problem, we only fix x degree of freedom at left edge
        dof[bc1, 0] = -1
        dof[bc2, 1] = -1

    else:
        dof[bcs, :] = -1

    # Assign variables to the nodes
    ndof = 0
    for i in range(nnodes):
        if dof[i,0] >= 0:
            dof[i,0] = ndof
            ndof += 1
        if dof[i,1] >= 0:
            dof[i,1] = ndof
            ndof += 1

    # Set up the force vector
    force = np.zeros(ndof)

    # Assign the force vector
    if prob == 'cantilever':

        if use_concentrated_force:

            south_east_corner_node = -1
            xpos = X[0, 0]
            ypos = X[0, 0]
            for i in range(nnodes):
                if X[i, 0] >= xpos and X[i, 1] <= ypos:
                    south_east_corner_node = i
                    xpos = X[i, 0]
                    ypos = X[i, 1]

            force[dof[south_east_corner_node, 1]] = -force_magnitude

        else:

            quads_on_edge = forest.getQuadsWithName('2')

            for q in quads_on_edge:
                if q.info < 2:
                    n1 = q.info
                    n2 = q.info + 2
                else:
                    n1 = q.info//2
                    n2 = n1 + 1

                v1 = dof[conn[q.face, n1], 1]
                v2 = dof[conn[q.face, n2], 1]

                force[v1] = -force_magnitude
                force[v2] = -force_magnitude

            force[:] /= np.count_nonzero(force)

    elif prob == 'michell':

        if use_concentrated_force:

            distance = Lx**2 + Ly**2
            xtarget = Lx
            ytarget = Ly / 2
            middle_node = -1
            for i in range(nnodes):
                xpos = X[i, 0]
                ypos = X[i, 1]
                if (xpos-xtarget)**2 + (ypos-ytarget)**2 <= distance:
                    middle_node = i
                    distance = (xpos-xtarget)**2 + (ypos-ytarget)**2

            force[dof[middle_node, 1]] = -force_magnitude

        else:

            quads_on_edge = forest.getQuadsWithName('3')

            for q in quads_on_edge:
                if q.info < 2:
                    n1 = q.info
                    n2 = q.info + 2
                else:
                    n1 = q.info//2
                    n2 = n1 + 1

                v1 = dof[conn[q.face, n1], 1]
                v2 = dof[conn[q.face, n2], 1]

                force[v1] = -force_magnitude
                force[v2] = -force_magnitude

            force[:] /= np.count_nonzero(force)

    elif prob == 'MBB':

        if use_concentrated_force:

            distance = Lx**2 + Ly**2
            xtarget = 0
            ytarget = Ly
            north_west_node = -1
            for i in range(nnodes):
                xpos = X[i, 0]
                ypos = X[i, 1]
                if (xpos-xtarget)**2 + (ypos-ytarget)**2 <= distance:
                    north_west_node = i
                    distance = (xpos-xtarget)**2 + (ypos-ytarget)**2

            force[dof[north_west_node, 1]] = -force_magnitude

        else:

            quads_on_edge = forest.getQuadsWithName('5')

            for q in quads_on_edge:
                if q.info < 2:
                    n1 = q.info
                    n2 = q.info + 2
                else:
                    n1 = q.info//2
                    n2 = n1 + 1

                v1 = dof[conn[q.face, n1], 1]
                v2 = dof[conn[q.face, n2], 1]

                force[v1] = -force_magnitude
                force[v2] = -force_magnitude

            force[:] /= np.count_nonzero(force)


    elif prob == 'lbracket':

        if use_concentrated_force:

            distance = Lx**2 + Ly**2
            xtarget = Lx
            ytarget = Ly*ratio2
            middle_node = -1
            for i in range(nnodes):
                xpos = X[i, 0]
                ypos = X[i, 1]
                if (xpos-xtarget)**2 + (ypos-ytarget)**2 <= distance:
                    middle_node = i
                    distance = (xpos-xtarget)**2 + (ypos-ytarget)**2

            force[dof[middle_node, 1]] = -force_magnitude

        else:

            quads_on_edge = forest.getQuadsWithName('3')

            for q in quads_on_edge:
                if q.info < 2:
                    n1 = q.info
                    n2 = q.info + 2
                else:
                    n1 = q.info//2
                    n2 = n1 + 1

                v1 = dof[conn[q.face, n1], 1]
                v2 = dof[conn[q.face, n2], 1]

                force[v1] = -force_magnitude
                force[v2] = -force_magnitude

            force[:] /= np.count_nonzero(force)

    nelems = len(conn)

    return nelems, nnodes, ndof, conn, X, dof, force


if __name__ == '__main__':
    pass
