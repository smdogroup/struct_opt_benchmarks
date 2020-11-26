import numpy as np
from mpi4py import MPI

try:
    from tmr import TMR
except:
    print('\n[Warning] Cannot import tmr, 3d unstructured mesh generator is unavailable!')

try:
    from egads4py import egads
except:
    print('\n[Warning] Cannot import egads4py, 3d unstructured mesh generator is unavailable!')

def create_3d_unstruct_mesh(n, AR, prob, ratio1, ratio2, hole_r, forced_portion,
                            force_magnitude, MBB_bc_portion, loaded_thickness, ly):

    comm = MPI.COMM_WORLD

    # Create the egads context
    ctx = egads.context()

    # dimensions
    Lx = ly*AR
    Ly = ly
    h  = ly

    parts = []

    if prob == 'lbracket':
        # Create box 1
        x0 = [0, 0, 0]
        x1 = [Lx*ratio1, Ly, h]
        B1 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])
        parts.append(ctx.makeTopology(egads.MODEL, children=[B1]))

        # Create box 2
        x0 = [Lx*ratio1, 0, 0]
        x1 = [Lx*(1-ratio1), Ly*ratio2, h]
        B2 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])

        # Create the cylinder cutout
        xc = Lx*(ratio1+1)*0.5
        yc = Ly*ratio2*0.5
        x0 = [xc, yc, 0]
        x1 = [xc, yc, h]
        r = hole_r*Ly*ratio2
        C1 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r])
        parts.append(B2.solidBoolean(C1, egads.SUBTRACTION))

    else:
        # Create the box
        x0 = [0, 0, 0]
        x1 = [Lx, Ly, h]
        B1 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])

        # Create the cylinder cutout
        xc = Lx/2
        yc = Ly/2
        x0 = [xc, yc, 0]
        x1 = [xc, yc, h]
        r = hole_r*Ly
        C1 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r])
        # parts.append(ctx.makeTopology(egads.MODEL, children=[B1]))
        parts.append(B1.solidBoolean(C1, egads.SUBTRACTION))

    # Create all of the models
    geos = []
    for p in parts:
        geos.append(TMR.ConvertEGADSModel(p))

    # Create the full list of vertices, edges, faces and volumes
    verts = []
    edges = []
    faces = []
    vols = []
    for geo in geos:
        verts.extend(geo.getVertices())
        edges.extend(geo.getEdges())
        faces.extend(geo.getFaces())
        vols.extend(geo.getVolumes())

    # Set all of the matching faces
    TMR.setMatchingFaces(geos)

    # Create the geometry
    geo = TMR.Model(verts, edges, faces, vols)

    # Create the new mesh
    mesh = TMR.Mesh(comm, geo)

    # Set the meshing options
    opts = TMR.MeshOptions()
    # opts.mesh_type_default = TMR.UNSTRUCTURED
    opts.write_mesh_quality_histogram = 1
    opts.triangularize_print_iter = 50000

    # Create the surface mesh
    htarget = ly / n
    mesh.mesh(htarget, opts)

    # Write the surface mesh to a file
    # mesh.writeToVTK('block.vtk', 'hex')

    # Create the model from the unstructured volume mesh
    model = mesh.createModelFromMesh()

    # Create the corresponding mesh topology from the mesh-model
    topo = TMR.Topology(comm, model)

    # Create the quad forest and set the topology of the forest
    forest = TMR.OctForest(comm)
    forest.setTopology(topo)
    forest.createTrees()
    forest.balance(1)

    # Create the nodes
    forest.createNodes()

    # Get the mesh connectivity
    conn = forest.getMeshConn()
    nnodes = np.max(conn) + 1
    nelems = len(conn)

    # Get the node locations
    X = forest.getPoints()

    # Set boundary conditions
    dof = - np.ones((nnodes, 3), dtype=int)
    geo_bc_tol = 1e-6

    if prob == 'cantilever' or prob == 'michell':
        bc_x_max = geo_bc_tol
        bc_x_min = - geo_bc_tol
        bc_y_max = Ly + geo_bc_tol
        bc_y_min = - geo_bc_tol
        bc_z_max = h + geo_bc_tol
        bc_z_min = - geo_bc_tol

        ndof = 0
        for i in range(nnodes):
            if (bc_x_min <= X[i, 0] <= bc_x_max and
                bc_y_min <= X[i, 1] <= bc_y_max and
                bc_z_min <= X[i, 2] <= bc_z_max):
                # This is the bc node
                pass
            else:
                dof[i, 0] = ndof
                ndof += 1
                dof[i, 1] = ndof
                ndof += 1
                dof[i, 2] = ndof
                ndof += 1

    elif prob == 'MBB':
        bc1_x_max = geo_bc_tol
        bc1_x_min = -geo_bc_tol
        bc1_y_max = Ly + geo_bc_tol
        bc1_y_min = - geo_bc_tol
        bc1_z_max = h + geo_bc_tol
        bc1_z_min = - geo_bc_tol

        bc2_x_max = Lx + geo_bc_tol
        bc2_x_min = Lx*(1-MBB_bc_portion) - geo_bc_tol
        bc2_y_max = geo_bc_tol + geo_bc_tol
        bc2_y_min = - geo_bc_tol
        bc2_z_max = h + geo_bc_tol
        bc2_z_min = - geo_bc_tol

        ndof = 0
        for i in range(nnodes):
            if (bc1_x_min <= X[i, 0] <= bc1_x_max and
                bc1_y_min <= X[i, 1] <= bc1_y_max and
                bc1_z_min <= X[i, 2] <= bc1_z_max):
                # This is the bc node
                dof[i, 1] = ndof
                ndof += 1
            elif (bc2_x_min <= X[i, 0] <= bc2_x_max and
                  bc2_y_min <= X[i, 1] <= bc2_y_max and
                  bc2_z_min <= X[i, 2] <= bc2_z_max):
                # This is also bc node
                dof[i, 0] = ndof
                ndof += 1
            else:
                dof[i, 0] = ndof
                ndof += 1
                dof[i, 1] = ndof
                ndof += 1
                dof[i, 2] = ndof
                ndof += 1

    elif prob == 'lbracket':
        bc_x_max = Lx*ratio1 + geo_bc_tol
        bc_x_min = - geo_bc_tol
        bc_y_max = Ly + geo_bc_tol
        bc_y_min = Ly - geo_bc_tol
        bc_z_max = h + geo_bc_tol
        bc_z_min = - geo_bc_tol

        ndof = 0
        for i in range(nnodes):
            if (bc_x_min <= X[i, 0] <= bc_x_max and
                bc_y_min <= X[i, 1] <= bc_y_max and
                bc_z_min <= X[i, 2] <= bc_z_max):
                # This is the bc node
                pass
            else:
                dof[i, 0] = ndof
                ndof += 1
                dof[i, 1] = ndof
                ndof += 1
                dof[i, 2] = ndof
                ndof += 1

    # Set loading
    force = np.zeros(ndof)
    geo_tol = 1e-6

    if prob == 'cantilever':
        load_x_max = Lx + geo_tol
        load_x_min = Lx - geo_tol
        load_y_max = Ly*forced_portion + geo_tol
        load_y_min = - geo_tol
        load_z_max = 0.5*h*(1+loaded_thickness)
        load_z_min = 0.5*h*(1-loaded_thickness)

        force_dir = 1
        force_scal = -1.0

    elif prob == 'MBB':
        load_x_max = Lx*MBB_bc_portion + geo_tol
        load_x_min = - geo_tol
        load_y_max = Ly + geo_tol
        load_y_min = Ly - geo_tol
        load_z_max = 0.5*h*(1+loaded_thickness) + geo_tol
        load_z_min = 0.5*h*(1-loaded_thickness) - geo_tol

        force_dir = 1
        force_scal = -1.0

    elif prob == 'michell':
        load_x_max = Lx + geo_tol
        load_x_min = Lx - geo_tol
        load_y_max = 0.5*Ly*(1+forced_portion) + geo_tol
        load_y_min = 0.5*Ly*(1-forced_portion) - geo_tol
        load_z_max = 0.5*h*(1+loaded_thickness) + geo_tol
        load_z_min = 0.5*h*(1-loaded_thickness) - geo_tol

        force_dir = 1
        force_scal = -1.0

    elif prob == 'lbracket':
        load_x_max = Lx + geo_tol
        load_x_min = Lx - geo_tol
        load_y_max = Ly*ratio2 + geo_tol
        load_y_min = Ly*ratio2*(1-forced_portion) - geo_tol
        load_z_max = 0.5*h*(1+loaded_thickness) + geo_tol
        load_z_min = 0.5*h*(1-loaded_thickness) - geo_tol

        force_dir = 1
        force_scal = -1.0

    nforce = 0
    for i in range(nnodes):
        if (load_x_min <= X[i, 0] <= load_x_max and
            load_y_min <= X[i, 1] <= load_y_max and
            load_z_min <= X[i, 2] <= load_z_max):
            force[dof[i, force_dir]] = force_scal*force_magnitude
            nforce += 1

    force /= nforce

    # for n in range(X.shape[0]):
    #     print('node:{:4d} x:{:6f} y:{:6f} z:{:6f}'.format(n, X[n, 0], X[n, 1], X[n, 2]))

    # for ne in range(conn.shape[0]):
    #     print('elem:{:4d} conn = [{:3d} {:3d} {:3d} {:3d}]'.format(ne, conn[ne,0], conn[ne,1], conn[ne,2], conn[ne,3]))

    return nelems, nnodes, ndof, conn, X, dof, force


if __name__ == '__main__':
    n  = 20
    AR = 1.0
    prob = 'lbracket'
    ratio1 = 0.4
    ratio2 = 0.4
    hole_r = 0.25
    forced_portion = 0.3
    MBB_bc_portion = 0.3
    ly = 1.0
    force_magnitude = 25.0
    loaded_thickness = 0.5
    create_3d_unstruct_mesh(n, AR, prob, ratio1, ratio2, hole_r, forced_portion,
                            force_magnitude, MBB_bc_portion, loaded_thickness, ly)