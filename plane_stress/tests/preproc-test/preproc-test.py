from plane_stress.preprocessors import preproc

n = 48
AR = 1.5
ratio1 = 0.4
ratio2 = 0.4
hole_radius = 0.25
nr0 = 32
outdir = 'pkls'
plot_mesh = True

force_magnitude = 25.0  # total force applied to structure
forced_portion = 0.2  # portion of the edge of which load is applied to
MBB_bc_portion = 0.1
ly = 1.0  # for this project we always set height of domain to be 1.0
density = 2700.0  # this is only used for frequency analysis
E = 70e3  # Young's modulus
nu = 0.3  # Poisson's ratio

for prob in ['cantilever', 'michell', 'MBB', 'lbracket']:
    for meshtype in ['structured', 'unstructured']:
        for use_hole in [False, True]:
            for use_concentrated_force in [False, True]:
                if use_hole is True and meshtype == 'structured':
                    pass
                else:
                    preproc(n, AR, prob, meshtype, ratio1, ratio2,
                        use_concentrated_force, use_hole, hole_radius,
                        nr0, outdir, plot_mesh, force_magnitude,
                        forced_portion, MBB_bc_portion, ly, density, E, nu)


