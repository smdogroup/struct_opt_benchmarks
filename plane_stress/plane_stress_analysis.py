
import plane_stress

class PlaneStressAnalysis:
    def __init__(self, rho, E, nu, ys, conn, bcs):

        # The variables associated with each element
        self.element_vars = None
        self.

        pass

    def applyFilter(self, x, xfilter):
        """
        Apply the filter to the design variables
        """
        pass

    def applyFilterTranspose(self, xfilter, x):
        """
        Apply the transpose operation to the filter
        """
        pass


    def assembleMat(self, x, mat_type='stiffness'):
        """
        Return an assembled stiffness, mass or geometric stiffness matrix
        with the specified matrix type.
        """

        pass

