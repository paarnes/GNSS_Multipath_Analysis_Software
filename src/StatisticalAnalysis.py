import numpy as np

class StatisticalAnalysis:
    """
    Class for compute statistical parameters
    """
    
    def __init__(self, A, l, N, h):
        """
        Initialize the class with:
        A : Design matrix (2D numpy array)
        l : Observations vector (1D numpy array)
        N : Normal equations matrix (2D numpy array)
        h : Adjustments vector (1D numpy array)
        """
        self.A = A
        self.l = l
        self.N = N
        self.h = h
        self.n, self.e = A.shape  # Number of observations and unknowns

    def compute_residuals(self):
        """Compute the residuals vector (V)."""
        V = self.A @ self.h - self.l
        return V

    def compute_sse(self, V):
        """Compute the Sum of Squared Errors (SSE)."""
        SSE = V.T @ V
        return SSE

    def compute_s0(self, SSE):
        """Compute the standard deviation of unit weight (S0)."""
        S0 = np.sqrt(SSE / (self.n - self.e))
        return S0

    def compute_covariance_matrix(self, Qxx, S0) -> np.ndarray:
        """Compute the Cxx covariance matrices based
        on cofactor matrix Qxx and standard deviation on unit weight"""
        Cxx = S0 ** 2 * Qxx
        return Cxx

    def compute_cofactor_matrix(self) -> np.ndarray:
        """Compute the Qxx and Cxx covariance matrices."""
        Qxx = np.linalg.inv(self.N)
        return Qxx

    def compute_cofactors(self, Qxx):
        """Compute cofactors qX, qY, qZ, qdT."""
        cof = np.diag(Qxx)
        qX, qY, qZ, qdT = cof[0], cof[1], cof[2], cof[3]
        return qX, qY, qZ, qdT

    def compute_dops(self, qX, qY, qZ, qdT):
        """Compute the DOPs (PDOP, TDOP, GDOP)."""
        PDOP = np.sqrt(qX + qY + qZ)
        TDOP = np.sqrt(qdT)
        GDOP = np.sqrt(PDOP ** 2 + TDOP ** 2)
        return PDOP, TDOP, GDOP

    def compute_standard_deviations(self, Cxx):
        """Compute standard deviations Sx, Sy, Sz, St."""
        std_devs = np.sqrt(np.diag(Cxx))
        Sx, Sy, Sz, St = std_devs[0], std_devs[1], std_devs[2], std_devs[3]
        return Sx, Sy, Sz, St

    def analyze(self):
        """
        Perform the complete analysis, computing all parameters, and return results.
        """
        # Step-by-step calculations
        V = self.compute_residuals()
        SSE = self.compute_sse(V)
        S0 = self.compute_s0(SSE)
        Qxx = self.compute_cofactor_matrix()
        Cxx = self.compute_covariance_matrix(Qxx, S0)
        qX, qY, qZ, qdT = self.compute_cofactors(Qxx)
        PDOP, TDOP, GDOP = self.compute_dops(qX, qY, qZ, qdT)
        Sx, Sy, Sz, St = self.compute_standard_deviations(Cxx)

        # Compile results into a dictionary
        return {
            "Residuals": V,
            "SSE": SSE,
            "S0": S0,
            "Qxx": Qxx,
            "Cxx": Cxx,
            "Cofactors": {"qX": qX, "qY": qY, "qZ": qZ, "qdT": qdT},
            "DOPs": {"PDOP": PDOP, "TDOP": TDOP, "GDOP": GDOP},
            "Standard Deviations": {"Sx": Sx, "Sy": Sy, "Sz": Sz, "St": St},
        }



