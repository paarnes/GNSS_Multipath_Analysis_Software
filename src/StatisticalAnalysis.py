import numpy as np

class StatisticalAnalysis:
    """
    Class for performing statistical analysis on observational and computational data.

    This class computes key statistical parameters such as residuals, sum of squared errors (SSE),
    standard deviation of unit weight (S0), cofactor and covariance matrices, Dilution of Precision (DOPs),
    and standard deviations of parameters. These metrics are commonly used in geodetic and positioning computations.

    Attributes:
    ----------
    A : np.ndarray
        Design matrix (shape: [n_observations, n_parameters]).
    l : np.ndarray
        Observations vector (shape: [n_observations]).
    N : np.ndarray
        Normal equations matrix (shape: [n_parameters, n_parameters]).
    h : np.ndarray
        Adjustments vector (shape: [n_parameters]).
    n : int
        Number of observations (rows in `A`).
    e : int
        Number of unknowns (columns in `A`).
    """

    def __init__(self, A: np.ndarray, l:np.ndarray, N:np.ndarray, h:np.ndarray):
        """
        Initialize the StatisticalAnalysis class.

        Parameters:
        ----------
        A : np.ndarray
            Design matrix.
        l : np.ndarray
            Observations vector.
        N : np.ndarray
            Normal equations matrix.
        h : np.ndarray
            Adjustments vector.
        """
        self.A = A
        self.l = l
        self.N = N
        self.h = h
        self.n, self.e = A.shape  # Number of observations and unknowns

    def compute_residuals(self):
        """
        Compute the residuals vector (V), representing the differences
        between observed and computed values.

        Returns:
        -------
        np.ndarray
            Residuals vector (shape: [n_observations]).
        """
        V = self.A @ self.h - self.l
        return V

    def compute_sse(self, V):
        """
        Compute the Sum of Squared Errors (SSE).

        Parameters:
        ----------
        V : np.ndarray
            Residuals vector.

        Returns:
        -------
        float
            Sum of squared errors.
        """
        SSE = V.T @ V
        return SSE

    def compute_s0(self, SSE):
        """
        Compute the standard deviation of unit weight (S0).

        Parameters:
        ----------
        SSE : float
            Sum of squared errors.

        Returns:
        -------
        float
            Standard deviation of unit weight.
        """
        S0 = np.sqrt(SSE / (self.n - self.e))
        return S0

    def compute_covariance_matrix(self, Qxx, S0) -> np.ndarray:
        """
        Compute the covariance matrix (Cxx) based on the cofactor matrix
        (Qxx) and standard deviation of unit weight (S0).

        Parameters:
        ----------
        Qxx : np.ndarray
            Cofactor matrix.
        S0 : float
            Standard deviation of unit weight.

        Returns:
        -------
        np.ndarray
            Covariance matrix (shape: [n_parameters, n_parameters]).
        """
        Cxx = S0 ** 2 * Qxx
        return Cxx

    def compute_cofactor_matrix(self) -> np.ndarray:
        """
        Compute the cofactor matrix (Qxx).

        Returns:
        -------
        np.ndarray
            Cofactor matrix (shape: [n_parameters, n_parameters]).
        """
        Qxx = np.linalg.inv(self.N)
        return Qxx

    def compute_cofactors(self, Qxx):
        """
        Extract cofactors for specific parameters: qX, qY, qZ, qdT.

        Parameters:
        ----------
        Qxx : np.ndarray
            Cofactor matrix.

        Returns:
        -------
        Tuple[float, float, float, float]
            Cofactors corresponding to X, Y, Z, and dT.
        """
        cof = np.diag(Qxx)
        qX, qY, qZ, qdT = cof[0], cof[1], cof[2], cof[3]
        return qX, qY, qZ, qdT

    def compute_dops(self, qX, qY, qZ, qdT):
        """
        Compute the Dilution of Precision (DOP) metrics: PDOP, TDOP, and GDOP.

        PDOP: Position (3D) dilution of precision
        TDOP: Time dilution of precision
        GDOP: Geometric dilution of precision

        Parameters:
        ----------
        qX : float
            Cofactor for X coordinate.
        qY : float
            Cofactor for Y coordinate.
        qZ : float
            Cofactor for Z coordinate.
        qdT : float
            Cofactor for clock bias.

        Returns:
        -------
        Tuple[float, float, float]
            Positional DOP (PDOP), Time DOP (TDOP), and Geometric DOP (GDOP).
        """
        PDOP = np.sqrt(qX + qY + qZ)
        TDOP = np.sqrt(qdT)
        GDOP = np.sqrt(PDOP ** 2 + TDOP ** 2)
        return PDOP, TDOP, GDOP

    def compute_standard_deviations(self, Cxx):
        """
        Compute standard deviations for the parameters: Sx, Sy, Sz, St.

        Parameters:
        ----------
        Cxx : np.ndarray
            Covariance matrix.

        Returns:
        -------
        Tuple[float, float, float, float]
            Standard deviations for X, Y, Z, and clock bias.
        """
        std_devs = np.sqrt(np.diag(Cxx))
        Sx, Sy, Sz, St = std_devs[0], std_devs[1], std_devs[2], std_devs[3]
        return Sx, Sy, Sz, St

    def run_statistical_analysis(self, n_decimals:int=3):
        """
        Perform a complete statistical analysis, compute all parameters,
        and return the results in a dictionary.

        Parameters:
        ----------
        n_decimals : int, optional
            Number of decimal places for rounding the results (default is 3).

        Returns:
        -------
        dict
            Dictionary containing the computed residuals, SSE, S0, Qxx, Cxx,
            cofactors, DOPs, and standard deviations.
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

        # Round all results to the specified number of decimals
        V = np.round(V, n_decimals)
        SSE = np.round(SSE, n_decimals)
        S0 = np.round(S0, n_decimals)
        Qxx = np.round(Qxx, n_decimals)
        Cxx = np.round(Cxx, n_decimals)
        qX, qY, qZ, qdT = map(lambda x: np.round(x, n_decimals), (qX, qY, qZ, qdT))
        PDOP, TDOP, GDOP = map(lambda x: np.round(x, n_decimals), (PDOP, TDOP, GDOP))
        Sx, Sy, Sz, St = map(lambda x: np.round(x, n_decimals), (Sx, Sy, Sz, St))

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



