# Draw 2D NACA four digit airfoils
# Author: Brian J. Burke
# July 7, 2023

import numpy as np


class Naca4:
    def __init__(self, series: int):
        """
        Instantiate the quadratic basis function surrogate class.

        Parameters
        ----------
        M (int) : number of points in one dimension
        objFun (function handle) : pointer to 4objective function
        xlim (float) : limits of the domain
        """

        self.series = series
        self.maxCamber = self.series / 1000