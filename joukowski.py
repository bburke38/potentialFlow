# Potential flow modeling using Joukowski transformation from cylindrical flow.
# Author: Brian J. Burke

import numpy as np


def JtoZ(zeta: complex, c=1) -> complex:
    """
    Joukowski transformation from zeta (cylinder) to z (airfoil).

    Parameters
    ----------
    zeta (complex) : coordinates in the complex plane
    c (int = 1) : radius of the cylinder (semi-chord of airfoil)
    """
    z = zeta + c**2 / zeta

    return z

def JtoZeta(z: complex, c=1) -> complex:
    """
    Inverse Joukowski transformation from z (airfoil) to zeta (cylinder).

    Parameters
    ----------
    z (complex) : coordinates in the complex plane
    c (int = 1) : radius of the cylinder (semi-chord of airfoil)
    """
    zeta = z/2 * (1 + np.sqrt(1 - 4*c**2/(z**2)))

    return zeta

def calcStreamfunction(zeta, U=1.0, alpha=(np.pi*10/180), a=1.1, xi0=-0.1, eta0=0.1):
    """
    Calculate the streamfunctions of the flow around a cylinder at an angle of attack.
    Can be mapped back to airfoil coordinates.

    Parameters
    ----------
    zeta (complex) : coordinates in the complex plane
    U (float = 1.0) : freestream velocity parameter
    alpha : angle of attack in radians
    a : map circle to ellipse (a > c)
    xi0 : horizontal shift for symmetric airfoil
    eta0 : vertical shift for circular arc airfoil
    """
    zeta0 = xi0 + eta0*1j
    zeta = zeta - zeta0

    # Calculate streamfunction for flow around a cylinder
    F = U * ( np.exp(-alpha*1j)*zeta + np.exp(alpha*1j)*a**2/zeta )

    # Add correction for circulation to move a stagnation point to the trailing edge
    gamma = 4*np.pi*U*a*np.sin(alpha)
    F = F + gamma*1j/(2*np.pi) * np.log(zeta/a)

    return np.imag(F)


def F_circle(x):
    z = x[0] + x[1]*1j

    F = calcStreamfunction(z)

    return F

def calcPsiFromZeta(xy, U=1.0, alpha=(np.pi*10/180), a=1.1, xi0=-0.1, eta0=0.1):
    """
    Calculate stream function for a real-valued set of input coordinates.

    Parameters
    ----------
    xy: set of coordinates in real space
    """
    z = xy[0] + xy[1]*1j

    zeta = JtoZeta(z)

    F = calcStreamfunction(zeta, U, alpha, a, xi0, eta0)

    return F
class JFoil:
    def __init__(self, a: int, xi0=0, eta0=0):

        self.a = a
        self.xi0 = xi0
        self.eta0 = eta0

        return
    
    def gen_circle(self):
        nu = np.linspace(0, 2*np.pi, 100)
        nui = nu * 1j
        zeta = self.a * np.exp(nui) + self.xi0 + self.eta0 * 1j

        return zeta
    
    def toReal(self):
        zeta = self.gen_circle()
        z = zeta + 1 / zeta

        return z
    
