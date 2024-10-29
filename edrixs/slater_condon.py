__all__ = ['get_gaunt', 'get_slater_condon', 'SlaterWF_radial']

import numpy as np
import math
from sympy.physics.wigner import gaunt
from .basis_transform import tmat_c2r, tmat_r2c, tmat_c2j, transform_utensor
from .utils import info_atomic_shell, case_to_shell_name, slater_integrals_name


def get_gaunt(l1, l2):
    """
    Calculate the Gaunt coefficents :math:`C_{l_1,l_2}(k,m_1,m_2)`

    .. math::

        C_{l_1,l_2}(k,m_1,m_2)=\\sqrt{\\frac{4\\pi}{2k+1}} \\int
        \\mathop{d\\phi} \\mathop{d\\theta} sin(\\theta)
        Y_{l_1}^{m_1\\star}(\\theta,\\phi) Y_{k}^{m_1-m_2}(\\theta,\\phi)
        Y_{l_2}^{m_2}(\\theta,\\phi)

    Parameters
    ----------
    l1: int
        The first quantum number of angular momentum.
    l2: int
        The second quantum number of angular momentum.

    Returns
    -------
    res: 3d float array
        The calculated Gaunt coefficents.

        The 1st index (:math:`= 0, 1, ..., l_1+l_2+1`) is the order :math:`k`.

        The 2nd index (:math:`= 0, 1, ... ,2l_1`) is the magnetic quantum
        number :math:`m_1` plus :math:`l_1`

        The 3nd index (:math:`= 0, 1, ... ,2l_2`) is the magnetic quantum
        number :math:`m_2` plus :math:`l_2`

    Notes
    -----
    It should be noted that :math:`C_{l_1,l_2}(k,m_1,m_2)` is
    nonvanishing only when

    :math:`k + l_1 + l_2 = \\text{even}`,

    and

    :math:`|l_1 -  l_2| \\leq k \\leq l_1 + l_2`.

    Please see Ref. [1]_ p. 10 for more details.

    References
    ----------
    .. [1] Sugano S, Tanabe Y and Kamimura H. 1970. Multiplets of
       Transition-Metal Ions in Crystals. Academic Press, New York and London.

    Examples
    --------
    >>> import edrixs

    Get gaunt coefficients between :math:`p`-shell and :math:`d`-shell

    >>> g = edrixs.get_gaunt(1, 2)

    """

    from sympy import N
    res = np.zeros((l1 + l2 + 1, 2 * l1 + 1, 2 * l2 + 1), dtype=np.float64)
    for k in range(l1 + l2 + 1):
        if not (np.mod(l1 + l2 + k, 2) == 0 and np.abs(l1 - l2) <= k <= l1 + l2):
            continue
        for i1, m1 in enumerate(range(-l1, l1 + 1)):
            for i2, m2 in enumerate(range(-l2, l2 + 1)):
                res[k, i1, i2] = (N(gaunt(l1, k, l2, -m1, m1 - m2, m2)) *
                                  (-1.0)**m1 * np.sqrt(4 * np.pi / (2 * k + 1)))
    return res

def SlaterWF_radial(r,Z,n,l):
        
        rho = 2.0 * r * Z / float(n)
        angular = ['s','p','d','f']
        orb = str(n) + angular[l]

        if orb == '1s':
            Rrad = 2.0 * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '2s':
            Rrad = 1.0 / (2.0 * math.sqrt(2.0)) * (2.0 - rho) * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '2p':
            Rrad = 1.0 / (2.0 * math.sqrt(6.0)) * rho * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '3s':
            Rrad = 1.0 / (9.0 * math.sqrt(3.0)) * (6.0 - 6.0 * rho + rho**2) * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '3p':
            Rrad = 1.0 / (9.0 * math.sqrt(6.0)) * rho * (4.0 - rho) * (Z**1.5) *np.exp(-rho/2.0)
        elif orb == '3d':
            Rrad = 1.0 / (9.0 * math.sqrt(30.0)) * (rho**2) * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '4s':
            Rrad = 1.0 / 96.0 * (24.0 - 36.0 * rho + 12.0 * (rho**2) - rho**3) * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '4p':
            Rrad = 1.0 / (32.0 * math.sqrt(15.0)) * rho * (20.0 - 10.0 * rho + rho**2) * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '4d':
            Rrad = 1.0 / (96.0 * math.sqrt(5.0)) * (rho**2) * (6.0 - rho) * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '4f':
            Rrad = 1.0 / (96.0 * math.sqrt(35.0)) * rho**3 * np.exp(-rho/2.0)
        elif orb == '5s':
            Rrad = 1.0 / (300.0 * math.sqrt(5.0)) * (120.0 - 240.0 * rho + 120.0 * (rho**2) - 20.0 * (rho**3) + rho**4) * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '5p':
            Rrad = 1.0 / (150.0 * math.sqrt(30.0)) * rho * (120.0 - 90.0 * rho + 18.0 * (rho**2) - rho**3) * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '5d':
            Rrad = 1.0 / (150.0 * math.sqrt(70.0)) * (rho**2) * (42.0 - 14.0 * rho + rho**2) * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '6s':
            Rrad = 1.0 / (2160.0 * math.sqrt(6.0)) * (720.0 - 1800.0 * rho + 1200.0 * (rho**2) - 300.0 * (rho**3) + 30.0 * (rho**4) - rho**5) * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '6p':
            Rrad = 1.0 / (432.0 * math.sqrt(210.0)) * rho * (840.0 - 840.0 * rho + 252.0 * (rho**2) - 28.0 * (rho**3) + (rho**4)) * (Z**1.5) * np.exp(-rho/2.0)
        elif orb == '6d':
            Rrad = 1.0 / (864.0 * math.sqrt(105.0)) * (rho**2) * (336.0 - 168.0 * rho + 24.0 * (rho**2) - (rho**3)) * (Z**1.5) * np.exp(-rho/2.0)
        else:
            Rrad = 0.0
            
        return Rrad


def get_slater_condon(n_arr, l_arr, Z_arr, k, r_arr):
    """
    Calculate the Slater-Condon coefficient :math:`F^{k}` using slater type orbitals

    Parameters
    ----------
    n_arr: numpy array of principle quantum numbers. 
    [n1,n2,n3,n4]
    l_arr: numpy array of orbital angular momentum quantum numbers. 
    [l1,l2,l3,l4]
    Z_arr: numpy array of effective charge. 
    [Z1,Z2,Z3,Z4]
    k : k to be computed with

    Returns
    -------
    Fk : a single integral value
    """

    R1 = SlaterWF_radial(r_arr, Z_arr[0], n_arr[0], l_arr[0])
    R2 = SlaterWF_radial(r_arr, Z_arr[1], n_arr[1], l_arr[1])
    R3 = SlaterWF_radial(r_arr, Z_arr[2], n_arr[2], l_arr[2])
    R4 = SlaterWF_radial(r_arr, Z_arr[3], n_arr[3], l_arr[3])

    dr = r_arr[1] - r_arr[0]
    nr = r_arr.shape[0]

    # [nr, nr']
    arr_1d = np.arange(nr)
    ir_arr = np.repeat(arr_1d[:,None], nr, axis=1)
    irp_arr = np.repeat(arr_1d[None,:], nr, axis=0)
    diff = ir_arr - irp_arr

    r_2d = np.repeat(r_arr[:,None], nr, axis=1)
    r_2dp = np.repeat(r_arr[None,:], nr, axis=0)

    kernel = np.where(diff>0, r_2dp**(k)/(r_2d**(k+1)), r_2d**(k)/(r_2dp**(k+1)))
    kernel[:nr, :nr] = np.diag((1.0 / r_arr))

    r_sq = r_arr**2
    Fk = np.einsum('r,p,r,p,rp,p,r->',r_sq,r_sq,R1,R2,kernel,R3,R4) * dr * dr
    
    return Fk