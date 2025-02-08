#!/usr/bin/env python

__all__ = ['combination', 'fock_bin', 'get_fock_bin_by_N', 'get_fock_half_N',
           'get_fock_full_N', 'get_fock_basis_by_NLz', 'get_fock_basis_by_NSz',
           'get_fock_basis_by_NJz', 'get_fock_basis_by_N_abelian',
           'get_fock_basis_by_N_LzSz', 'write_fock_dec_by_N', 'write_fock_dec_by_N_constrainedN1N2',
           'write_fock_dec_by_N_constrainedN1N2_multi', 'write_fock_dec_by_N_constrainedN1N2N3',
           'write_fock_dec_by_N_general']

import numpy as np
import itertools


def combination(n, m):
    """
    Calculate the combination :math:`C_{n}^{m}`,

    .. math::

        C_{n}^{m} = \\frac{n!}{m!(n-m)!}.

    Parameters
    ----------
    n: int
       Number n.
    m: int
        Number m.

    Returns
    -------
    res: int
        The calculated result.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.combination(6, 2)
    15

    """

    if m > n or n < 0 or m < 0:
        print("wrong number in combination")
        return
    if m == 0 or n == m:
        return 1

    largest = max(m, n - m)
    smallest = min(m, n - m)
    numer = 1.0
    for i in range(largest + 1, n + 1):
        numer *= i

    denom = 1.0
    for i in range(1, smallest + 1):
        denom *= i

    res = int(numer / denom)
    return res


def fock_bin(n, k):
    """
    Return all the possible :math:`n`-length binary
    where :math:`k` of :math:`n` digitals are set to 1.

    Parameters
    ----------
    n: int
        Binary length :math:`n`.
    k: int
        How many digitals are set to be 1.

    Returns
    -------
    res: list of int-lists
        A list of list containing the binary digitals.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.fock_bin(4, 2)
    [[1, 1, 0, 0],
     [1, 0, 1, 0],
     [1, 0, 0, 1],
     [0, 1, 1, 0],
     [0, 1, 0, 1],
     [0, 0, 1, 1]]

    """

    if n == 0:
        return [[0]]

    res = []
    for bits in itertools.combinations(list(range(n)), k):
        s = [0] * n
        for bit in bits:
            s[bit] = 1
        res.append(s)
    return res


def get_fock_bin_by_N(*args):
    """
    Get binary form to represent a Fock state.

    Parameters
    ----------
    args: ints
        args[0]: number of orbitals for 1st-shell,

        args[1]: number of occupancy for 1st-shell,

        args[2]: number of orbitals for 2nd-shell,

        args[3]: number of occupancy for 2nd-shell,

        ...

        args[ :math:`2N-2`]: number of orbitals for :math:`N` th-shell,

        args[ :math:`2N-1`]: number of occupancy for :math:`N` th-shell.

    Returns
    -------
    result: list of int list
        The binary form of Fock states.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.get_fock_bin_by_N(4, 2)
    [[1, 1, 0, 0],
     [1, 0, 1, 0],
     [1, 0, 0, 1],
     [0, 1, 1, 0],
     [0, 1, 0, 1],
     [0, 0, 1, 1]]

    >>> edrixs.get_fock_bin_by_N(4, 2, 2, 1)
    [[1, 1, 0, 0, 1, 0],
     [1, 0, 1, 0, 1, 0],
     [1, 0, 0, 1, 1, 0],
     [0, 1, 1, 0, 1, 0],
     [0, 1, 0, 1, 1, 0],
     [0, 0, 1, 1, 1, 0],
     [1, 1, 0, 0, 0, 1],
     [1, 0, 1, 0, 0, 1],
     [1, 0, 0, 1, 0, 1],
     [0, 1, 1, 0, 0, 1],
     [0, 1, 0, 1, 0, 1],
     [0, 0, 1, 1, 0, 1]]

    """

    n = len(args)

    if n % 2 != 0:
        print("Error: number of arguments is not even")
        return

    if n == 2:
        return fock_bin(args[0], args[1])
    else:
        result = []
        res1 = fock_bin(args[0], args[1])
        res2 = get_fock_bin_by_N(*args[2:])
        for ifock in res2:
            for jfock in res1:
                result.append(jfock + ifock)
        return result


def get_fock_half_N(N):
    res = [[] for i in range(N + 1)]
    for i in range(2**N):
        occu = bin(i).count('1')
        res[occu].append(i)
    return res


def get_fock_full_N(norb, N):
    """
    Get the decimal digitals to represent Fock states.

    Parameters
    ----------
    norb: int
        Number of orbitals.
    N: int
        Number of occupancy.

    Returns
    -------
    res: list of int
        The decimal digitals to represent Fock states.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.fock_bin(4,2)
    [[1, 1, 0, 0],
     [1, 0, 1, 0],
     [0, 1, 1, 0],
     [1, 0, 0, 1],
     [0, 1, 0, 1],
     [0, 0, 1, 1]]

    >>> import edrixs
    >>> edrixs.get_fock_full_N(4,2)
    [3, 5, 6, 9, 10, 12]

    """

    res = []
    half_N = get_fock_half_N(norb // 2)
    for m in range(norb // 2 + 1):
        n = N - m
        if n >= 0 and n <= norb // 2:
            res.extend([i * 2**(norb // 2) + j for i in half_N[m] for j in half_N[n]])
    return res


def get_fock_basis_by_NLz(norb, N, lz_list):
    """
    Get decimal digitals to represent Fock states, use good quantum number:

    - orbital angular momentum :math:`L_{z}`

    Parameters
    ----------
    norb: int
        Number of orbitals.
    N: int
        Number of total occupancy.
    lz_list: list of int
        Quantum number :math:`l_{z}` for each orbital.

    Returns
    -------
    res: dict
        A dictionary containing the decimal digitals, the key is good
        quantum numbers :math:`L_{z}`, the value is a list of int.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.get_fock_basis_by_NLz(6, 2, [-1, -1, 0, 0, 1, 1])
    {
     -2: [3],
     -1: [5, 6, 9, 10],
      0: [12, 17, 18, 33, 34],
      1: [20, 36, 24, 40],
      2: [48]
    }
    """

    res = get_fock_basis_by_N_abelian(norb, N, lz_list)
    return res


def get_fock_basis_by_NSz(norb, N, sz_list):
    """
    Get decimal digitals to represent Fock states, use good quantum number:

    - spin angular momentum :math:`S_{z}`

    Parameters
    ----------
    norb: int
        Number of orbitals.
    N: int
        Number of total occupancy.
    sz_list: list of int
        Quantum number :math:`s_{z}` for each orbital.

    Returns
    -------
    res: dict
        A dictionary containing the decimal digitals, the key is good quantum
        numbers :math:`S_{z}`, the value is a list of int.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.get_fock_basis_by_NSz(6, 2, [1, -1, 1, -1, 1, -1])
    {
     -2: [10, 34, 40],
     -1: [],
      0: [3, 6, 9, 12, 18, 33, 36, 24, 48],
      1: [],
      2: [5, 17, 20]
    }
    """

    res = get_fock_basis_by_N_abelian(norb, N, sz_list)
    return res


def get_fock_basis_by_NJz(norb, N, jz_list):
    """
    Get decimal digitals to represent Fock states, use good quantum number:

    - total angular momentum :math:`J_{z}`

    Parameters
    ----------
    norb: int
        Number of orbitals.
    N: int
        Number of total occupancy.
    jz_list: list of int
        Quantum number :math:`j_{z}` for each orbital.

    Returns
    -------
    res: dict
        A dictionary containing the decimal digitals, the key is good quantum
        numbers :math:`j_{z}`, the value is a list of int.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.get_fock_basis_by_NJz(6, 2, [-1, 1, -3, -1, 1, 3])
    {
     -6: [],
     -5: [],
     -4: [5, 12],
     -3: [],
     -2: [6, 9, 20],
     -1: [],
      0: [3, 10, 17, 36, 24],
      1: [],
      2: [18, 33, 40],
      3: [],
      4: [34, 48],
      5: [],
      6: []
    }
    """

    res = get_fock_basis_by_N_abelian(norb, N, jz_list)
    return res


def get_fock_basis_by_N_abelian(norb, N, a_list):
    """
    Get decimal digitals to represent Fock states, use some Abelian good quantum number.

    Parameters
    ----------
    norb: int
        Number of orbitals.
    N: int
        Number of total occupancy.
    a_list: list of int
        Quantum number of the Abelian symmetry for each orbital.

    Returns
    -------
    basis: dict
        A dictionary containing the decimal digitals, the key is good quantum numbers,
        the value is a list of int.
    """

    result = get_fock_full_N(norb, N)
    min_a, max_a = min(a_list) * N, max(a_list) * N
    basis = {}
    for i in range(min_a, max_a + 1):
        basis[i] = []
    for n in result:
        a = sum([a_list[i] for i in range(0, n.bit_length()) if (n >> i & 1)])
        basis[a].append(n)
    return basis


def get_fock_basis_by_N_LzSz(norb, N, lz_list, sz_list):
    """
    Get decimal digitals to represent Fock states, use good quantum number:

    - orbital angular momentum :math:`L_{z}`
    - spin angular momentum :math:`S_{z}`

    Parameters
    ----------
    norb: int
        Number of orbitals.
    N: int
        Number of total occupancy.
    lz_list: list of int
        Quantum number :math:`l_{z}` for each orbital.
    sz_list: list of int
        Quantum number :math:`s_{z}` for each orbital.

    Returns
    -------
    basis: dict
        A dictionary containing the decimal digitals, the key is a tuple containing good quantum
        numbers ( :math:`l_{z}`, :math:`s_{z}`), the value is a list of int.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.get_fock_basis_by_N_LzSz(6, 2, [-1, -1, 0, 0, 1, 1], [1, -1, 1, -1, 1, -1])
    {
     (-2, -2): [],
     (-2, -1): [],
      (-2, 0): [3],
      (-2, 1): [],
      (-2, 2): [],
     (-1, -2): [10],
     (-1, -1): [],
      (-1, 0): [6, 9],
      (-1, 1): [],
      (-1, 2): [5],
      (0, -2): [34],
      (0, -1): [],
       (0, 0): [12, 18, 33],
       (0, 1): [],
       (0, 2): [17],
      (1, -2): [40],
      (1, -1): [],
       (1, 0): [36, 24],
       (1, 1): [],
       (1, 2): [20],
      (2, -2): [],
      (2, -1): [],
       (2, 0): [48],
       (2, 1): [],
       (2, 2): []
    }
    """
    result = get_fock_full_N(norb, N)
    min_Lz, max_Lz = min(lz_list) * N, max(lz_list) * N
    min_Sz, max_Sz = min(sz_list) * N, max(sz_list) * N
    basis = {}
    for i in range(min_Lz, max_Lz + 1):
        for j in range(min_Sz, max_Sz + 1):
            basis[(i, j)] = []
    for n in result:
        Lz, Sz = np.sum([[lz_list[i], sz_list[i]] for i in range(0, n.bit_length())
                         if (n >> i & 1)], axis=0)
        basis[(Lz, Sz)].append(n)
    return basis


def write_fock_dec_by_N(N, r, fname='fock_i.in'):
    """
    Get decimal digitals to represent Fock states, sort them by
    ascending order and then write them to file.

    Parameters
    ----------
    N: int
       Number of orbitals.
    r: int
        Number of occuancy.
    fname: string
        File name.

    Returns
    -------
    ndim: int
        The dimension of the Hilbert space

    Examples
    --------
    >>> import edrixs
    >>> edrixs.write_fock_dec_by_N(4, 2, 'fock_i.in')
    file fock_i.in looks like
    15
    3
    5
    6
    9
    10
    12
    17
    18
    20
    24
    33
    34
    36
    40
    48

    where, the first line is the total numer of Fock states,
    and the following lines are the Fock states in decimal form.
    """

    res = get_fock_full_N(N, r)
    res.sort()
    ndim = len(res)
    f = open(fname, 'w')
    print(ndim, file=f)
    for item in res:
        print(item, file=f)
    f.close()
    return ndim

def product_extend_general(res_all, N_all):
    """
    Combined the Fock states with "1" & "2"!
    order : 1 --> 2 --> 3 --> 4 ....
    
    Parameters:
    ------------
    res_all [res1_all, res2_all, res3_all...]
    N_all [N1, N2, N3]
    
    res1_all [#_vnoccu, #1_combinations]
    res2_all [#_vnoccu, #2_combinations]
    ...
    N1
    N2
    N3
    ...
    
    # output :
    res_all [#_vnocci x #1 x #2]
    """
    n_case = len(res_all)
    nr = len(res_all[0])
    
    # Summing over all cases to prevent memory-demanding high-dimensional array
    for icase in range(n_case-1):
        if icase == 0:
            res_now = res_all[icase]
            N_now = N_all[icase]
        else:
            res_now = res_cur
            N_now = np.sum(N_all[:icase+1])

        res_next = res_all[icase+1]
        N_next = N_all[icase+1]

        print(icase, len(res_now), len(res_next))

        scale_ = int(2**N_next)
        res_cur = []  #[nr]
        for ir in range(nr):
            res_int = []
            ndim_now = len(res_now[ir])
            ndim_next = len(res_next[ir])
            for i in range(ndim_now):
                for j in range(ndim_next):
                    res = int(res_now[ir][i]*scale_ + res_next[ir][j])
                    res_int.append(res)
            res_cur.append(res_int)

    #res_all = res_cur.flatten()    
    res_all = [x for xs in res_cur for x in xs]
    return res_all

def product_extend(res1_all, res2_all, N1, N2):
    """
    Combined the Fock states with "1" & "2"!
    order : 1 --> 2
    
    Parameters:
    ------------
    res1_all [#_vnoccu, #1_combinations]
    res2_all [#_vnoccu, #2_combinations]
    N1
    N2
    
    # output :
    res_all [#_vnocci x #1 x #2]
    """
    nr = len(res1_all)
    scale_ = int(2**N2)
    res_all = []
    
    for ir in range(nr):
        ndim1 = len(res1_all[ir])
        ndim2 = len(res2_all[ir])
        for i in range(ndim1):
            for j in range(ndim2):
                res = int(res1_all[ir][i]*scale_ + res2_all[ir][j])
                res_all.append(res)

    return res_all

def product_extend2(res1_all, res2_all, N1, N2):
    """
    Combined the Fock states with "1" & "2"!
    order : 1 --> 2
    
    Parameters:
    ------------
    res1_all [#_vnoccu, #1_combinations]
    res2_all [#_vnoccu, #2_combinations]
    N1
    N2
    
    # output :
    res_all [#_vnocci, #1 x #2]
    """
    nr = len(res1_all)
    scale_ = int(2**N2)
    res_all = []
    
    for ir in range(nr):
        ndim1 = len(res1_all[ir])
        ndim2 = len(res2_all[ir])
        res_r = []
        for i in range(ndim1):
            for j in range(ndim2):
                res = int(res1_all[ir][i]*scale_ + res2_all[ir][j])
                res_r.append(res)
        res_all.append(res_r)

    return res_all

def product_extend_multisites(res_imp_all, Nimp):
    """
    Combined the Fock states with "1" -> "2"!
    
    Parameters:
    ------------
    res_imp_all     : [Nsites, #_noccu, #_combinations*]
    Nimp            : [Nsites]

    Output:
    res_all         : [#_occu, #_combinations*]
    """
    nsites = len(res_imp_all)
    nr = len(res_imp_all[0])
    
    N_sum = np.sum(np.array(Nimp))
    N_sum_arr = np.array([ np.sum(np.array(Nimp[i:])) for i in range(nsites) ])
    scale_arr = 2**N_sum / np.power(2, N_sum_arr)   # [nsites]
    
    res_all = []
    
    cur = nsites-1
    res_now = res_imp_all[cur]
    N_now = Nimp[cur]
    for isite in range(nsites-1):
        res_comb = product_extend2(res_imp_all[cur-1], res_now,\
                Nimp[cur-1], N_now)
        
        cur = cur - 1
        N_now = N_now + Nimp[cur]
        res_now = res_comb
    
    #print(len(res_now), len(res_now[1]))

    #for ir in range(nr):
    #    for isite in range(nsites):     # Sum over sites
    #        if isite == 0:
    #            comb = np.array(res_imp_all[isite][ir])*scale_arr[isite]  # [#_comb*]
    #        else:
    #            comb = comb +\
    #                    np.array(res_imp_all[isite][ir])*scale_arr[isite]  # [#_comb*]
    #            
    #    res_all.append(comb)

    return res_now, N_sum

def write_fock_dec_by_N_constrainedN1N2N3(N1, r1_range, N2, r2_range, N3, r3_range, fname='fock_i.in'):
    """
    Get decimal digitals to represent Fock states, sort them by
    ascending order and then write them to file.

    r1_range / r2_range / r3_range should have the same number of elements.

    But now we want to constrain the occupation of Fe / Oxygen to some certain number
    of electrons...! "1" is moved to the left and larger than "2".

    Parameters
    ----------
    N: int
       Number of orbitals.
    r: int
        Number of occuancy.
    fname: string
        File name.

    Returns
    -------
    ndim: int
        The dimension of the Hilbert space

    Examples
    --------
    >>> import edrixs
    >>> edrixs.write_fock_dec_by_N(4, 2, 'fock_i.in')
    file fock_i.in looks like
    15
    3
    5
    6
    9
    10
    12
    17
    18
    20
    24
    33
    34
    36
    40
    48

    where, the first line is the total numer of Fock states,
    and the following lines are the Fock states in decimal form.
    """

    # Input output order : N1 -> N2
    # fock order : N2 -> N1

    # Bugs still exist in for the 3rd case here...!

    # N1, r1
    nr1 = int(r1_range.shape[0])
    res1_all = []
    for ir1, r1 in enumerate(r1_range):
        res1_ = get_fock_full_N(N1, r1)
        res1_all.append(res1_)

    # N2, r2
    nr2 = int(r2_range.shape[0])
    res2_all = []
    for ir2, r2 in enumerate(r2_range):
        res2_ = get_fock_full_N(N2, r2)
        res2_all.append(res2_)

    # N3, r3
    nr3 = int(r3_range.shape[0])
    res3_all = []
    for ir3, r3 in enumerate(r3_range):
        res3_ = get_fock_full_N(N3, r3)
        res3_all.append(res3_)
    
    # Fock basis writing order : [3 -> 2 -> 1]
    res_combined_32 = product_extend(res3_all, res2_all, N3, N2)
    res_combined_32.sort()

    res_combined = product_extend(res_combined_32, res1_all, N3+N2, N1)
    res_combined.sort()

    ndim = len(res_combined)
    f = open(fname, 'w')
    print(ndim, file=f)
    for item in res_combined:
        print(item, file=f)
    f.close()
    return ndim

def write_fock_dec_by_N_constrainedN1N2(N1, r1_range, N2, r2_range, fname='fock_i.in'):
    """
    Get decimal digitals to represent Fock states, sort them by
    ascending order and then write them to file.

    r1_range / r2_range should have the same number of elements.

    But now we want to constrain the occupation of Fe / Oxygen to some certain number
    of electrons...! "1" is moved to the left and larger than "2".

    Parameters
    ----------
    N: int
       Number of orbitals.
    r: int
        Number of occuancy.
    fname: string
        File name.

    Returns
    -------
    ndim: int
        The dimension of the Hilbert space

    Examples
    --------
    >>> import edrixs
    >>> edrixs.write_fock_dec_by_N(4, 2, 'fock_i.in')
    file fock_i.in looks like
    15
    3
    5
    6
    9
    10
    12
    17
    18
    20
    24
    33
    34
    36
    40
    48

    where, the first line is the total numer of Fock states,
    and the following lines are the Fock states in decimal form.
    """

    # Input output order : N1 -> N2
    # fock order : N2 -> N1

    # N1, r1
    nr1 = int(r1_range.shape[0])
    res1_all = []
    for ir1, r1 in enumerate(r1_range):
        res1_ = get_fock_full_N(N1, r1)
        res1_all.append(res1_)

    # N2, r2
    nr2 = int(r2_range.shape[0])
    res2_all = []
    for ir2, r2 in enumerate(r2_range):
        res2_ = get_fock_full_N(N2, r2)
        res2_all.append(res2_)
    
    # Fock basis writing order : [2 -> 1]
    res_combined = product_extend(res2_all, res1_all, N2, N1)
    res_combined.sort()
    ndim = len(res_combined)
    f = open(fname, 'w')
    print(ndim, file=f)
    for item in res_combined:
        print(item, file=f)
    f.close()
    return ndim

def write_fock_dec_by_N_constrainedN1N2_multi(N_imp, rimp_range, N_baths, rbaths_range, fname='fock_i.in'):
    """
    Get decimal digitals to represent Fock states, sort them by
    ascending order and then write them to file.

    r1_range / r2_range should have the same number of elements.

    But now we want to constrain the occupation of Fe / Oxygen to some certain number
    of electrons...! "1" is moved to the left and larger than "2".

    Parameters
    ----------
    N_imp[nsites]: int
       Number of orbitals.
    rimp_range [nsites, possilbe occupation numbers] : int
        Number of occuancy.
    N_baths: int
       Number of orbitals.
    rbaths_range [possible occupation number] : int
        Number of occuancy.
    fname: string
        File name.

    Returns
    -------
    ndim: int
        The dimension of the Hilbert space

    Examples
    --------
    >>> import edrixs
    >>> edrixs.write_fock_dec_by_N(4, 2, 'fock_i.in')
    file fock_i.in looks like
    15
    3
    5
    6
    9
    10
    12
    17
    18
    20
    24
    33
    34
    36
    40
    48

    where, the first line is the total numer of Fock states,
    and the following lines are the Fock states in decimal form.
    """

    # Input output order : N1 -> N2
    # fock order : N2 -> N1

    print("Building constrained basis...!")
    # N1, r1 -> Impurity
    nsites = int(N_imp.shape[0])
    resimp_all_sites = []

    for isite in range(nsites):
        print("Site ", isite)
        nr_imp = int(rimp_range.shape[1])
        resimp_all = []
        for irimp, rimp in enumerate(rimp_range[isite,:]):
            print("rimp", rimp)
            res_imp = get_fock_full_N(N_imp[isite], rimp)       # Array of 
            resimp_all.append(res_imp)

        resimp_all_sites.append(resimp_all)
    
    #print(resimp_all)
    # resimp_all_sites [nsite, noccu, n_combination*]

    # N2, r2
    nr_baths = int(rbaths_range.shape[0])
    resbaths_all = []
    for irbaths, rbaths in enumerate(rbaths_range):
        print("baths :", irbaths)
        print("N, r", N_baths, rbaths)
        res_baths = get_fock_full_N(N_baths, rbaths)
        resbaths_all.append(res_baths)

    # resbaths_all[noccu, n_combination*]

    # Impurity - Combining several impurities
    #last = nsites - 1
    #res_now = resimp_all_sites[last]        # still 2D
    #N_now = N_imp[last]
    
    resimp_all_comb, N_imps = product_extend_multisites(resimp_all_sites, N_imp)
    #print(N_imps)
    #print(len(resimp_all_comb), len(resimp_all_comb[1]))
    #print(len(resbaths_all), len(resbaths_all[1]))

    #for isite in range(nsites-1):
    #    print("combining site ", isite)
    #    res_combined = product_extend(resimp_all_sites[last-1], res_now, N_imp[last-1], N_now)
    #    last = last - 1
    #    res_now = res_combined
    #    N_now = N_now + N_imp[last-1]
    #
    #print(len(res_now))

    # Baths
    print("combining baths!")
    res_combined = product_extend(resbaths_all, resimp_all_comb, N_baths, N_imps)

    print(len(res_combined))

    res_combined.sort()
    ndim = len(res_combined)
    f = open(fname, 'w')
    print(ndim, file=f)
    for item in res_combined:
        print(item, file=f)
    f.close()
    return ndim

def write_fock_dec_by_N_general(N_imp, rimp_range, N_bath, rbath_range, fname='fock_i.in'):
    """
    Get decimal digitals to represent Fock states, sort them by
    ascending order and then write them to file.

    r1_range / r2_range should have the same number of elements.

    General case for multi-impurities and multi-baths.

    Parameters
    ----------
    N_imp [nimp]: int
       Number of orbitals.
    rimp_range [nimp, possilbe occupation numbers] : int
        Number of occuancy.
    N_bath [nbaths]: int
       Number of orbitals.
    rbath_range [nbaths, possible occupation number] : int
        Number of occuancy.
    fname: string
        File name.

    Returns
    -------
    ndim: int
        The dimension of the Hilbert space

    Examples
    --------
    >>> import edrixs
    >>> edrixs.write_fock_dec_by_N(4, 2, 'fock_i.in')
    file fock_i.in looks like
    15
    3
    5
    6
    9
    10
    12
    17
    18
    20
    24
    33
    34
    36
    40
    48

    where, the first line is the total numer of Fock states,
    and the following lines are the Fock states in decimal form.
    """

    # Input output order : N1 -> N2
    # fock order : N2 -> N1

    res_all = []
    N_all = []

    print("Building constrained basis...!")
    # N1, r1 -> Impurity
    nsites = int(N_imp.shape[0])
    resimp_all_sites = []

    for isite in range(nsites):
        print("Site ", isite)
        nr_imp = int(rimp_range.shape[1])
        resimp_all = []
        for irimp, rimp in enumerate(rimp_range[isite,:]):
            print("rimp", rimp)
            res_imp = get_fock_full_N(N_imp[isite], rimp)       # Array of 
            resimp_all.append(res_imp)

        resimp_all_sites.append(resimp_all)
        N_all.append(N_imp[isite])
        res_all.append(resimp_all)

    # N2, r2 -> baths
    nbaths = int(N_bath.shape[0]) 
    resbaths_all_sites = []

    for ibath in range(nbaths):
        print("Baths ", ibath)
        nr_bath = int(rbath_range.shape[1])
        resbaths_all = []
        for irbath, rbath in enumerate(rbath_range[ibath,:]):
            print("rbath", rbath)
            res_bath = get_fock_full_N(N_bath[ibath], rbath)
            resbaths_all.append(res_bath)

        resbaths_all_sites.append(resbaths_all)
        N_all.append(N_bath[ibath])
        res_all.append(resbaths_all)

    # Combining all impurities + all baths
    print("combining baths!")
    print(len(res_all), len(res_all[0]), len(res_all[1]), len(res_all[2]))
    print(N_all)
    # Reverse order
    res_all = res_all[::-1]
    N_all = N_all[::-1]
    res_combined = product_extend_general(res_all, N_all)

    res_combined.sort()
    ndim = len(res_combined)
    f = open(fname, 'w')
    print(ndim, file=f)
    for item in res_combined:
        print(item, file=f)
    f.close()
    return ndim
