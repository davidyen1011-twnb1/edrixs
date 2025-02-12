__all__ = ['ed_siam_fort_general', 'xas_siam_fort_general', 'rixs_siam_fort_general']

import numpy as np
import scipy

from .iostream import (
    write_tensor, write_emat, write_umat, write_config, read_poles_from_file
)
from .angular_momentum import (
    get_sx, get_sy, get_sz, get_lx, get_ly, get_lz, rmat_to_euler, get_wigner_dmat
)
from .photon_transition import (
    get_trans_oper, quadrupole_polvec, dipole_polvec_xas, dipole_polvec_rixs, unit_wavevector
)
from .coulomb_utensor import get_umat_slater, get_umat_slater_3shells
from .manybody_operator import two_fermion, four_fermion
from .fock_basis import get_fock_bin_by_N, write_fock_dec_by_N,\
        write_fock_dec_by_N_constrainedN1N2,\
        write_fock_dec_by_N_constrainedN1N2_multi,\
        write_fock_dec_by_N_general
from .basis_transform import cb_op2, tmat_r2c, cb_op
from .utils import info_atomic_shell, slater_integrals_name, boltz_dist
from .rixs_utils import scattering_mat
from .plot_spectrum import get_spectra_from_poles, merge_pole_dicts
from .soc import atom_hsoc


def ed_siam_fort_general(comm, c_name, *, static_core_pot=0, c_level=0,
                 c_soc=0, ext_B=None, hopping=None, hopping_n=None, on_which='spin', do_ed=1, ed_solver=2, neval=1,
                 nvector=1, ncv=3, idump=False, maxiter=1000, eigval_tol=1e-8, min_ndim=1000,
                 umat_input_i=None, umat_input_n=None, folder="./", v_norb=None, c_norb=None,\
                 b_norb=None, v_orbl=None, v_noccu_imp=None, v_noccu_baths=None):
    """
    Find the ground state of the initial Hamiltonian of a Single Impuirty Anderson Model (SIAM),
    and also prepare input files, *hopping_i.in*, *hopping_n.in*, *coulomb_i.in*, *coulomb_n.in*
    for following XAS and RIXS calculations.

    For multi-impurities and multi-baths for convenient constrained basis computation.

    Parameters
    ----------
    comm: MPI_Comm
        MPI Communicator
    static_core_pot: float
        Static core hole potential.
    c_name : string
        Name of the core shell, options are 's', 'p', 'd', 'f'.
    c_level: float
        Energy level of core shell.
    c_soc: float
        Spin-orbit coupling strength of core electrons.
    ext_B: tuple of three float numbers
        Vector of external magnetic field with respect to global :math:`xyz`-axis.

        They will be set to zero if not provided.
    on_which: string
        Apply Zeeman exchange field on which sector. Options are 'spin', 'orbital' or 'both'.
    hopping: 2d complex array
        General hopping matrix when siam_type=1, including imp_mat and hybridization functions,
        for siam_type=1 and the initial configurations.
    hopping_n: 2d complex array
        General hopping matrix when siam_type=1, including imp_mat and hybridization functions,
        for siam_type=1 and the intermediate configurations. If hopping_n=None,
        hopping will be used.
    do_ed: int
        - 1: Only do ED for given occupancy number *v_noccu*, return eigenvalues and
          density matrix, write eigenvectors to files eigvec.n
    ed_solver: int
        Type of ED solver, options can be 0, 1, 2

        - 0: use Lapack to fully diagonalize Hamiltonian to get all the eigenvalues.

        - 1: use standard Lanczos algorithm to find only a few lowest eigenvalues,
          no re-orthogonalization has been applied, so it is not very accurate.

        - 2: use parallel version of Arpack library to find a few lowest eigenvalues,
          it is accurate and is the recommeded choice in real calculations of XAS and RIXS.
    neval: int
        Number of eigenvalues to be found. For ed_solver=2, the value should not be too small,
        neval > 10 is usually a safe value.
    nvector: int
        Number of eigenvectors to be found and written into files.
    ncv: int
        Used for ed_solver=2, it should be at least ncv > neval + 2. Usually, set it a little
        bit larger than neval, for example, set ncv=200 when neval=100.
    idump: logical
        Whether to dump the eigenvectors to files "eigvec.n", where n means the n-th vectors.
    maxiter: int
        Maximum number of iterations in finding all the eigenvalues, used for ed_solver=1, 2.
    eigval_tol: float
        The convergence criteria of eigenvalues, used for ed_solver=1, 2.
    min_ndim: int
        The minimum dimension of the Hamiltonian when the ed_solver=1, 2 can be used, otherwise,
        ed_solver=1 will be used.
    umat_input_i : rank 4 tensor, interaction tensor Umat for initial states
    umat_input_n : rank 4 tensor, interaction tensor Umat for intermediate states
    v_orbl [n_imp] : l quantum number for each impurity sites
    v_norb [n_imp] : impurity valence orbitals
    c_norb [n_imp] : impurity core orbitals
    b_norb [n_baths] : bath orbitals
    v_noccu_imp [n_imp, *possible occupations] : occupancy of impurity sites
    v_noccu_baths [n_baths, *possible occupations] : occupancy of bath sites
    

    Returns
    -------
    eval_i: 1d float array
        Eigenvalues of initial Hamiltonian.
    denmat: 2d complex array
        Density matrix.
    noccu_gs: int
        Occupancy of the ground state.
    """
    try:
        from .fedrixs import ed_fsolver
    except:
        from fedrixs import ed_fsolver

    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()
    if rank == 0:
        print("edrixs >>> Running ED ...", flush=True)

    ntot_v = np.sum(v_norb) + np.sum(b_norb)    
    ntot_c = np.sum(c_norb)
    ntot = ntot_v + np.sum(c_norb)                  # v_imp + c_imp + v_bath
    ntot_imp = np.sum(v_norb) + np.sum(c_norb)      # v_imp + c_imp
    n_imp = len(v_norb)                             # number of impurity sites

    umat_i = np.zeros((ntot, ntot, ntot, ntot), dtype=complex)
    umat_n = np.zeros((ntot, ntot, ntot, ntot), dtype=complex)
    
    # Index for impurity valence + core (Bath is not correlated here)
    # [Imp1, Imp2,...bath1, bath2..., core1, core2...]
    last_v = 0
    last_c = 0
    indx = []
    for i in range(n_imp):
        indx.append(list(range(last_v, last_v+v_norb[i])) + [ntot_v + i for i in range(last_c, last_c+c_norb[i])])
        last_v += v_norb[i]
        last_c += c_norb[i]

    # Initial states
    for imp in range(n_imp):
        for i in range(v_norb[imp]+c_norb[imp]):
            for j in range(v_norb[imp]+c_norb[imp]):
                for k in range(v_norb[imp]+c_norb[imp]):
                    for m in range(v_norb[imp]+c_norb[imp]):
                        umat_i[indx[imp][i], indx[imp][j], indx[imp][k], indx[imp][m]] = umat_input_i[i, j, k, m]
    
    # Intermediate states
    for imp in range(n_imp):
        for i in range(v_norb[imp]+c_norb[imp]):
            for j in range(v_norb[imp]+c_norb[imp]):
                for k in range(v_norb[imp]+c_norb[imp]):
                    for m in range(v_norb[imp]+c_norb[imp]):
                        umat_n[indx[imp][i], indx[imp][j], indx[imp][k], indx[imp][m]] = umat_input_n[i, j, k, m]

    if rank == 0:
        write_umat(umat_i, folder+'coulomb_i.in')
        write_umat(umat_n, folder+'coulomb_n.in')

    emat_i = np.zeros((ntot, ntot), dtype=complex)
    emat_n = np.zeros((ntot, ntot), dtype=complex)

    emat_i[0:ntot_v, 0:ntot_v] += hopping
    emat_n[0:ntot_v, 0:ntot_v] += hopping_n

    # atomic SOC
    if c_name in ['p', 'd', 'f']:
        last_c = 0
        for imp in range(n_imp):
            emat_n[ntot_v+last_c:ntot_v+last_c+c_norb[imp], ntot_v+last_c:ntot_v+last_c+c_norb[imp]] +=\
                atom_hsoc(c_name, c_soc)
            last_c += c_norb[imp]

    # Static core potential
    last_v = 0
    for imp in range(n_imp):
        emat_n[last_v:last_v+v_norb[imp], last_v:last_v+v_norb[imp]] -= np.eye(v_norb[imp]) * static_core_pot
        last_v += v_norb[imp]

    # Zeeman terms
    last_v = 0
    for imp in range(n_imp):
        # L & S operator for each impurity
        lx, ly, lz = get_lx(v_orbl[imp], True), get_ly(v_orbl[imp], True), get_lz(v_orbl[imp], True)
        sx, sy, sz = get_sx(v_orbl[imp]), get_sy(v_orbl[imp]), get_sz(v_orbl[imp])
        if ext_B is not None:
            # For single site case
            if on_which.strip() == 'spin':
                zeeman = ext_B[0] * (2 * sx) + ext_B[1] * (2 * sy) + ext_B[2] * (2 * sz)
            elif on_which.strip() == 'orbital':
                # There should be a scale factor somewhere, init?
                zeeman = ext_B[0] * lx + ext_B[1] * ly + ext_B[2] * lz
            elif on_which.strip() == 'both':
                # A factor somewhere?
                zeeman = ext_B[0] * (lx + 2 * sx) + ext_B[1] * (ly + 2 * sy) + ext_B[2] * (lz + 2 * sz)
            else:
                raise Exception("Unknown value of on_which", on_which)
            
            emat_i[last_v:last_v+v_norb[imp], last_v:last_v+v_norb[imp]] += zeeman
            emat_n[last_v:last_v+v_norb[imp], last_v:last_v+v_norb[imp]] += zeeman
        
        last_v += v_norb[imp]

    # Perform ED 
    if do_ed == 1:

        v_noccu = v_noccu_imp[:,0] + np.sum(v_noccu_baths[:,0])     # [n_imp]
        # Shift the core level
        eval_shift = c_level * c_norb / v_noccu                     # [n_imp]
        last_v = 0
        last_c = 0
        for imp in range(n_imp):
            emat_i[last_v:last_v+v_norb[imp], last_v:last_v+v_norb[imp]] += np.eye(v_norb[imp]) * eval_shift[imp]
            emat_n[ntot_v+last_c:ntot_v+last_c+c_norb[imp], ntot_v+last_c:ntot_v+last_c+c_norb[imp]] += np.eye(c_norb[imp]) * c_level[imp]
            last_v += v_norb[imp]
            last_c += c_norb[imp]

        if rank == 0:
            # Write hopping files
            write_emat(emat_i, folder+'hopping_i.in')
            write_emat(emat_n, folder+'hopping_n.in')
            write_config(
                directory=folder, ed_solver=ed_solver, num_val_orbs=ntot_v, neval=neval, nvector=nvector, ncv=ncv,
                idump=idump, maxiter=maxiter, min_ndim=min_ndim, eigval_tol=eigval_tol
            )
            
            # Using Constrained basis
            print("..................................... ")
            print("Using constrained basis. ")
            print("..................................... ")
            write_fock_dec_by_N_general(v_norb, v_noccu_imp,\
                        b_norb, v_noccu_baths, folder+"fock_i.in")
        
            print("edrixs >>> do_ed=1, perform ED at noccu: ", v_noccu, flush=True)
        
        # Run ED from here
        comm.Barrier()
        ed_fsolver(fcomm, rank, size, folder)
        comm.Barrier()
        data = np.loadtxt(folder+'eigvals.dat', ndmin=2)
        eval_i = np.zeros(neval, dtype=float)
        eval_i[0:neval] = data[0:neval, 1]
        data = np.loadtxt(folder+'denmat.dat', ndmin=2)
        tmp = (nvector, ntot_v, ntot_v)
        denmat = data[:, 3].reshape(tmp) + 1j * data[:, 4].reshape(tmp)
        return eval_i, denmat, v_noccu

    else:
        return None, None, None
        raise Exception("Unknown case of do_ed ", do_ed)

def xas_siam_fort_general(comm, ominc, v_name, c_name,*, gamma_c=0.1, thin=1.0, phi=0, pol_type=None,
                num_gs=1, nkryl=200, temperature=1.0,
                loc_axis=None, scatter_axis=None, folder="./",\
                v_norb=None, c_norb=None, b_norb=None, v_orbl=None,\
                v_noccu_imp=None, v_noccu_baths=None, v_noccu_imp_n=None, v_noccu_baths_n=None):
    """
    Calculate XAS for single impurity Anderson model (SIAM) with Fortran solver.

    Parameters
    ----------
    comm: MPI_comm
        MPI communicator.
    v_name : string for valence shell ['s', 'p', 't2g', 'd', 'f']
    c_name : string for core shell ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']
    ominc: 1d float array
        Incident energy of photon.
    gamma_c: a float number or a 1d float array with the same shape as ominc.
        The core-hole life-time broadening factor. It can be a constant value
        or incident energy dependent.
    v_noccu: int
        Total occupancy of valence shells.
    thin: float number
        The incident angle of photon (in radian).
    phi: float number
        Azimuthal angle (in radian), defined with respect to the
        :math:`x`-axis of the local scattering axis: scatter_axis[:,0].
    pol_type: list of tuples
        Type of polarization, options can be:

        - ('linear', alpha), linear polarization, where alpha is the angle between the
          polarization vector and the scattering plane in radians.

        - ('left', 0), left circular polarization.

        - ('right', 0), right circular polarization.

        - ('isotropic', 0). isotropic polarization.

        It will set pol_type=[('isotropic', 0)] if not provided.
    num_gs: int
        Number of initial states used in XAS calculations.
    nkryl: int
        Maximum number of poles obtained.
    temperature: float number
        Temperature (in K) for boltzmann distribution.
    loc_axis: 3*3 float array
        The local axis with respect to which local orbitals are defined.

        - x: local_axis[:,0],

        - y: local_axis[:,1],

        - z: local_axis[:,2].

        It will be an identity matrix if not provided.
    scatter_axis: 3*3 float array
        The local axis defining the scattering geometry. The scattering plane is defined in
        the local :math:`zx`-plane.

        - local :math:`x`-axis: scatter_axis[:,0]

        - local :math:`y`-axis: scatter_axis[:,1]

        - local :math:`z`-axis: scatter_axis[:,2]

        It will be set to an identity matrix if not provided.
    v_orbl [n_imp] : l quantum number for each impurity sites
    v_norb [n_imp] : impurity valence orbitals
    c_norb [n_imp] : impurity core orbitals
    b_norb [n_baths] : bath orbitals
    v_noccu_imp: int array, 
        Number of total occupancy of impurity
    v_noccu_baths: int array, 
        Number of total occupancy of baths
    v_noccu_imp_n: int array, 
        Number of total occupancy of impurity (For intermediate states)
    v_noccu_baths_n: int array, 
        Number of total occupancy of baths (For intermediate states)

    Returns
    -------
    xas: 2d array, shape=(len(ominc), len(pol_type))
        The calculated XAS spectra. The first dimension is for ominc, and the second dimension
        if for different polarizations.
    poles: list of dict, shape=(len(pol_type), )
        The calculated XAS poles for different polarizations.
    """
    try:
        from .fedrixs import xas_fsolver
    except:
        from fedrixs import xas_fsolver

    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()

    ntot_v = np.sum(v_norb) + np.sum(b_norb)    
    ntot_c = np.sum(c_norb)
    ntot_b = np.sum(b_norb)
    ntot = ntot_v + np.sum(c_norb)                  # v_imp + c_imp + v_bath
    ntot_imp = np.sum(v_norb) + np.sum(c_norb)      # v_imp + c_imp
    n_imp = len(v_norb)                             # number of impurity sites

    # Polarization
    if pol_type is None:
        pol_type = [('isotropic', 0)]
    if loc_axis is None:
        loc_axis = np.eye(3)
    else:
        loc_axis = np.array(loc_axis)
    if scatter_axis is None:
        scatter_axis = np.eye(3)
    else:
        scatter_axis = np.array(scatter_axis)

    if rank == 0:
        print("edrixs >>> Running XAS ...", flush=True)
        write_config(directory=folder,num_val_orbs=ntot_v, num_core_orbs=ntot_c,
                     num_gs=num_gs, nkryl=nkryl)

        # Using Constrained basis
        print("..................................... ")
        print("Using constrained basis. ")
        print("..................................... ")
        write_fock_dec_by_N_general(v_norb, v_noccu_imp,\
                    b_norb, v_noccu_baths, folder+"fock_i.in")
        write_fock_dec_by_N_general(v_norb, v_noccu_imp_n,\
                    b_norb, v_noccu_baths_n, folder+"fock_n.in")

    case = v_name + c_name
    tmp = get_trans_oper(case)
    npol, n, m = tmp.shape
    tmp_g = np.zeros((npol, n, m), dtype=complex)
    trans_mat = np.zeros((npol, ntot, ntot), dtype=complex)
    # Transform the transition operators to global-xyz axis
    # dipolar transition
    if npol == 3:
        for i in range(3):
            for j in range(3):
                tmp_g[i] += loc_axis[i, j] * tmp[j]
    # quadrupolar transition
    elif npol == 5:
        alpha, beta, gamma = rmat_to_euler(loc_axis)
        wignerD = get_wigner_dmat(4, alpha, beta, gamma)
        rotmat = np.dot(np.dot(tmat_r2c('d'), wignerD), np.conj(np.transpose(tmat_r2c('d'))))
        for i in range(5):
            for j in range(5):
                tmp_g[i] += rotmat[i, j] * tmp[j]
    else:
        raise Exception("Have NOT implemented this case: ", npol)
    #trans_mat[:, 0:v_norb, ntot_v:ntot] = tmp_g

    last_v = 0
    last_c = 0
    for i in range(n_imp):
        trans_mat[:, last_v:last_v+v_norb[i], ntot_v+last_c:ntot_v+last_c+c_norb[i]] = tmp_g
        last_v += v_norb[i]
        last_c += c_norb[i]

    n_om = len(ominc)
    gamma_core = np.zeros(n_om, dtype=float)
    if np.isscalar(gamma_c):
        gamma_core[:] = np.ones(n_om) * gamma_c
    else:
        gamma_core[:] = gamma_c

    # loop over different polarization
    xas = np.zeros((n_om, len(pol_type)), dtype=float)
    poles = []
    comm.Barrier()
    for it, (pt, alpha) in enumerate(pol_type):
        if pt.strip() == 'left' or pt.strip() == 'right' or pt.strip() == 'linear':
            if rank == 0:
                print("edrixs >>> Loop over for polarization: ", it, pt, flush=True)
                kvec = unit_wavevector(thin, phi, scatter_axis, 'in')
                polvec = np.zeros(npol, dtype=complex)
                pol = dipole_polvec_xas(thin, phi, alpha, scatter_axis, pt)
                if npol == 3:  # Dipolar transition
                    polvec[:] = pol
                if npol == 5:  # Quadrupolar transition
                    polvec[:] = quadrupole_polvec(pol, kvec)

                trans = np.zeros((ntot, ntot), dtype=complex)
                for i in range(npol):
                    trans[:, :] += trans_mat[i] * polvec[i]
                write_emat(trans, folder+'transop_xas.in')

            # call XAS solver in fedrixs
            comm.Barrier()
            xas_fsolver(fcomm, rank, size, folder)
            comm.Barrier()

            # Read the spectrum from 'xas_poles.' files...
            file_list = [folder+'xas_poles.' + str(i+1) for i in range(num_gs)]
            pole_dict = read_poles_from_file(file_list)
            poles.append(pole_dict)
            xas[:, it] = get_spectra_from_poles(pole_dict, ominc, gamma_core, temperature)
        elif pt.strip() == 'isotropic':
            pole_dicts = []
            for k in range(npol):
                if rank == 0:
                    print("edrixs >>> Loop over for polarization: ", it, pt, flush=True)
                    print("edrixs >>> Isotropic, component: ", k, flush=True)
                    write_emat(trans_mat[k], folder+'transop_xas.in')

                # call XAS solver in fedrixs
                comm.Barrier()
                xas_fsolver(fcomm, rank, size, folder)
                comm.Barrier()

                # Read the spectrum from 'xas_poles.' files...
                file_list = [folder+'xas_poles.' + str(i+1) for i in range(num_gs)]
                pole_tmp = read_poles_from_file(file_list)
                xas[:, it] += get_spectra_from_poles(pole_tmp, ominc, gamma_core, temperature)
                pole_dicts.append(pole_tmp)
            xas[:, it] = xas[:, it] / npol
            poles.append(merge_pole_dicts(pole_dicts))
        else:
            raise Exception("Unknown polarization type: ", pt)

    return xas, poles


def rixs_siam_fort_general(comm, ominc, eloss, v_name, c_name,*, gamma_c=0.1, gamma_f=0.1,
                   thin=1.0, thout=1.0, phi=0, pol_type=None, num_gs=1,
                   nkryl=200, linsys_max=1000, linsys_tol=1e-10, temperature=1.0,
                   loc_axis=None, scatter_axis=None, folder="./",\
                   v_norb=None, c_norb=None, b_norb=None, v_orbl=None,\
                    v_noccu_imp=None, v_noccu_baths=None, v_noccu_imp_n=None, v_noccu_baths_n=None):
    """
    Calculate RIXS for single impurity Anderson model with Fortran solver.

    Parameters
    ----------
    comm: MPI_comm
        MPI communicator.
    v_name : string for valence shell ['s', 'p', 't2g', 'd', 'f']
    c_name : string for core shell ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']
    ominc: 1d float array
        Incident energy of photon.
    eloss: 1d float array
        Energy loss.
    gamma_c: a float number or a 1d float array with same shape as ominc.
        The core-hole life-time broadening factor. It can be a constant value
        or incident energy dependent.
    gamma_f: a float number or a 1d float array with same shape as eloss.
        The final states life-time broadening factor. It can be a constant value
        or energy loss dependent.
    thin: float number
        The incident angle of photon (in radian).
    thout: float number
        The scattered angle of photon (in radian).
    phi: float number
        Azimuthal angle (in radian), defined with respect to the
        :math:`x`-axis of scattering axis: scatter_axis[:,0].
    pol_type: list of 4-elements-tuples
        Type of polarizations. It has the following form:

        (str1, alpha, str2, beta)

        where, str1 (str2) can be 'linear', 'left', 'right', and alpha (beta) is
        the angle (in radians) between the linear polarization vector and the scattering plane.

        It will set pol_type=[('linear', 0, 'linear', 0)] if not provided.
    num_gs: int
        Number of initial states used in RIXS calculations.
    nkryl: int
        Maximum number of poles obtained.
    linsys_max: int
        Maximum iterations of solving linear equations.
    linsys_tol: float
        Convergence for solving linear equations.
    temperature: float number
        Temperature (in K) for boltzmann distribution.
    loc_axis: 3*3 float array
        The local axis with respect to which local orbitals are defined.

        - x: local_axis[:,0],

        - y: local_axis[:,1],

        - z: local_axis[:,2].

        It will be an identity matrix if not provided.
    scatter_axis: 3*3 float array
        The local axis defining the scattering geometry. The scattering plane is defined in
        the local :math:`zx`-plane.

        - local :math:`x`-axis: scatter_axis[:,0]

        - local :math:`y`-axis: scatter_axis[:,1]

        - local :math:`z`-axis: scatter_axis[:,2]

        It will be set to an identity matrix if not provided.
    v_orbl [n_imp] : l quantum number for each impurity sites
    v_norb [n_imp] : impurity valence orbitals
    c_norb [n_imp] : impurity core orbitals
    b_norb [n_baths] : bath orbitals
    v_noccu_imp: int array, 
        Number of total occupancy of impurity
    v_noccu_baths: int array, 
        Number of total occupancy of baths
    v_noccu_imp_n: int array, 
        Number of total occupancy of impurity (For intermediate states)
    v_noccu_baths_n: int array, 
        Number of total occupancy of baths (For intermediate states)

    Returns
    -------
    rixs: 3d float array, shape=(len(ominc), len(eloss), len(pol_type))
        The calculated RIXS spectra. The 1st dimension is for the incident energy,
        the 2nd dimension is for the energy loss and the 3rd dimension is for
        different polarizations.
    poles: 2d list of dict, shape=(len(ominc), len(pol_type))
        The calculated RIXS poles. The 1st dimension is for incident energy, and the
        2nd dimension is for different polarizations.
    """
    try:
        from .fedrixs import rixs_fsolver
    except:
        from fedrixs import rixs_fsolver

    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()

    # v_name, c_name
    ntot_v = np.sum(v_norb) + np.sum(b_norb)    
    ntot_c = np.sum(c_norb)
    ntot_b = np.sum(b_norb)
    ntot = ntot_v + np.sum(c_norb)                  # v_imp + c_imp + v_bath
    ntot_imp = np.sum(v_norb) + np.sum(c_norb)      # v_imp + c_imp
    n_imp = len(v_norb)      

    if pol_type is None:
        pol_type = [('linear', 0, 'linear', 0)]
    if loc_axis is None:
        loc_axis = np.eye(3)
    else:
        loc_axis = np.array(loc_axis)
    if scatter_axis is None:
        scatter_axis = np.eye(3)
    else:
        scatter_axis = np.array(scatter_axis)

    if rank == 0:
        print("edrixs >>> Running RIXS ...", flush=True)
        # Using Constrained basis
        print("..................................... ")
        print("Using constrained basis. ")
        print("..................................... ")
        write_fock_dec_by_N_general(v_norb, v_noccu_imp,\
                    b_norb, v_noccu_baths, folder+"fock_i.in")
        write_fock_dec_by_N_general(v_norb, v_noccu_imp_n,\
                    b_norb, v_noccu_baths_n, folder+"fock_n.in")
        write_fock_dec_by_N_general(v_norb, v_noccu_imp,\
                    b_norb, v_noccu_baths, folder+"fock_f.in")

        case = v_name + c_name
        tmp = get_trans_oper(case)
        npol, n, m = tmp.shape
        tmp_g = np.zeros((npol, n, m), dtype=complex)
        trans_mat = np.zeros((npol, ntot, ntot), dtype=complex)
        # Transform the transition operators to global-xyz axis
        # dipolar transition
        if npol == 3:
            for i in range(3):
                for j in range(3):
                    tmp_g[i] += loc_axis[i, j] * tmp[j]
        # quadrupolar transition
        elif npol == 5:
            alpha, beta, gamma = rmat_to_euler(loc_axis)
            wignerD = get_wigner_dmat(4, alpha, beta, gamma)
            rotmat = np.dot(np.dot(tmat_r2c('d'), wignerD), np.conj(np.transpose(tmat_r2c('d'))))
            for i in range(5):
                for j in range(5):
                    tmp_g[i] += rotmat[i, j] * tmp[j]
        else:
            raise Exception("Have NOT implemented this case: ", npol)
        #trans_mat[:, 0:v_norb, ntot_v:ntot] = tmp_g

        last_v = 0
        last_c = 0
        for i in range(n_imp):
            trans_mat[:, last_v:last_v+v_norb[i], ntot_v+last_c:ntot_v+last_c+c_norb[i]] = tmp_g
            last_v += v_norb[i]
            last_c += c_norb[i]

    n_om = len(ominc)
    neloss = len(eloss)
    gamma_core = np.zeros(n_om, dtype=float)
    if np.isscalar(gamma_c):
        gamma_core[:] = np.ones(n_om) * gamma_c
    else:
        gamma_core[:] = gamma_c
    gamma_final = np.zeros(neloss, dtype=float)
    if np.isscalar(gamma_f):
        gamma_final[:] = np.ones(neloss) * gamma_f
    else:
        gamma_final[:] = gamma_f

    # loop over different polarization
    rixs = np.zeros((n_om, neloss, len(pol_type)), dtype=float)
    poles = []
    comm.Barrier()
    # loop over different polarization
    for iom, omega in enumerate(ominc):
        if rank == 0:
            write_config(
                directory=folder, num_val_orbs=ntot_v, num_core_orbs=c_norb,
                omega_in=omega, gamma_in=gamma_core[iom],
                num_gs=num_gs, nkryl=nkryl, linsys_max=linsys_max,
                linsys_tol=linsys_tol
            )
        poles_per_om = []
        # loop over polarization
        for ip, (it, alpha, jt, beta) in enumerate(pol_type):
            if rank == 0:
                print(flush=True)
                print("edrixs >>> Calculate RIXS for incident energy: ", omega, flush=True)
                print("edrixs >>> Polarization: ", ip, flush=True)
                polvec_i = np.zeros(npol, dtype=complex)
                polvec_f = np.zeros(npol, dtype=complex)
                ei, ef = dipole_polvec_rixs(thin, thout, phi, alpha, beta,
                                            scatter_axis, (it, jt))
                # dipolar transition
                if npol == 3:
                    polvec_i[:] = ei
                    polvec_f[:] = ef
                # quadrupolar transition
                elif npol == 5:
                    ki = unit_wavevector(thin, phi, scatter_axis, direction='in')
                    kf = unit_wavevector(thout, phi, scatter_axis, direction='out')
                    polvec_i[:] = quadrupole_polvec(ei, ki)
                    polvec_f[:] = quadrupole_polvec(ef, kf)
                else:
                    raise Exception("Have NOT implemented this type of transition operators")
                trans_i = np.zeros((ntot, ntot), dtype=complex)
                trans_f = np.zeros((ntot, ntot), dtype=complex)
                for i in range(npol):
                    trans_i[:, :] += trans_mat[i] * polvec_i[i]
                write_emat(trans_i, folder+'transop_rixs_i.in')
                for i in range(npol):
                    trans_f[:, :] += trans_mat[i] * polvec_f[i]
                write_emat(np.conj(np.transpose(trans_f)), folder+'transop_rixs_f.in')

            # call RIXS solver in fedrixs
            comm.Barrier()
            rixs_fsolver(fcomm, rank, size, folder)
            comm.Barrier()

            file_list = [folder+'rixs_poles.' + str(i+1) for i in range(num_gs)]
            pole_dict = read_poles_from_file(file_list)
            poles_per_om.append(pole_dict)
            rixs[iom, :, ip] = get_spectra_from_poles(pole_dict, eloss,
                                                      gamma_final, temperature)

        poles.append(poles_per_om)

    return rixs, poles

