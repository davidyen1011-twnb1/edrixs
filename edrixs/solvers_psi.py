__all__ = ['ed_siam_fort_general']

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
    v_norb_multi: int array, 
        Number of total orbitals of impurity+baths(except for the last one), only implemented for do_ed=1
    multibaths : Truee / False. Whether activate multibaths calculations.
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

    constrained_basis = False

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
    emat_n[0:ntot_v, 0:ntot_v] -= np.eye(v_norb) * static_core_pot

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
            emat_i[last_v:last_v+v_norb[imp], last_v:last_v+v_norb[imp]] += np.eye(v_norb) * eval_shift
            emat_n[ntot_v+last_c:ntot_v+last_c+c_norb[imp], ntot_v+last_c:ntot_v+last_c+c_norb[imp]] += np.eye(c_norb) * c_level
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

