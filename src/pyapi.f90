subroutine ed_fsolver(comm, my_id, num_procs, folder)
    use m_control, only: master, origin_myid, origin_nprocs, origin_comm
    use m_control, only: myid, nprocs, new_comm, ndim_i
    use m_global, only: dealloc_fock_i
    use mpi

    implicit none

    integer, intent(in) :: comm
    integer, intent(in) :: my_id
    integer, intent(in) :: num_procs
    character(*), intent(in), optional :: folder

    integer :: ierror
    integer :: color
    integer :: key
    character(len=50) :: folder_

    if (present(folder)) then
        folder_ = folder
    else
        folder_ = "./"
    endif

!F2PY intent(in) comm
!F2PY intent(in) my_id
!F2PY intent(in) num_procs
    origin_comm = comm
    origin_myid = my_id
    origin_nprocs = num_procs

    call config(folder_)  
    ! read fock to know the dimension of the Hamiltonian
    call read_fock_i(folder_)
    call dealloc_fock_i()
    if (ndim_i < origin_nprocs) then
        if (origin_myid==master) then
            print *, " fedrixs >>> Warning: number of CPU processors ", origin_nprocs, &
                     "is larger than ndim_i: ", ndim_i
            print *, " fedrixs >>> Only ", ndim_i, " processors will really work!"
        endif
        if (origin_myid < ndim_i) then
            color = 1 
            key = origin_myid
        else
            color = 2
            key = origin_myid
        endif
        call MPI_COMM_SPLIT(origin_comm, color, key, new_comm, ierror)
        call MPI_COMM_RANK(new_comm, myid, ierror)
        call MPI_COMM_SIZE(new_comm, nprocs, ierror)
    else
        myid = origin_myid
        nprocs = origin_nprocs
        new_comm = origin_comm
    endif

    if (origin_myid < ndim_i) then
        if (present(folder)) then
            call ed_driver(folder)
        else
            call ed_driver()
        endif
    endif

    call MPI_BARRIER(origin_comm, ierror)
    return
end subroutine ed_fsolver

subroutine ed_fsolver_intermediate(comm, my_id, num_procs, folder)
    use m_control, only: master, origin_myid, origin_nprocs, origin_comm
    use m_control, only: myid, nprocs, new_comm, ndim_n_nocore, num_core_orbs
    use m_global, only: dealloc_fock_n
    use mpi

    implicit none

    integer, intent(in) :: comm
    integer, intent(in) :: my_id
    integer, intent(in) :: num_procs
    character(*), intent(in), optional :: folder

    integer :: ndim_n
    integer :: ierror
    integer :: color
    integer :: key
    character(len=50) :: folder_

    if (present(folder)) then
        folder_ = folder
    else
        folder_ = "./"
    endif

    ndim_n = ndim_n_nocore * num_core_orbs

!F2PY intent(in) comm
!F2PY intent(in) my_id
!F2PY intent(in) num_procs
    origin_comm = comm
    origin_myid = my_id
    origin_nprocs = num_procs

    call config(folder_)  
    ! read fock to know the dimension of the Hamiltonian
    call read_fock_n(folder_)
    call dealloc_fock_n()
    if (ndim_n < origin_nprocs) then
        if (origin_myid==master) then
            print *, " fedrixs >>> Warning: number of CPU processors ", origin_nprocs, &
                     "is larger than ndim_n: ", ndim_n
            print *, " fedrixs >>> Only ", ndim_n, " processors will really work!"
        endif
        if (origin_myid < ndim_n) then
            color = 1 
            key = origin_myid
        else
            color = 2
            key = origin_myid
        endif
        call MPI_COMM_SPLIT(origin_comm, color, key, new_comm, ierror)
        call MPI_COMM_RANK(new_comm, myid, ierror)
        call MPI_COMM_SIZE(new_comm, nprocs, ierror)
    else
        myid = origin_myid
        nprocs = origin_nprocs
        new_comm = origin_comm
    endif

    if (origin_myid < ndim_n) then
        if (present(folder)) then
            call ed_driver_intermediate(folder)
        else
            call ed_driver_intermediate()
        endif
    endif

    call MPI_BARRIER(origin_comm, ierror)
    return
end subroutine ed_fsolver_intermediate

subroutine xas_fsolver(comm, my_id, num_procs, folder)
    use m_control, only: master, origin_myid, origin_nprocs, origin_comm
    use m_control, only: myid, nprocs, new_comm, ndim_n, ndim_i, ndim_n_nocore, num_core_orbs
    use m_global, only: dealloc_fock_i, dealloc_fock_n
    use mpi

    implicit none

    integer, intent(in) :: comm
    integer, intent(in) :: my_id
    integer, intent(in) :: num_procs
    character(*), intent(in), optional :: folder

    integer :: ierror
    integer :: color
    integer :: key
    integer :: min_dim
    
    character(len=50) :: folder_

    if (present(folder)) then
        folder_ = folder
    else
        folder_ = "./"
    endif

!F2PY intent(in) comm
!F2PY intent(in) my_id
!F2PY intent(in) num_procs
    origin_comm = comm
    origin_myid = my_id
    origin_nprocs = num_procs

    call config(folder_)  
    call read_fock_i(folder_)
    call dealloc_fock_i()
    call read_fock_n(folder_)
    call dealloc_fock_n()
    ndim_n = ndim_n_nocore * num_core_orbs
    min_dim = min(ndim_i, ndim_n)
    if (min_dim < origin_nprocs) then
        if (origin_myid==master) then
            print *, " fedrixs >>> Warning: number of CPU processors ", origin_nprocs, &
                     "is larger than min(ndim_i, ndim_n): ", ndim_i, ndim_n
            print *, " fedrixs >>> Only ", min_dim, " processors will really work!"
        endif
        if (origin_myid < min_dim) then
            color = 1 
            key = origin_myid
        else
            color = 2
            key = origin_myid
        endif
        call MPI_COMM_SPLIT(origin_comm, color, key, new_comm, ierror)
        call MPI_COMM_RANK(new_comm, myid, ierror)
        call MPI_COMM_SIZE(new_comm, nprocs, ierror)
    else
        myid = origin_myid
        nprocs = origin_nprocs
        new_comm = origin_comm
    endif

    if (origin_myid < min_dim) then
        if (present(folder)) then
            call xas_driver(folder)
        else
            call xas_driver()
        endif
    endif

    call MPI_BARRIER(origin_comm, ierror)

    return
end subroutine xas_fsolver

subroutine rixs_fsolver(comm, my_id, num_procs, folder)
    use m_control, only: master, origin_myid, origin_nprocs, origin_comm
    use m_control, only: myid, nprocs, new_comm, ndim_n, ndim_i, ndim_f, ndim_n_nocore, num_core_orbs
    use m_global, only: dealloc_fock_i, dealloc_fock_n, dealloc_fock_f
    use mpi

    implicit none

    integer, intent(in) :: comm
    integer, intent(in) :: my_id
    integer, intent(in) :: num_procs
    character(*), intent(in), optional :: folder

    integer :: ierror
    integer :: color
    integer :: key
    integer :: min_dim
    
    character(len=50) :: folder_

    if (present(folder)) then
        folder_ = folder
    else
        folder_ = "./"
    endif

!F2PY intent(in) comm
!F2PY intent(in) my_id
!F2PY intent(in) num_procs

    origin_comm = comm
    origin_myid = my_id
    origin_nprocs = num_procs

    call config(folder_)  
    call read_fock_i(folder_)
    call dealloc_fock_i()
    call read_fock_n(folder_)
    call dealloc_fock_n()
    call read_fock_f(folder_)
    call dealloc_fock_f()
    ndim_n = ndim_n_nocore * num_core_orbs
    min_dim = min(ndim_i, ndim_n, ndim_f)
    if (min_dim < origin_nprocs) then
        if (origin_myid==master) then
            print *, " fedrixs >>> Warning: number of CPU processors ", origin_nprocs, &
                     "is larger than min(ndim_i, ndim_n, ndim_f): ", ndim_i, ndim_n, ndim_f
            print *, " fedrixs >>> Only ", min_dim, " processors will really work!"
        endif
        if (origin_myid < min_dim) then
            color = 1 
            key = origin_myid
        else
            color = 2
            key = origin_myid
        endif
        call MPI_COMM_SPLIT(origin_comm, color, key, new_comm, ierror)
        call MPI_COMM_RANK(new_comm, myid, ierror)
        call MPI_COMM_SIZE(new_comm, nprocs, ierror)
    else
        myid = origin_myid
        nprocs = origin_nprocs
        new_comm = origin_comm
    endif

    if (origin_myid < min_dim) then
        if (present(folder)) then
            call rixs_driver(folder)
        else
            call rixs_driver()
        endif
    endif

    call MPI_BARRIER(origin_comm, ierror)
    return
end subroutine rixs_fsolver

subroutine opavg_fsolver(comm, my_id, num_procs)
    use m_control, only: master, origin_myid, origin_nprocs, origin_comm
    use m_control, only: myid, nprocs, new_comm, ndim_i
    use m_global, only: dealloc_fock_i
    use mpi

    implicit none

    integer, intent(in) :: comm
    integer, intent(in) :: my_id
    integer, intent(in) :: num_procs

    integer :: ierror
    integer :: color
    integer :: key

!F2PY intent(in) comm
!F2PY intent(in) my_id
!F2PY intent(in) num_procs
    origin_comm = comm
    origin_myid = my_id
    origin_nprocs = num_procs

    call config()  
    call read_fock_i()
    call dealloc_fock_i()
    if (ndim_i < origin_nprocs) then
        if (origin_myid==master) then
            print *, " fedrixs >>> Warning: number of CPU processors ", origin_nprocs, &
                     "is larger than ndim_i: ", ndim_i
            print *, " fedrixs >>> Only ", ndim_i, " processors will really work!"
        endif
        if (origin_myid < ndim_i) then
            color = 1 
            key = origin_myid
        else
            color = 2
            key = origin_myid
        endif
        call MPI_COMM_SPLIT(origin_comm, color, key, new_comm, ierror)
        call MPI_COMM_RANK(new_comm, myid, ierror)
        call MPI_COMM_SIZE(new_comm, nprocs, ierror)
    else
        myid = origin_myid
        nprocs = origin_nprocs
        new_comm = origin_comm
    endif
    print *, " fedrixs >>> ", origin_myid, origin_nprocs, myid, nprocs

    if (origin_myid < ndim_i) then
        call opavg_driver()
    endif

    call MPI_BARRIER(origin_comm, ierror)

    return
end subroutine opavg_fsolver
