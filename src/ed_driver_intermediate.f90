subroutine ed_driver_intermediate(folder)
    use m_constants
    use m_control 
    use m_types
    use m_global
    use m_lanczos
    use mpi

    ! ...........................................
    !       ED for the intermediate states.
    ! ...........................................
    
    implicit none

    character(*), intent(in), optional :: folder

    integer :: i,j,k
    integer :: icfg, jcfg
    integer :: mblock
    integer :: nblock
    integer :: nloc
    integer :: needed(nprocs,nprocs)
    integer :: end_indx(2,2,nprocs)
    integer :: ierror 
    integer :: info
    integer :: actual_step
    integer :: nconv
    integer :: request
    integer :: stat(MPI_STATUS_SIZE)
    integer :: sgn
    integer :: tot_sign
    integer(dp) :: old, new
    integer(dp) :: num_of_nonzeros
    integer, external :: binary_search

    real(dp) :: rtemp
    real(dp) :: time_begin
    real(dp) :: time_end
    real(dp) :: time_used
    real(dp), allocatable :: eigvals(:)
    real(dp), allocatable :: eigval_full(:)

    complex(dp) :: omega
    complex(dp) :: denmat(num_val_orbs, num_val_orbs, nvector)
    complex(dp) :: denmat_mpi(num_val_orbs, num_val_orbs, nvector)

    complex(dp), allocatable :: ham_full(:,:)
    complex(dp), allocatable :: ham_full_mpi(:,:)
    complex(dp), allocatable :: eigvecs_full(:,:)

    complex(dp), allocatable :: eigvecs(:,:)
    complex(dp), allocatable :: eigvecs_mpi(:)

    character(len=100) :: folder_
    character(len=100) :: fname
    !character(len=20) :: fname
    character(len=10) :: char_I

    if (present(folder)) then
        folder_ = folder
    else
        folder_ = "./"
    endif

    time_begin = 0.0_dp
    time_end   = 0.0_dp
    time_used  = 0.0_dp
    call cpu_time(time_begin)

    if (myid == master) then
        print *, "--------------------------------------------"
        print *, " fedrixs >>> ED begin (intermediate) ... "
        print *
    endif
    call read_hopping_n(folder_)
    call read_coulomb_n(folder_)
    call read_fock_n(folder_)
    
    ! ................................................
    ! ndim_n = ndim_n_nocore * num_core_orbs 
    ! This is set in subroutine xas_fsolver, but
    !       now I just want to do ED.
    !   dimension is counted with just 1 core-hole
    ! ................................................
    ndim_n = ndim_n_nocore * num_core_orbs

    ! min_ndim : parameter related to the parallelization
    if (ndim_n < min_ndim ) then
        if (myid==master) then
            print *, " fedrixs >>> ndim_n:", ndim_n, " is smaller than min_ndim:", min_ndim
            print *, " fedrixs >>> set ed_solver = 0, use full-diagonalization !"
            print *
        endif
        ed_solver = 0
    endif

    if(ndim_n < neval) then
        neval = ndim_n
    endif
    if(ndim_n < nvector) then
        nvector=ndim_n
    endif

    if (myid == master) then
        write(mystd,"(a20, i15)")   "num_val_orbs:  ", num_val_orbs
        write(mystd,"(a20, i15)")   "ed_solver:     ", ed_solver
        write(mystd,"(a20, i15)")   "neval:         ", neval
        write(mystd,"(a20, i15)")   "nvector:       ", nvector
        write(mystd,"(a20, i15)")   "maxiter:       ", maxiter
        write(mystd,"(a20, i15)")   "min_ndim:      ", min_ndim
        write(mystd,"(a20, i15)")   "ncv:           ", ncv
        write(mystd,"(a20, i15)")   "nhopp_n:       ", nhopp_n
        write(mystd,"(a20, i15)")   "ncoul_n:       ", ncoul_n
        write(mystd,"(a20, i15)")   "ndim_n_nocore: ", ndim_n_nocore
        write(mystd,"(a20, i15)")   "ndim_n:        ", ndim_n
        write(mystd,"(a20, i15)")   "nprocs:        ", nprocs

        write(mystd,"(a20, e15.2)") "eigval_tol:    ", eigval_tol
        write(mystd,"(a20, l15)")   "idump:         ", idump
        print *
    endif

    if (myid==master) then
        print *, " fedrixs >>> Build Hamiltonian ..."
    endif
    call partition_task(nprocs, ndim_n, ndim_n, end_indx)
    rtemp = 1.0_dp
    omega = czero
    mblock = nprocs
    nblock = nprocs
    nloc = end_indx(2,2,myid+1)-end_indx(1,2,myid+1) + 1
    call alloc_ham_csr(nblock)
    call build_ham_n(ndim_n_nocore, fock_n, num_val_orbs, num_core_orbs, nblock, end_indx, &
        nhopp_n, hopping_n, ncoul_n, coulomb_n, omega, rtemp, ham_csr)
    call MPI_BARRIER(new_comm, ierror)
    call get_needed_indx(nblock, ham_csr, needed)
    call get_number_nonzeros(nblock, ham_csr, num_of_nonzeros)
    call dealloc_fock_i()
    call cpu_time(time_end)
    time_used = time_used + time_end - time_begin
    if (myid==master) then
        print *, " fedrixs >>> Number of nonzero elements of the Hamiltonian", num_of_nonzeros
        print *, " fedrixs >>> Done ! Time used: ", time_end - time_begin, "  seconds"
        print *
        print *, " fedrixs >>> Diagonalize Hamiltonian to find a few lowest states ..."
        print *
    endif
    time_begin = time_end

    allocate(eigvals(neval))
    ! full diagonalization
    if (ed_solver == 0) then
        allocate(ham_full(ndim_n, ndim_n))
        allocate(ham_full_mpi(ndim_n, ndim_n))
        allocate(eigvecs_full(ndim_n, ndim_n))
        allocate(eigval_full(ndim_n))
        allocate(eigvecs(ndim_n,nvector))
        ham_full = czero
        ham_full_mpi = czero
        eigvecs_full = czero
        eigval_full = zero
        do i=1,nblock
            call csr_to_full(ndim_n, ndim_n, ham_csr(i), ham_full_mpi)
        enddo
        call MPI_ALLREDUCE(ham_full_mpi, ham_full, size(ham_full_mpi), MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
        if (myid==master) then
            call full_diag_ham(ndim_n, ham_full, eigval_full, eigvecs_full) 
        endif
        call MPI_BARRIER(new_comm, ierror)
        call MPI_BCAST(eigval_full,   size(eigval_full),  MPI_DOUBLE_PRECISION, master, new_comm, ierror)
        call MPI_BCAST(eigvecs_full,  size(eigvecs_full), MPI_DOUBLE_COMPLEX, master, new_comm, ierror)
        eigvals = eigval_full(1:neval)
        eigvecs = eigvecs_full(:,1:nvector)
        deallocate(ham_full)
        deallocate(ham_full_mpi)
        deallocate(eigvecs_full)
        deallocate(eigval_full)
    ! ordinary Lanczos method
    elseif (ed_solver == 1) then
        eigvals = zero
        allocate(eigvecs(nloc,nvector))
        eigvecs = czero
        call diag_ham_lanczos(nblock, end_indx, needed, nloc, neval, maxiter, eigval_tol, ham_csr, eigvals, nvector, eigvecs)
    ! arpack is used
    elseif (ed_solver==2) then
        eigvals = zero
        allocate(eigvecs(nloc,nvector))
        eigvecs = czero
        call diag_ham_arpack(nblock, end_indx, needed, nloc, neval, ham_csr, eigvals, nvector, &
                             eigvecs, ncv, maxiter, eigval_tol, info, actual_step, nconv)
        if (myid==master) then
            print *, "    arpack return info:  ", info
            print *, "    actual arnoldi step:  ", actual_step
            print *, "    number of converged Ritz values:  ", nconv
            print *
        endif
    endif
    call MPI_BARRIER(new_comm, ierror)
    call dealloc_ham_csr(nblock)
    call cpu_time(time_end)
    time_used = time_used + time_end - time_begin
    if (myid==master) then
        print *
        print *, " fedrixs >>> Done ! Time used: ", time_end-time_begin, "  seconds"
        print *
    endif
    time_begin=time_end

    !fname="eigvec."//trim(adjustl(char_I))
    fname=trim(folder_)//"eigvals_n.dat"
    call write_lowest_eigvals(fname, neval, eigvals)

    if (myid==master) then
        print *, " fedrixs >>> Calculate the density matrix ... "
    endif

    !call read_fock_n(folder_)
    denmat = czero
    denmat_mpi = czero
    allocate(eigvecs_mpi(ndim_n))    
    do k=1, nvector
        eigvecs_mpi = czero
        if (ed_solver==0) then
            eigvecs_mpi = eigvecs(:,k)
        else
            ! send the data
            do i=1,nprocs
               if (myid+1 /= i) then
                   call MPI_ISEND(eigvecs(:,k), nloc, MPI_DOUBLE_COMPLEX, i-1, &
                          i*(10*nprocs)+myid+1, new_comm, request, ierror)
               endif
            enddo
            eigvecs_mpi(end_indx(1,2,myid+1):end_indx(2,2,myid+1)) = eigvecs(:,k)

            ! receive data
            do i=1,nprocs
               if (myid+1 /= i) then
                   call MPI_RECV(eigvecs_mpi(end_indx(1,2,i):end_indx(2,2,i)), &
                        end_indx(2,2,i)-end_indx(1,2,i)+1, MPI_DOUBLE_COMPLEX, &
                       i-1, (myid+1)*(10*nprocs)+i, new_comm, stat, ierror)
               endif
            enddo
        endif
        call MPI_BARRIER(new_comm, ierror)

        !do icfg=end_indx(1,1,myid+1), end_indx(2,1,myid+1)
        !    do j=1,num_val_orbs
        !        if (btest(fock_n(icfg), j-1)) then
        !            denmat_mpi(j,j,k) = denmat_mpi(j,j,k) + conjg(eigvecs_mpi(icfg)) * eigvecs_mpi(icfg) 
        !        endif
        !    enddo

        !    if ( abs(eigvecs_mpi(icfg)) < 1E-10 ) cycle
        !    do i=1,num_val_orbs-1
        !        do j=i+1,num_val_orbs
        !            if (.not. btest(fock_n(icfg), i-1)) cycle
        !            old = fock_n(icfg)

        !            tot_sign = 1
        !            call make_newfock('-', i, old, new, sgn) 
        !            tot_sign = tot_sign * sgn
        !            old = new

        !            if (btest(old, j-1)) cycle
        !            call make_newfock('+', j,  old, new, sgn) 
        !            tot_sign = tot_sign * sgn

        !            jcfg = binary_search(ndim_n, fock_n, new)
                    
         !           if (jcfg == -1) cycle
         !           denmat_mpi(i,j,k) = denmat_mpi(i,j,k) + conjg(eigvecs_mpi(icfg)) * eigvecs_mpi(jcfg) * tot_sign
         !           denmat_mpi(j,i,k) = denmat_mpi(j,i,k) + conjg(eigvecs_mpi(jcfg)) * eigvecs_mpi(icfg) * tot_sign
         !       enddo
         !   enddo
        !enddo

        ! dump eigenvectors 
        if ( idump ) then 
            write(char_I, '(i5)') k
            fname=trim(folder_)//"eigvec_n."//trim(adjustl(char_I))
            if (myid==master) then
                call write_eigvecs(fname, ndim_n, eigvecs_mpi, eigvals(k))
            endif
        endif
    enddo  ! over k=1, nvector

    call MPI_BARRIER(new_comm, ierror)
    !call MPI_ALLREDUCE(denmat_mpi, denmat, size(denmat_mpi), MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
    
    !fname=trim(folder_)//"denmat_n.dat"
    !call write_denmat(fname, nvector, num_val_orbs, denmat)

    call dealloc_fock_n()
    call dealloc_hopping_n()
    call dealloc_coulomb_n()
    deallocate(eigvecs)
    !deallocate(eigvecs_mpi)
    deallocate(eigvals)

    call cpu_time(time_end)
    time_used = time_used + time_end - time_begin
    if (myid==master) then
        print *, " fedrixs >>> Done ! Time used: ", time_end-time_begin, "  seconds"
        print *
        print *, " fedrixs >>> ED end ! Total time used: ", time_used, "  seconds"
        print *
    endif

    return
end subroutine ed_driver_intermediate
