  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_mix(next_dkn_in, dkn_in, dkn_out, residues, overlap, &
       iter, ientry)

    !==========================================================================!
    ! Selects method for kernel DIIS.                                          !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in November 2010.                         !
    !==========================================================================!

    use comms, only: pub_on_root, comms_abort
    use constants, only: stdout
    use rundat, only: pub_kernel_diis_type, pub_kernel_diis_max, &
         pub_kernel_diis_liter
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_copy

    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: dkn_in(:,:)
    type(SPAM3), intent(inout) :: dkn_out(:,:)
    type(SPAM3), intent(inout) :: residues(:,:)
    type(SPAM3), intent(in   ) :: overlap
    integer,     intent(in   ) :: iter
    integer,     intent(in   ) :: ientry


    integer :: is

    select case(pub_kernel_diis_type)
    case('D')
       ! ars: do diagonalisation only
       do is=1, pub_cell%num_spins
          call sparse_copy(next_dkn_in(is), dkn_out(is, ientry))
       end do

    case('L')
       ! ars: do linear mixing only
       call kernel_diis_linear_mixing(next_dkn_in, dkn_in, &
            dkn_out, ientry)

    case('P')
       !ars: do Pulay mixing. Optional linear mixing for the first iterations.
       if(iter.le.pub_kernel_diis_liter) then
          call kernel_diis_linear_mixing(next_dkn_in, dkn_in, &
               dkn_out, ientry)
       else
          call kernel_diis_pulay_mixing(next_dkn_in, &
               dkn_out, residues, overlap, ientry)
       end if

    case default
       if(pub_on_root) write(stdout,*) "Unknown mixing method. ONETEP stops."
       call comms_abort

    end select

  end subroutine kernel_diis_mix


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~ Linear mixing routines ~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine kernel_diis_linear_mixing(next_dkn_in, dkn_in, dkn_out, ientry)

    !==========================================================================!
    ! Performs linear mixing of density kernels.                               !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in June 2010.                             !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, VERBOSE
    use rundat, only: pub_kernel_diis_max, pub_kernel_diis_c_in, &
         pub_kernel_diis_c_out, pub_output_detail
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3

    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: dkn_in(:,:)
    type(SPAM3), intent(inout) :: dkn_out(:,:)
    integer,     intent(in   ) :: ientry



    if(pub_on_root.and.pub_output_detail.ge.VERBOSE) then
       write(stdout,'(/,a)') "... PERFORMING LINEAR MIXING ..."
       write(stdout,'(a9,f6.4,a10,f6.4/)') &
            " C_in = ", pub_kernel_diis_c_in, ", C_out = ", pub_kernel_diis_c_out
    end if


    ! ars: mix kernels linearly
    call kernel_diis_linear_mix_kernels(next_dkn_in, dkn_out, &
         dkn_in, ientry)


  end subroutine kernel_diis_linear_mixing



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_linear_mix_kernels(next_dkn_in, dkn_out, dkn_in, ientry)

    !==================================================================!
    ! This subroutine performs a linear mixing of the density kernels. !
    !------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in June 2010.                      !
    !==================================================================!

    use constants, only: DP
    use rundat, only: pub_kernel_diis_c_in, pub_kernel_diis_c_out, &
         pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_axpy, sparse_scale

    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: dkn_out(:,:)
    type(SPAM3), intent(inout) :: dkn_in(:,:)
    integer,     intent(in   ) :: ientry

    integer :: is


    do is = 1, pub_cell%num_spins
       ! ars: K_in(n+1) = 0
       call sparse_scale(next_dkn_in(is), 0.0_DP)
       ! ars: K_in(n+1) = c_out*K_out(n)
       call sparse_axpy(next_dkn_in(is), dkn_out(is,ientry), &
            pub_kernel_diis_c_out)
       ! ars: K_in(n+1) = c_out*K_out(n) + c_in*K_in(n)
       call sparse_axpy(next_dkn_in(is), dkn_in(is,ientry), &
            pub_kernel_diis_c_in)
    end do

  end subroutine kernel_diis_linear_mix_kernels



  !_____________________________________________________________________________
  !_____________________________________________________________________________


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~ Pulay mixing routines ~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine kernel_diis_pulay_mixing(next_dkn_in,dkn_out,residues,overlap,ientry)

    !==========================================================================!
    ! Performs Pulay mixing of density kernels.                                !
    ! Based on Pulay, Chem. Phys. Lett. 73, 2, 1980.                           !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in June 2010.                             !
    ! Updated by Alvaro Ruiz Serrano in October 2010 to combine linear+Pulay   !
    ! mixing schemes.                                                          !
    !==========================================================================!


    use comms, only: pub_on_root
    use constants, only: DP, stdout, VERBOSE
    use rundat, only: pub_kernel_diis_max, pub_output_detail
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check


    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: dkn_out(:,:)
    type(SPAM3), intent(inout) :: residues(:,:)
    type(SPAM3), intent(in)    :: overlap
    integer,     intent(in)    :: ientry


    ! ars: local variables
    real(kind=DP), allocatable :: Bmat(:,:,:)!ars: linear system to be solved
    real(kind=DP), allocatable :: coeffs(:,:)!ars: Pulay mixing coeffs
    integer :: ierr


    if(pub_on_root.and.pub_output_detail.ge.VERBOSE) then
       write(stdout,'(/,a)') "... PERFORMING PULAY MIXING ..."
       write(stdout,'(a)') "Coefficients and Lagrange multiplier:"
    end if

    ! ars: allocate Bmat and coeffs
    allocate(Bmat(1:pub_cell%num_spins,1:ientry+1,1:ientry+1),stat=ierr)
    call utils_alloc_check('kernel_diis_mix_density', 'Bmat', ierr)
    allocate(coeffs(1:pub_cell%num_spins,1:ientry+1), stat=ierr)
    call utils_alloc_check('kernel_diis_mix_density', 'coeffs', ierr)

    ! ars: calculate B matrix.
    write(*,*) 'pulay calc Bmatrix, ientry : ', ientry
    call kernel_diis_pulay_calc_Bmatrix(Bmat, residues, overlap, ientry)

    ! ars: solve linear system Bd=0 (including Lagrange mult)
    write(*,*) 'pulay find coeffs, ientry : ', ientry
    call kernel_diis_pulay_find_coeffs(coeffs, Bmat, ientry)
    if(pub_output_detail.ge.VERBOSE) &
         call kernel_diis_pulay_coeffs_summ(coeffs, ientry)

    ! ars: mix kernel
    write(*,*) 'pulay mix kerneles, ientry : ', ientry
    call kernel_diis_pulay_mix_kernels(next_dkn_in, dkn_out, coeffs, ientry)
 
    write(*,*) 'pulay deallocate'

    ! ars: deallocate Bmat and coeffs
    deallocate(coeffs, stat=ierr)
    call utils_dealloc_check('kernel_diis_mix_density', 'coeffs', ierr)
    deallocate(Bmat, stat=ierr)
    call utils_dealloc_check('kernel_diis_mix_density', 'Bmat', ierr)


  end subroutine kernel_diis_pulay_mixing

  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_pulay_calc_Bmatrix(Bmat,residues, overlap, ientry)

    !==============================================================!
    ! This subroutine calculates matrix B of the system of linear  !
    ! equations to be solved in order to mix the density kernels   !
    ! according to the Pulay method.                               !
    !--------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in April 2010.                 !
    !==============================================================!

    use constants, only: stdout, DP
    use comms, only: pub_on_root, pub_my_node_id, comms_barrier
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_transpose, sparse_product,&
         sparse_create, sparse_destroy, sparse_trace
    implicit none
    real(kind=DP), intent(inout) :: Bmat(:,:,:)
    type(SPAM3),   intent(inout) :: residues(:,:)
    type(SPAM3),   intent(in)    :: overlap
    integer,       intent(in)    :: ientry
    real(kind=DP)                :: Bmat_tmp(ientry,ientry)
    integer                      :: table(ientry*ientry*num_spins,3),kkk,kk,iini,istep
    integer                      :: is, ii, jj
    type(SPAM3)                  :: spam3buffer1, spam3buffer2, spam3buffer3

    call sparse_create(spam3buffer1, residues(1,1), overlap)
    call sparse_create(spam3buffer2, residues(1,1))
    call sparse_create(spam3buffer3, spam3buffer1, spam3buffer2)

    kkk=0
    do jj = 1, ientry
     do ii = 1, ientry
      do is = 1, num_spins
       kkk=kkk+1
       table(kkk,1)=jj
       table(kkk,2)=ii
       table(kkk,3)=is
      enddo
     enddo
    enddo

    iini=0; istep=1; Bmat=0.0

    if(.not. use_gpu_some_only)then
      iini=mpi_rank-1; istep=mpi_size
    else
      if(mpi_rank/=1.and.wont_compute_residues) then
       write(*,*) 'ERROR : WRITING KERNEL WITH 1 GPU : MPI_RANK = ' , mpi_rank
       write(*,*) 'WAS NOT SUPPOSED TO BE HERE'
       stop
      endif
      if(mpi_rank/=1) goto 8888
    endif 

    do kk=iini+1,kkk,istep

             if(pub_my_node_id==0) write(*,'(a,i4,2i6)') 'MPI_RANK SCANNING : ', mpi_rank,kk,kkk

             jj = table(kk,1); ii = table(kk,2); is = table(kk,3);

             if(pub_my_node_id==0.and.verboseall) write(*,*) 'ii , jj , is : ',ii,jj,is

#ifndef debugMIXv3
             if(use_gpu_onlyme) then
               call matmul_in_dense_format_r(spam3buffer1, residues(is,ii), overlap)
             else
#endif
               ! ars: spam3buffer1 = R_i x S
               call sparse_product(spam3buffer1,residues(is,ii),overlap)
#ifndef debugMIXv3
             endif
#endif
             ! ars: spam3buffer2 = (R_j)^t
             call sparse_transpose(spam3buffer2,residues(is,jj))

#ifndef debugMIXv3
             if(use_gpu_onlyme)then
               call matmul_in_dense_format_r(spam3buffer3,spam3buffer1,spam3buffer2)
             else
#endif
             ! ars: spam3buffer3 = (R_i) x S x (R_j)^t
               call sparse_product(spam3buffer3, spam3buffer1, spam3buffer2)
#ifndef debugMIXv3
             endif
#endif
             ! ars: Bmat(is,ii,jj) = tr[(R_i) x S x (R_j)^t x S]
             Bmat(is,ii,jj) = sparse_trace(spam3buffer3,overlap)

             if(pub_my_node_id==0.and.verboseall) write(*,'(a,2f20.10)') 'B() = ', Bmat(is,ii,jj)

    end do
 
    if(.not.wont_compute_residues)then
     if(mpi_size/=1)then
        8888 continue
          is=1
          Bmat_tmp(1:ientry,1:ientry)=Bmat(is,1:ientry,1:ientry)
          write(*,*) 'MPI_RANK CALLING NFS SUM_B : ', mpi_rank,mpi_size
          write(*,*) 'IENTRY SIZE                : ', ientry
          call sync_nodes_by_shared_nfs_mat(Bmat_tmp,'sumBup',mpi_rank,mpi_size,ientry,longdelay=.true.,safe=.true.)
          write(*,*) 'SYNC DONE, RANK            : ', mpi_rank,mpi_size
          Bmat(is,1:ientry,1:ientry)=Bmat_tmp(1:ientry,1:ientry)
        if(num_spins==2)then
          is=2
          Bmat_tmp(1:ientry,1:ientry)=Bmat(is,1:ientry,1:ientry)
          write(*,*) 'MPI_RANK CALLING NFS SUM_B : ', mpi_rank,mpi_size
          write(*,*) 'IENTRY SIZE                : ', ientry
          call sync_nodes_by_shared_nfs_mat(Bmat_tmp,'sumBdown',mpi_rank,mpi_size,ientry,longdelay=.true.,safe=.true.)
          write(*,*) 'SYNC DONE, RANK            : ', mpi_rank,mpi_size
          Bmat(is,1:ientry,1:ientry)=Bmat_tmp(1:ientry,1:ientry)
        endif
     endif
    endif

    if(num_spins/=pub_cell%num_spins)then
       if(pub_cell%num_spins==1)then
         write(*,*) 'ERROR ,   num_spins = ', num_spins
         write(*,*) 'pub_cell%num_spins  = ', pub_cell%num_spins
         stop
       endif
       Bmat(2,:,:)=Bmat(1,:,:)
     endif

    ! ars: destroy SPAM3 buffers
    call sparse_destroy(spam3buffer1)
    call sparse_destroy(spam3buffer2)
    call sparse_destroy(spam3buffer3)

    ! ars: add space for Langrange multiplier
    Bmat(:,1:ientry,ientry+1) = -1.0_DP
    Bmat(:,ientry+1,1:ientry) = -1.0_DP
    Bmat(:,ientry+1,ientry+1) = 0.0_DP

  end subroutine

  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_pulay_find_coeffs(coeffs,Bmat,ientry)

    !==============================================================!
    ! This subroutine finds the coefficients for the Pulay mixing  !
    ! after solving a system of linear equations Bd=0.             !
    !--------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in May 2010.                   !
    !==============================================================!

    use comms, only: pub_my_node_id, pub_on_root, comms_barrier
    use constants, only: DP, stdout
    use simulation_cell, only: pub_cell


    implicit none

    ! ars: arguments
    real(kind=DP), intent(inout) :: coeffs(:,:)
    real(kind=DP), intent(in   ) :: Bmat(:,:,:)
    integer      , intent(in   ) :: ientry

    ! ars: library wrapper variables
    integer :: INFO, LDA, LDB, M, N, NRHS
    integer :: IPIV(1:ientry+1)
    real(kind=DP) :: A(ientry+1,ientry+1)
    real(kind=DP) :: B(ientry+1)
    character(LEN=1) :: TRANS

    ! ars: local variables
    integer :: is



    ! ars: set up parameters for LAPACK DGETRF
    M = ientry+1
    N = ientry+1
    LDA = ientry+1
    ! ars: set up parameters for LAPACK DGETRS
    TRANS ='N'
    NRHS = 1
    LDB = ientry+1

    ! ars: solve linear system
    do is = 1, pub_cell%num_spins

       ! ars: call LAPACK DGETRF
       A = Bmat(is,:,:)
!#ifdef debug
       call write_array( A , ' matrix Bcoef '//trim(adjustl(tostring(ientry+1))) )
!#endif
       call DGETRF(M, N, A, LDA, IPIV, INFO)
       if((INFO.ne.0).and.pub_on_root) write(stdout,*) &
            "Error calling DGETRF. INFO|node = ", INFO, pub_my_node_id

       ! ars: call LAPACK DGETRS
       B(1:ientry) = 0.0_DP
       B(ientry+1) = -1.0_DP
       call DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       if((INFO.ne.0).and.pub_on_root) write(stdout,*) &
            "Error calling DGETRS. INFO|node = ", INFO, pub_my_node_id

       ! ars: set coeffs(is) before exit
       coeffs(is,:) = B(:)

    end do

    call comms_barrier

    if(pub_dmft_impose_same_coeffs.and.pub_cell%num_spins==2) coeffs(2,:)=coeffs(1,:)

  end subroutine kernel_diis_pulay_find_coeffs

  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_pulay_coeffs_summ(coeffs, ientry)

    !==================================================================!
    ! This subroutine prints a summary of the coefficients that are    !
    ! used during the density kernel mixing.                           !
    !------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in July 2010.                      !
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: DP, stdout
    use simulation_cell, only: pub_cell

    implicit none

    real(kind=DP), intent(in   ) :: coeffs(:,:)
    integer,       intent(in   ) :: ientry

    integer :: is, it

    if(pub_on_root) then
       do it = 1, ientry+1
          do is = 1, pub_cell%num_spins

             if(it.le.ientry) then
                write(stdout,'(a6,i1,a1,i2,a4,f20.12)') &
                     "Coeff(",is,",",it,") = ", coeffs(is,it)
             else
                write(stdout,'(a7,i1,a6,f20.12,/)') &
                     "Lambda(",is,")   = ", coeffs(is,it)
             end if

          end do
       end do
    end if

  end subroutine kernel_diis_pulay_coeffs_summ



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_pulay_mix_kernels(next_dkn_in, dkn_out, coeffs, ientry)

    !==================================================================!
    ! This subroutine mixes the kernels according to the Pulay method. !
    !------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in May 2010.                       !
    !==================================================================!

    use constants, only: DP
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_axpy, sparse_scale

    implicit none

    type(SPAM3),   intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: dkn_out(:,:)
    real(kind=DP), intent(in   ) :: coeffs(:,:)
    integer,       intent(in   ) :: ientry

    integer :: is, counter

    ! ars: K_in(n+1) = 0
    do is = 1, pub_cell%num_spins
       call sparse_scale(next_dkn_in(is),0.0_DP)
    end do

    ! ars: K_in(n+1) = sum_i K_out(i) * d_i
    do counter = 1, ientry
       do is = 1, pub_cell%num_spins
          call sparse_axpy(next_dkn_in(is),dkn_out(is,counter),coeffs(is,counter))
       end do
    end do


  end subroutine

  !_____________________________________________________________________________
  !_____________________________________________________________________________

  !~~~~~~~~~~~~~~~~~~
  !~~ Common routines
  !~~~~~~~~~~~~~~~~~~
  subroutine kernel_diis_sparse_create(dkn_in, dkn_out, next_dkn_in, residues, denskern)

    !===================================================================!
    ! This subroutine creates the SPAM3 structures for dkn_in,          !
    ! dkn_out, next_dkn_in, residues and hamiltonian_array              !
    !-------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010.                     !
    !===================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_create, SPAM3

    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: dkn_in(:,:)
    type(SPAM3), intent(inout) :: dkn_out(:,:)
    type(SPAM3), intent(inout) :: residues(:,:)
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3)                :: denskern(pub_cell%num_spins)

    ! ars: local variables
    integer :: is, iter


    ! ars: create SPAM3 arrays
    do iter = 1, pub_kernel_diis_max
       do is = 1, pub_cell%num_spins
!#ifndef debugMIXv1
!          dkn_in(is,iter)%structure = 'K'
!          call sparse_create(dkn_in(is,iter))
!          dkn_out(is,iter)%structure = 'K'
!          call sparse_create(dkn_out(is,iter))
!          residues(is,iter)%structure = 'K'
!          call sparse_create(residues(is,iter))
!#else
          dkn_in(is,iter)%structure = 'K'
          call sparse_create(dkn_in(is,iter),denskern(1))
          dkn_out(is,iter)%structure = 'K'
          call sparse_create(dkn_out(is,iter),denskern(1))
          residues(is,iter)%structure = 'K'
          call sparse_create(residues(is,iter),denskern(1))
!#endif
       end do
    end do

    do is=1, pub_cell%num_spins
       next_dkn_in%structure = 'K'
       call sparse_create(next_dkn_in(is))
    end do


  end subroutine



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_sparse_destroy(dkn_in, dkn_out, next_dkn_in, residues)

    !======================================================================!
    ! This subroutine destroys the SPAM3 structures for for dkn_in, !
    ! dkn_out, next_dkn_in, residues and hamiltonian_array.   !
    !----------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in April 2010.                         !
    !======================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_destroy, SPAM3


    implicit none


    ! ars: arguments
    type(SPAM3), intent(inout) :: dkn_in(:,:)
    type(SPAM3), intent(inout) :: dkn_out(:,:)
    type(SPAM3), intent(inout) :: residues(:,:)
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)

    ! ars: local variables
    integer :: is, iter


    ! ars: destroy SPAM3 arrays
    do iter = 1, pub_kernel_diis_max
       do is = 1, pub_cell%num_spins
          call sparse_destroy(dkn_in(is,iter))
          call sparse_destroy(dkn_out(is,iter))
          call sparse_destroy(residues(is,iter))
       end do
    end do

    do is=1, pub_cell%num_spins
       call sparse_destroy(next_dkn_in(is))
    end do


  end subroutine kernel_diis_sparse_destroy




  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_init(residues, next_dkn_in, dkn_out, dkn_in, denskern)

    !====================================================!
    ! This subroutine initialises the first element of   !
    ! dkn_in by copying the latest density kernel !
    ! onto it (coming from the last NGWF iteration).     !
    !----------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in April 2010.       !
    !====================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, DP
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_copy, sparse_scale

    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: residues(:,:)
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: dkn_out(:,:)
    type(SPAM3), intent(inout) :: dkn_in(:,:)
    type(SPAM3), intent(in   ) :: denskern(pub_cell%num_spins)

    ! ars: local variables
    integer :: is, it


    ! ars: initialise all SPAM3 matrices to zero
    do it = 1, pub_kernel_diis_max
       do is = 1, pub_cell%num_spins
          call sparse_scale(residues(is,it), 0.0_DP)
          call sparse_scale(dkn_in(is,it), 0.0_DP)
          call sparse_scale(dkn_out(is,it), 0.0_DP)
       end do
    end do

    do is = 1, pub_cell%num_spins
       call sparse_scale(next_dkn_in(is), 0.0_DP)
    end do

    ! ars: initialise first ientry of dkn_in
    do is = 1, pub_cell%num_spins
       call sparse_copy(dkn_in(is,1),denskern(is))
    end do

  end subroutine kernel_diis_init



  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_residue_inout(residue,kern_out,kern_in)

    !====================================================================!
    ! This subroutine calculates the density kernel residue after each   !
    ! kernel DIIS iteration. R_i = K^out_i - K^in_i.                     !
    !--------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010.                      !
    !====================================================================!

    use comms, only: pub_on_root
    use constants, only: DP
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_copy, sparse_axpy, SPAM3


    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: residue(pub_cell%num_spins)
    type(SPAM3), intent(in   ) :: kern_out(pub_cell%num_spins)
    type(SPAM3), intent(in   ) :: kern_in(pub_cell%num_spins)
    integer                    :: is, nspins

    do is = 1,pub_cell%num_spins
       call sparse_copy(residue(is),kern_out(is))
       ! calculate [R = K^out - K^in]
       call sparse_axpy(residue(is),kern_in(is),-1.0_DP)
    end do

  end subroutine

  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_shift(shifted, dkn_in, dkn_out, residues, iter)

    !================================================================!
    ! This subroutine shifts the matrices in arrays in case iter has !
    ! reached pub_kernel_diis_max value.                             !
    !----------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in May 2010.                     !
    !================================================================!

    use constants, only: DP, stdout
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_copy

    implicit none

    ! ars: arguments
    logical,     intent(  out) :: shifted
    type(SPAM3), intent(inout) :: dkn_in(:,:)
    type(SPAM3), intent(inout) :: dkn_out(:,:)
    type(SPAM3), intent(inout) :: residues(:,:)
    integer,     intent(in   ) :: iter

    ! ars: local variables
    integer :: scntr, is

    shifted = .false.

    ! ars: move old entries one position back
    if (iter.ge.pub_kernel_diis_max) then
       shifted = .true.
       do scntr = 1, pub_kernel_diis_max-1
          do is = 1, pub_cell%num_spins
             call sparse_copy(dkn_in(is,scntr), dkn_in(is, scntr+1))
             call sparse_copy(dkn_out(is,scntr), dkn_out(is, scntr+1))
             call sparse_copy(residues(is,scntr), residues(is, scntr+1))
          end do
       end do
    end if

  end subroutine kernel_diis_shift



  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_find_ientry(ientry, iter, shifted)

    !==============================================================!
    ! This subroutine finds the next ientry to work with.           !
    !--------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in May 2010.                   !
    !==============================================================!

    use comms, only: comms_barrier
    use rundat, only: pub_kernel_diis_max

    implicit none

    ! ars: arguments
    integer, intent(  out) :: ientry
    integer, intent(in   ) :: iter
    logical, intent(in   ) :: shifted


    ! ars: find next position
    if (shifted) then
       ientry = pub_kernel_diis_max
    else
       ientry = iter+1
    end if

    call comms_barrier

  end subroutine kernel_diis_find_ientry

  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !_____________________________________________________________________________
