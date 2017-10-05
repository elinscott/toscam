! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                  Electronic transport module                   !
!                                                                !
!----------------------------------------------------------------!
! Written by Simon M.-M. Dubois in April 2010                    !
!================================================================!


module transport

  use constants, only: dp
  use dense, only: DEM

  implicit none

  private

  public :: etrans_calculate

  ! Type definition 
  type, public :: block_info
     character(len=1)   :: type    ! L/R = left/right electrode, D = full device
     real(kind=dp)      :: efermi  ! Characteristic Fermi Energy
     integer            :: atms(4) ! Range of atoms (original order)
     integer            :: orbs(4) ! Range of orbitals (original order)
  end type block_info

  type, public :: block_hs
     integer                     :: norb  
     type(DEM)                   :: s00, s01
     type(DEM),allocatable       :: h00(:), h01(:)
  end type block_hs

  type, public :: epath
     character(len=2)   :: type      ! TR = transmission/dos
     integer            :: nep       ! total number of points
     integer            :: seg(4)    ! (nreal, ncircle, nhline, nvline)

     integer, allocatable          :: nepnode(:)
     integer, allocatable          :: nepptr(:)

     complex(kind=dp), allocatable :: ep(:)
     complex(kind=dp), allocatable :: ew(:)
  end type epath

  ! Public variables
  integer, public, save :: pub_nleads
  type(block_info), public, save :: info_device
  type(block_info), allocatable, public, save  :: info_leads(:)

  integer, save :: idrain, isource
  type(block_hs), save :: hs_device
  type(block_hs), allocatable, save  :: hs_leads(:)
 
contains

!====================================================================!
!====================================================================!

  subroutine etrans_calculate(ham,rep,ngwf_basis,elements)

    !================================================================!
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in April 2010                    !
    !================================================================!

    use constants, only: dp, pi, two_pi, stdout, stderr, cmplx_0, &
         cmplx_1, cmplx_i, hartree_in_evs
    use comms, only: pub_my_node_id, pub_on_root, comms_barrier, &
         pub_total_num_nodes, pub_root_node_id, comms_send, comms_recv, &
         comms_abort
    use dense, only: DEM, dense_create, dense_destroy
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use rundat, only: pub_rootname, pub_etrans_same_leads, &
         pub_etrans_source, pub_etrans_drain, task
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit
  
    implicit none

    ! Arguments
    type(NGWF_HAM), intent(in)   :: ham
    type(NGWF_REP), intent(in)   :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(ELEMENT), intent(in)    :: elements(pub_cell%nat)

    ! Local variables
    type(epath)    :: etrc
    complex(kind=DP)  :: energy

    type(DEM), allocatable       :: greenf(:)
    type(DEM), allocatable       :: self(:)
    real(kind=DP), allocatable   :: dos(:,:), trc(:,:) 

    integer  :: nspin, norb
    integer  :: ieloc, ienergy, is, iorb, jorb, inode
    integer  :: ierr, info, output_unit
    character(len=80)   :: filename
   
    call timer_clock('etrans_calculate',1)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering etrans_calculate'
#endif

    ! If required, initialise device setup
    if (task .ne. 'NEGF') then
       if (pub_on_root) write(stdout,'(a)') 'etrans : Initialise transport setup !  '
       call transport_init_etrans_setup(ham,rep,ngwf_basis,elements)
    endif 

    nspin = pub_cell%num_spins
    norb  = hs_device%norb

    ! Initialise energy contour for transmission coefficient
    if (pub_on_root) write(stdout,'(a)') 'etrans : Initialise energy contour ! '
    call transport_init_energy_contour(etrc,'trans',info_device%efermi)

    allocate(trc(etrc%nep,nspin),stat=ierr)
    call utils_alloc_check('etrans_calculate','trc',ierr)
    allocate(dos(etrc%nep,nspin),stat=ierr)
    call utils_alloc_check('etrans_calculate','dos',ierr)
    trc = 0.0_dp
    dos = 0.0_dp

    ! Allocate self-energies and Green's functions
    allocate(self(nspin),stat=ierr)
    call utils_alloc_check('etrans_calculate','self',ierr)
    allocate(greenf(nspin),stat=ierr)
    call utils_alloc_check('etrans_calculate','greenf',ierr)

    do is = 1, nspin
       call dense_create(self(is),norb,norb,.true.)
       call dense_create(greenf(is),norb,norb,.true.)
    enddo

    if (pub_on_root) write(stdout,'(a)') &
         'etrans : Compute DOS and transmission coefficients ! '

    ! Main loop over the energy
    energy_loop: do ieloc = 1, etrc%nepnode(pub_my_node_id+1)
                                    
       ienergy = etrc%nepptr(pub_my_node_id+1) + ieloc  
       energy = etrc%ep(ienergy) 

       write(filename,'(a,a,i4.4)') trim(adjustl(pub_rootname)),'.e',ienergy
       call transport_greenf_calculate(energy,greenf,self,filename)

       do is = 1, nspin

          ! Compute the DOS (i.e. greenf*overlap)
          do iorb = 1 , norb
             do jorb = 1, norb
                dos(ienergy,is) = dos(ienergy,is) - &
                   dimag(greenf(is)%zmtx(iorb,jorb)*hs_device%s00%dmtx(jorb,iorb))
            enddo
          enddo
          dos(ienergy,is) = dos(ienergy,is) / PI
    
          ! Calculate the transmission matrix 
          call transport_trc_calculate(self(is),greenf(is),trc(ienergy,is))

       enddo
    
    enddo energy_loop

    deallocate(Greenf,stat=ierr)
    call utils_dealloc_check('etrans_calculate','greenf',ierr)
    deallocate(self,stat=ierr)
    call utils_dealloc_check('etrans_calculate','self',ierr)

    call comms_barrier

    ! Communicate data
    do inode = 0, pub_total_num_nodes-1
       if (inode == pub_root_node_id) then
          cycle
       elseif (inode == pub_my_node_id) then
          call comms_send(pub_root_node_id,trc(etrc%nepptr(inode+1)+1: &
             etrc%nepptr(inode+1)+etrc%nepnode(inode+1),1:nspin), &
             etrc%nepnode(inode+1)*nspin,pub_my_node_id)
          call comms_send(pub_root_node_id,dos(etrc%nepptr(inode+1)+1: &
             etrc%nepptr(inode+1)+etrc%nepnode(inode+1),1:nspin), &
             etrc%nepnode(inode+1)*nspin,pub_my_node_id)
       elseif (pub_on_root) then
          call comms_recv(inode,trc(etrc%nepptr(inode+1)+1: &
             etrc%nepptr(inode+1)+etrc%nepnode(inode+1),1:nspin), &
             etrc%nepnode(inode+1)*nspin, inode)
          call comms_recv(inode,dos(etrc%nepptr(inode+1)+1: &
             etrc%nepptr(inode+1)+etrc%nepnode(inode+1),1:nspin), &
             etrc%nepnode(inode+1)*nspin, inode)
       endif
    enddo
    call comms_barrier

  
    ! Write transmission coefficients to file 
    if (pub_on_root) then
          output_unit = utils_unit()
          write(filename,*) trim(adjustl(pub_rootname))//'.TRC'
          open(unit=output_unit, form="formatted", &
               file=filename, action="write")

          do ienergy = 1, etrc%nep
             write(output_unit,*) real(etrc%ep(ienergy))*HARTREE_IN_EVS, &
                (trc(ienergy,is),is=1,nspin), (dos(ienergy,is),is=1,nspin)
          enddo
          close(output_unit)
    end if


    deallocate(trc,stat=ierr)
    call utils_dealloc_check('etrans_calculate','trc',ierr)
    deallocate(dos,stat=ierr)
    call utils_dealloc_check('etrans_calculate','dos',ierr)

    call timer_clock('etrans_calculate',2)
 
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving etrans_calculate'
#endif

    return

  end subroutine etrans_calculate

!====================================================================!
!====================================================================!

  subroutine transport_init_etrans_setup(ham,rep,ngwf_basis,elements)

    !================================================================!
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in November 2011                 !
    !================================================================!

    use comms, only: pub_on_root, comms_bcast, comms_barrier, &
         pub_total_num_nodes, comms_abort
    use constants, only: dp, stdout, stderr, hartree_in_evs
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
         dense_eigensolve
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use parallel_strategy, only: pub_distr_atom
    use rundat, only: pub_etrans_same_leads, pub_etrans_bulk, &
         pub_etrans_write_setup, pub_rootname, &
         pub_etrans_source, pub_etrans_drain
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_get_block
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(in)   :: ham
    type(NGWF_REP), intent(in)   :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(ELEMENT), intent(in)    :: elements(pub_cell%nat)

    ! Local variables
    type(DEM)     :: eigs_dens          ! dense matrix for eigenvectors
    type(DEM)     :: hamiltonian_dens   ! dense matrix for hamiltonian and density kernel
    type(DEM)     :: overlap_dens       ! dense matrix for overlap
    real(kind=DP) :: efermi             ! Fermi energy
    real(kind=DP), allocatable, dimension(:,:) :: eigen_en     ! hamiltonian eigenvalues

    integer    :: nspin, num, norb, orb_count
    integer    :: iat, jat, ilead, jlead, is
    integer    :: iostart, iostop, jostart, jostop
    integer    :: ierr
    character(len=80) :: filename   
 
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering transport_init_etrans_setup'
#endif

    ! Perliminaries
    nspin = pub_cell%num_spins
    isource = pub_etrans_source
    idrain = pub_etrans_drain

    !--------------------------------------------------------------------------!
    ! Identify the range of atoms/orbitals corresponding to the device setup
    !--------------------------------------------------------------------------!
    orb_count = 1
    do iat = 1, pub_cell%nat

       ! total device
       if (iat == info_device%atms(1)) info_device%orbs(1) = orb_count
       if (iat == info_device%atms(2)) info_device%orbs(2) = orb_count + &
           ngwf_basis%num_on_atom(pub_distr_atom(iat)) - 1

       ! indiviual leads
       do ilead = 1, pub_nleads
          do jlead = 1, 2
             if (iat == info_leads(ilead)%atms(1)) info_leads(ilead)%orbs(1) = orb_count
             if (iat == info_leads(ilead)%atms(2)) info_leads(ilead)%orbs(2) = orb_count + &
                 ngwf_basis%num_on_atom(pub_distr_atom(iat)) - 1
             if (iat == info_leads(ilead)%atms(3)) info_leads(ilead)%orbs(3) = orb_count
             if (iat == info_leads(ilead)%atms(4)) info_leads(ilead)%orbs(4) = orb_count + &
                 ngwf_basis%num_on_atom(pub_distr_atom(iat)) - 1
          enddo
       enddo

       orb_count = orb_count + ngwf_basis%num_on_atom(pub_distr_atom(iat))
    enddo


    ! Allocate dense matrices
    allocate(hs_leads(pub_nleads),stat=ierr)
    call utils_alloc_check('transport_init_etrans_setup','hs_leads',ierr)

    do ilead = 1, pub_nleads
       hs_leads(ilead)%norb = info_leads(ilead)%orbs(2) - info_leads(ilead)%orbs(1) + 1
    enddo

    hs_device%norb = info_device%orbs(2) - info_device%orbs(1) + 1

    do ilead = 1, pub_nleads
       norb = hs_leads(ilead)%norb
       allocate(hs_leads(ilead)%h00(nspin),stat=ierr)
       call utils_alloc_check('transport_init_etrans_setup','leads%h00',ierr)
       allocate(hs_leads(ilead)%h01(nspin),stat=ierr)
       call utils_alloc_check('transport_init_etrans_setup','leads%h01',ierr)
      
       call dense_create(hs_leads(ilead)%s00,norb,norb)
       call dense_create(hs_leads(ilead)%s01,norb,norb)
       do is = 1, nspin
          call dense_create(hs_leads(ilead)%h00(is),norb,norb)
          call dense_create(hs_leads(ilead)%h01(is),norb,norb)
       enddo
    enddo

    norb = hs_device%norb
    allocate(hs_device%h00(nspin),stat=ierr)
    call utils_alloc_check('transport_init_etrans_setup','device%h00',ierr)
    allocate(hs_device%h01(nspin),stat=ierr)
    call utils_alloc_check('transport_init_etrans_setup','device%h01',ierr)
    
    call dense_create(hs_device%s00,norb,norb)
    call dense_create(hs_device%s01,norb,norb)
    do is = 1, nspin
       call dense_create(hs_device%h00(is),norb,norb)
       call dense_create(hs_device%h01(is),norb,norb)
    enddo

    !--------------------------------------------------------------------------!
    ! Fill in dense matrices associated with device
    !--------------------------------------------------------------------------!
    num = ngwf_basis%num
    allocate(eigen_en(num,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('transport_init_etrans_setup','eigen_en',ierr)

    call dense_create(hamiltonian_dens,num,num)
    call dense_create(overlap_dens,num,num)
    call dense_create(eigs_dens,num,num)

    iostart = info_device%orbs(1)
    iostop = info_device%orbs(2)
    do is = 1, pub_cell%num_spins
       call internal_get_submatrix(1,pub_cell%nat,1,pub_cell%nat, &
               rep%overlap,overlap_dens%dmtx,num,num,ngwf_basis)
       hs_device%s00%dmtx = overlap_dens%dmtx(iostart:iostop,iostart:iostop)
       call internal_get_submatrix(1,pub_cell%nat,1,pub_cell%nat, &
               ham%ham(is),hamiltonian_dens%dmtx,num,num,ngwf_basis)
       hs_device%h00(is)%dmtx = hamiltonian_dens%dmtx(iostart:iostop,iostart:iostop)
       call dense_eigensolve(num,eigen_en(:,is),hamiltonian_dens, &
            overlap_dens,1,eigs_dens)
    enddo

    call dense_destroy(hamiltonian_dens)
    call dense_destroy(overlap_dens)
    call dense_destroy(eigs_dens)

    !--------------------------------------------------------------------------!
    ! Compute Fermi energy
    !--------------------------------------------------------------------------!
    efermi = 0.0_DP
    do is = 1, pub_cell%num_spins
       if (rep%n_occ(is)<ngwf_basis%num.and.rep%n_occ(is)>0) then
          efermi = efermi + 0.5_DP*(eigen_en(rep%n_occ(is)+1,is)+ &
               eigen_en(rep%n_occ(is),is))
       else if (rep%n_occ(is) > 0) then
          efermi = efermi + eigen_en(rep%n_occ(is),is) + tiny(1.0_DP)
       else
          efermi = efermi + eigen_en(1,is) + tiny(1.0_DP)
       end if
    end do
    info_device%efermi = efermi / pub_cell%num_spins

    deallocate(eigen_en,stat=ierr)
    call utils_dealloc_check('transport_init_etrans_setup','eigen_en',ierr)

    !--------------------------------------------------------------------------!
    ! Fill in dense matrices associated with leads
    !--------------------------------------------------------------------------!
    do ilead = 1, pub_nleads

       iostart = info_leads(ilead)%orbs(1) - info_device%orbs(1) + 1 
       iostop = info_leads(ilead)%orbs(2) - info_device%orbs(1) + 1 
       jostart = info_leads(ilead)%orbs(3) - info_device%orbs(1) + 1 
       jostop = info_leads(ilead)%orbs(4) - info_device%orbs(1) + 1 


       if (iostart .lt. jostart) info_leads(ilead)%type = 'L'
       if (iostart .gt. jostart) info_leads(ilead)%type = 'R'
       hs_leads(ilead)%s00%dmtx = hs_device%s00%dmtx(iostart:iostop,iostart:iostop)
       hs_leads(ilead)%s01%dmtx = hs_device%s00%dmtx(iostart:iostop,jostart:jostop)

       do is = 1, pub_cell%num_spins
          hs_leads(ilead)%h00(is)%dmtx = &
                 hs_device%h00(is)%dmtx(iostart:iostop,iostart:iostop)
          hs_leads(ilead)%h01(is)%dmtx = &
                 hs_device%h00(is)%dmtx(iostart:iostop,jostart:jostop)
       enddo

       info_leads(ilead)%efermi = info_device%efermi

    enddo
    
    !--------------------------------------------------------------------------!
    ! Get rid of the periodic boundary conditions
    !--------------------------------------------------------------------------!
    do ilead = 1, pub_nleads

       iostart = info_leads(ilead)%orbs(1) - info_device%orbs(1) + 1 
       iostop = info_leads(ilead)%orbs(2) - info_device%orbs(1) + 1 

       do jlead = 1, pub_nleads

          if (ilead == jlead) cycle
          jostart = info_leads(jlead)%orbs(1) - info_device%orbs(1) + 1 
          jostop = info_leads(jlead)%orbs(2) - info_device%orbs(1) + 1 

          hs_device%s00%dmtx(iostart:iostop,jostart:jostop) = 0.0_dp
          do is = 1, nspin
             hs_device%h00(is)%dmtx(iostart:iostop,jostart:jostop) = 0.0_dp
          enddo

       enddo
    enddo

    !--------------------------------------------------------------------------!
    ! Write info to output file
    !--------------------------------------------------------------------------!
    if (pub_on_root) then
       do ilead = 1, pub_nleads
          write(stdout,'(/,a,i2.2,a)') 'Lead (',ilead,') :'
          write(stdout,'(a,i5.5,a,i5.5,a)') &
           '     atoms [',info_leads(ilead)%atms(1),',',info_leads(ilead)%atms(2),']'
          write(stdout,'(a,i5.5,a,i5.5,a)') &
           '  orbitals [',info_leads(ilead)%orbs(1),',',info_leads(ilead)%orbs(2),']'
       enddo
       write(stdout,'(/,a)') 'Total device :'
       write(stdout,'(a,i5.5,a,i5.5,a)') &
        '     atoms [',info_device%atms(1),',',info_device%atms(2),']'
       write(stdout,'(a,i5.5,a,i5.5,a)') &
        '  orbitals [',info_device%orbs(1),',',info_device%orbs(2),']'
       write(stdout,'(/,a,f18.6,a)') 'Fermi Energy = ',efermi*HARTREE_IN_EVS,' eV'

       if (pub_etrans_write_setup) then
          write(filename,'(a,a)') trim(adjustl(pub_rootname)),'_device.HST'
          call write_trc_block(hs_device,'D',filename)
          do ilead = 1, pub_nleads
             write(filename,'(a,a,i1,a)') trim(adjustl(pub_rootname)),'_lead',ilead,'.HST'
             call write_trc_block(hs_leads(ilead),info_leads(ilead)%type,filename)
          enddo
       endif

    endif

    call comms_barrier

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving transport_init_etrans_setup'
#endif


  end subroutine transport_init_etrans_setup

!====================================================================!
!====================================================================!

  subroutine internal_get_submatrix(row_start,row_stop,col_start, &
                col_stop,msparse,mdense,nrow,ncol,ngwf_basis)

    !================================================================!
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in November 2011                 !
    !================================================================!

    use comms, only: pub_my_node_id, comms_bcast, pub_on_root
    use constants, only: dp, stdout
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_orig_atom, pub_distr_atom, pub_node_of_atom
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_get_block
    use dense, only: DEM
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in)          :: row_start, row_stop, nrow
    integer, intent(in)          :: col_start, col_stop, ncol
    type(SPAM3), intent(in)      :: msparse
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=dp), intent(inout)     :: mdense(nrow,ncol)

    ! Internal
    real(kind=DP),allocatable   :: buffer(:,:)
    integer    :: ierr
    integer    :: iat, jat, iorb, jorb, ispin
    integer    :: iat_dist, jat_dist
    integer    :: iat_norb, jat_norb
    integer    :: iorb_loc, jorb_loc

    iorb_loc = 0
    do iat = col_start, col_stop
       iat_dist = pub_distr_atom(iat)
       iat_norb = ngwf_basis%num_on_atom(iat_dist)

       allocate(buffer(nrow,iat_norb),stat=ierr)
       call utils_alloc_check('internal_get_submatrix','buffer',ierr)

       if (pub_my_node_id == pub_node_of_atom(iat)) then

          jorb_loc = 0
          do jat = row_start, row_stop
             jat_dist = pub_distr_atom(jat)
             jat_norb = ngwf_basis%num_on_atom(jat_dist)

             call sparse_get_block(buffer(jorb_loc+1:jorb_loc+jat_norb,:),msparse,jat_dist,iat_dist)

             jorb_loc = jorb_loc + jat_norb
          enddo

       endif

       call comms_bcast(pub_node_of_atom(iat),buffer)
       mdense(:,iorb_loc+1:iorb_loc+iat_norb) = buffer(:,:)

       deallocate(buffer,stat=ierr)
       call utils_dealloc_check('internal_get_submatrix','buffer',ierr)

       iorb_loc = iorb_loc + iat_norb
    enddo

  end subroutine internal_get_submatrix

!====================================================================!
!====================================================================!

  subroutine transport_init_energy_contour(econt,contype,eref)

    !================================================================!
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in November 2011                 !
    !================================================================!

    use comms, only: pub_on_root, comms_bcast, comms_barrier, &
         pub_total_num_nodes, comms_abort
    use constants, only: dp, stdout, stderr, hartree_in_evs
    use rundat, only: pub_etrans_ecmplx, pub_etrans_enum , &
         pub_etrans_emin, pub_etrans_emax
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none
 
    ! Arguments
    type(epath), intent(inout)   :: econt
    character(len=*), intent(in)    :: contype
    real(kind=DP), intent(in)       :: eref

    ! Local variables
    real(kind=DP) :: estart, estop, estep
    integer       :: ipt, ierr, inode, eblock

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering transport_init_energy_contour'
#endif

    if (contype .eq. 'trans') then

       econt%type = 'TR'
       econt%seg(1) = pub_etrans_enum 
       econt%seg(2:4) = 0 

       econt%nep = sum(econt%seg(:))

       allocate(econt%ep(econt%nep),stat=ierr)
       call utils_alloc_check('transport_init_energy_contour','econt%ep',ierr)
       allocate(econt%ew(econt%nep),stat=ierr)
       call utils_alloc_check('transport_init_energy_contour','econt%ew',ierr)

       estart = pub_etrans_emin + eref 
       estop  = pub_etrans_emax + eref
       estep  = (estop-estart)/(econt%nep-1)
       do ipt = 1, econt%nep
          econt%ew(ipt) = 1.0_dp
          econt%ep(ipt) = cmplx(estart+real(ipt-1,kind=DP)*estep, &
                pub_etrans_ecmplx,kind=DP)
       enddo 

       ! Distribute energy contour
       allocate(econt%nepnode(pub_total_num_nodes),stat=ierr)
       call utils_alloc_check('transport_init_energy_contour','econt%nepnode',ierr)
       allocate(econt%nepptr(pub_total_num_nodes),stat=ierr)
       call utils_alloc_check('transport_init_energy_contour','econt%nepnode',ierr)

       eblock = ceiling(real(econt%nep)/pub_total_num_nodes)
       do inode =  1, pub_total_num_nodes
          econt%nepptr(inode) = (inode-1)*eblock
          econt%nepnode(inode) = min(eblock,econt%nep-econt%nepptr(inode))
       enddo
    
       ! Write this information into the main output file 
       if (pub_on_root) then
          write(stdout,'(/,a,f8.4,a,f8.4,a)') 'Energy range for transmission calculation = [', &
               estart*hartree_in_evs,',',estop*hartree_in_evs,'] eV'
       endif

    endif

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving transport_init_energy_contour'
#endif

  end subroutine transport_init_energy_contour

!====================================================================!
!====================================================================!

  subroutine transport_greenf_calculate(energy,greenf,self,ioname)

    !================================================================!
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in April 2010                    !
    !================================================================!

    use constants, only: dp, pi, two_pi, stdout, stderr, &
         cmplx_0, cmplx_1, cmplx_i
    use comms, only: pub_my_node_id, pub_on_root, comms_abort
    use dense, only: DEM
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit
  
    implicit none

    ! Arguments
    type(DEM), intent(inout) :: greenf(pub_cell%num_spins)
    type(DEM), intent(inout) :: self(pub_cell%num_spins)
    complex(kind=dp), intent(in) :: energy
    character(len=80), intent(in) :: ioname

    ! Local variables
    complex(kind=dp), allocatable :: inv_greenf(:,:)
    integer, allocatable :: ipiv(:)

    integer :: norb
    integer :: ierr, is, iorb 

    call timer_clock('greenf_calculate',1)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering transport_greenf_calculate'
#endif

    ! Compute the self energies
    call compute_self(energy,self)

    ! Compose inv_greenf and compute greenf
    norb = hs_device%norb
    allocate(inv_greenf(norb,norb),stat=ierr)
    call utils_alloc_check('greenf_calculate','inv_greenf',ierr)

    do is = 1, pub_cell%num_spins
       inv_greenf = energy*cmplx(hs_device%s00%dmtx(:,:),kind=dp)       &
                            -cmplx(hs_device%h00(is)%dmtx(:,:),kind=dp)
       inv_greenf = inv_greenf + self(is)%zmtx(:,:)          

       greenf(is)%zmtx = cmplx_0
       do iorb = 1, norb
          greenf(is)%zmtx(iorb,iorb) = cmplx_1
       enddo

       allocate(ipiv(norb),stat=ierr)
       call utils_alloc_check('greenf_calculate','ipiv',ierr)
       
       call ZGESV(norb,norb,inv_greenf,norb,ipiv,greenf(is)%zmtx,norb,ierr)
       if (ierr .ne. 0) then
          write(stderr,*) 'Error in trc_calculate:  '
          write(stderr,*) 'Failure of zgesv with info=', ierr,' !'
          call comms_abort
       end if

       deallocate(ipiv,stat=ierr)
       call utils_dealloc_check('greenf_calculate','ipiv',ierr)
    enddo

    deallocate(inv_greenf,stat=ierr)
    call utils_dealloc_check('greenf_calculate','inv_greenf',ierr)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving transport_greenf_calculate'
#endif

    call timer_clock('greenf_calculate',2)

  end subroutine transport_greenf_calculate


!====================================================================!
!====================================================================!

  subroutine transport_trc_calculate(self,greenf,trc)

    !================================================================!
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in April 2010                    !
    !================================================================!

    use constants, only: dp, pi, two_pi, stdout, stderr, cmplx_0, &
         cmplx_1, cmplx_i
    use comms, only: pub_my_node_id, pub_on_root, comms_barrier, &
         pub_total_num_nodes, pub_root_node_id, comms_send, comms_recv, &
         comms_abort
    use dense, only: DEM
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit
  
    implicit none

    ! Arguments
    type(DEM), intent(in) :: greenf
    type(DEM), intent(in) :: self
    real(kind=dp), intent(out)  :: trc

    ! Local variables
    complex(kind=dp), allocatable :: gammal(:,:), gammar(:,:)
    complex(kind=dp), allocatable :: trcl(:,:), trcr(:,:)
    integer :: norb, norbr, norbl
    integer :: lshift, rshift
    integer :: iorb, jorb, ierr


#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering transport_trc_calculate'
#endif

    norb = hs_device%norb
    norbl = hs_leads(isource)%norb
    norbr = hs_leads(idrain)%norb
    lshift = info_leads(isource)%orbs(1) - info_device%orbs(1)
    rshift = info_leads(idrain)%orbs(1) - info_device%orbs(1)

    allocate(gammal(norbl,norbl),stat=ierr)
    call utils_alloc_check('transport_trc_calculate','gammal',ierr)
    allocate(gammar(norbr,norbr),stat=ierr)
    call utils_alloc_check('transport_trc_calculate','gammar',ierr)
    allocate(trcl(norbl,norbr),stat=ierr)
    call utils_alloc_check('transport_trc_calculate','trcl',ierr)
    allocate(trcr(norbr,norbl),stat=ierr)
    call utils_alloc_check('transport_trc_calculate','trcr',ierr)

    gammal = cmplx_i*(self%zmtx(lshift+1:lshift+norbl,lshift+1:lshift+norbl) - &
             dconjg(transpose(self%zmtx(lshift+1:lshift+norbl,lshift+1:lshift+norbl))))
    gammar = cmplx_i*(self%zmtx(rshift+1:rshift+norbr,rshift+1:rshift+norbr) - &
             dconjg(transpose(self%zmtx(rshift+1:rshift+norbr,rshift+1:rshift+norbr))))

    ! Calculate the transmission matrix 
    call zgemm('N','N',norbl,norbr,norbl,cmplx_1,gammal,norbl, &
               greenf%zmtx(lshift+1:lshift+norbl,rshift+1:rshift+norbr), &
               norbr,cmplx_0,trcl,norbl)
    call zgemm('N','C',norbr,norbl,norbr,cmplx_1,gammar,norbr, &
               greenf%zmtx(lshift+1:lshift+norbl,rshift+1:rshift+norbr), &
               norbl,cmplx_0,trcr,norbr)

    trc = 0.0_dp
    do iorb=1,norbl
       do jorb=1,norbr
          trc = trc + real(trcl(iorb,jorb)*trcr(jorb,iorb))
      enddo
    enddo

    deallocate(gammal,stat=ierr)
    call utils_alloc_check('transport_trc_calculate','gammal',ierr)
    deallocate(gammar,stat=ierr)
    call utils_alloc_check('transport_trc_calculate','gammar',ierr)
    deallocate(trcl,stat=ierr)
    call utils_alloc_check('transport_trc_calculate','trcl',ierr)
    deallocate(trcr,stat=ierr)
    call utils_alloc_check('transport_trc_calculate','trcr',ierr)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving transport_trc_calculate'
#endif

  end subroutine transport_trc_calculate

!====================================================================!
!====================================================================!


    subroutine compute_self(energy_cmp,self)

    !================================================================!
    !                                                                !
    ! This subroutine compute the lead self-energies by means of     !
    ! iterative improvements of the bulk transfer matrix             !
    ! see Lopez-Sancho et al., J. Phys. F: Met.Phys. 15, 851 (1985)  !
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in April 2010                    !
    !================================================================!

    use constants, only: dp, pi, stdout, stderr, cmplx_0, cmplx_1, &
           cmplx_i
    use comms, only: comms_abort, pub_on_root
    use rundat, only: pub_etrans_same_leads, pub_etrans_bulk
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit

    implicit none

    ! Arguments
    complex(kind=dp), intent(in) :: energy_cmp
    type(DEM), intent(inout)     :: self(pub_cell%num_spins)

    ! Internal variables
    complex(kind=dp), allocatable   :: tmp1(:,:)
    complex(kind=dp), allocatable   :: transf(:,:), transf_bar(:,:)
    complex(kind=dp), allocatable   :: self1(:,:), self2(:,:)

    integer  :: norb, shift
    integer  :: is, il
    integer  :: ierr

    real(kind=dp), parameter    :: conv_tol=1.0e-7_dp

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering compute_self'
#endif

    !====================================================!
    ! Compute : Self_L = (e*S01-H01)^{dag} * TRANSF_bar  !
    !           Self_R = (e*S01-H01) * TRANSF            !
    !====================================================!

    do is = 1, pub_cell%num_spins

       self(is)%zmtx = cmplx_0

       if (pub_etrans_same_leads .or. pub_etrans_bulk) then

          norb = hs_leads(1)%norb
          allocate(transf(norb,norb),stat=ierr)
          call utils_alloc_check('compute_self','transf',ierr)
          allocate(transf_bar(norb,norb),stat=ierr)
          call utils_alloc_check('compute_self','transf_bar',ierr)
          allocate(self1(norb,norb),stat=ierr)
          call utils_alloc_check('compute_self','self1',ierr)
          allocate(self2(norb,norb),stat=ierr)
          call utils_alloc_check('compute_self','self2',ierr)
          
          ! Compute the transfer matrices
          call compute_transfer(transf,transf_bar,hs_leads(1)%h00(is)%dmtx, &
                  hs_leads(1)%h01(is)%dmtx,hs_leads(1)%s00%dmtx, &
                  hs_leads(1)%s01%dmtx,energy_cmp,norb,conv_tol,50)

          ! Compute the left and right self energies
          allocate(tmp1(norb,norb),stat=ierr)
          call utils_alloc_check('compute_self','tmp1',ierr)

          tmp1 =  energy_cmp*cmplx(hs_leads(1)%s01%dmtx,kind=dp) &
                            -cmplx(hs_leads(1)%h01(is)%dmtx,kind=dp)
         
          call ZGEMM('C','N',norb,norb,norb,cmplx_1,tmp1, &
                     norb,transf_bar,norb,cmplx_0,self1,norb)
          
          call ZGEMM('N','N',norb,norb,norb,cmplx_1,tmp1, &
                     norb,transf,norb,cmplx_0,self2,norb)

          deallocate(tmp1,stat=ierr)
          call utils_dealloc_check('compute_self','tmp1',ierr)

          do il = 1, pub_nleads
             shift = info_leads(il)%orbs(1) - info_device%orbs(1)
             if (info_leads(il)%type .eq. info_leads(1)%type) then
                self(is)%zmtx(shift+1:shift+norb,shift+1:shift+norb) = self1
             else
                self(is)%zmtx(shift+1:shift+norb,shift+1:shift+norb) = self2
             endif
          enddo      
                     
          deallocate(transf,stat=ierr)
          call utils_dealloc_check('compute_self','transf',ierr)
          deallocate(transf_bar,stat=ierr)
          call utils_dealloc_check('compute_self','transf_bar',ierr)
          deallocate(self1,stat=ierr)
          call utils_dealloc_check('compute_self','self1',ierr)
          deallocate(self2,stat=ierr)
          call utils_dealloc_check('compute_self','self2',ierr)

       else 

          do il = 1, pub_nleads

             norb = hs_leads(il)%norb
             allocate(transf(norb,norb),stat=ierr)
             call utils_alloc_check('compute_self','transf',ierr)
             allocate(transf_bar(norb,norb),stat=ierr)
             call utils_alloc_check('compute_self','transf_bar',ierr)
             allocate(self1(norb,norb),stat=ierr)
             call utils_alloc_check('compute_self','self1',ierr)
             
             ! Compute the transfer matrices
             call compute_transfer(transf,transf_bar,hs_leads(il)%h00(is)%dmtx, &
                     hs_leads(il)%h01(is)%dmtx,hs_leads(il)%s00%dmtx, &
                     hs_leads(il)%s01%dmtx,energy_cmp,norb,conv_tol,50)

             ! Compute the left and right self energies
             allocate(tmp1(norb,norb),stat=ierr)
             call utils_alloc_check('compute_self','tmp1',ierr)
             
             tmp1 =  energy_cmp*cmplx(hs_leads(il)%s01%dmtx,kind=dp) &
                               -cmplx(hs_leads(il)%h01(is)%dmtx,kind=dp)
        
             call ZGEMM('C','N',norb,norb,norb,cmplx_1,tmp1,norb, &
                    transf_bar,norb,cmplx_0,self1,norb)
             
             deallocate(tmp1,stat=ierr)
             call utils_dealloc_check('compute_self','tmp1',ierr)

             shift = info_leads(il)%orbs(1) - info_device%orbs(1)
             self(is)%zmtx(shift+1:shift+norb,shift+1:shift+norb) = self1

             deallocate(transf,stat=ierr)
             call utils_dealloc_check('compute_self','transf',ierr)
             deallocate(transf_bar,stat=ierr)
             call utils_dealloc_check('compute_self','transf_bar',ierr)
             deallocate(self1,stat=ierr)
             call utils_dealloc_check('compute_self','self1',ierr)

          enddo

       endif
    enddo

    return

  end subroutine compute_self


!====================================================================!
!====================================================================!

  subroutine compute_transfer(transf,transf_bar,h_00,h_01,s_00,s_01, &
                              energy,norb,convtol,maxiter) 

    !================================================================!
    !                                                                !
    ! This subroutine compute the lead self-energies by means of     !
    ! iterative improvements of the bulk transfer matrix             !
    ! see Lopez-Sancho et al., J. Phys. F: Met.Phys. 15, 851 (1985)  !
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in April 2010                    !
    !================================================================!

    use constants, only: dp, cmplx_0, cmplx_1, cmplx_i, stdout, stderr
    use comms, only: pub_on_root, comms_abort, pub_my_node_id


    implicit none

    ! Arguments
    integer, intent(in)           :: norb
    integer, intent(in)           :: maxiter
    complex(kind=dp), intent(in)  ::  energy
    complex(kind=dp), intent(inout) ::  transf(norb,norb)
    complex(kind=dp), intent(inout) ::  transf_bar(norb,norb)
    real(kind=dp), intent(in)     :: h_00(norb,norb)
    real(kind=dp), intent(in)     :: h_01(norb,norb)
    real(kind=dp), intent(in)     :: s_00(norb,norb)
    real(kind=dp), intent(in)     :: s_01(norb,norb)
    real(kind=dp), intent(in)     :: convtol

    ! Internal variables
    real(kind=dp)    :: conver,conver2,transfnorm,transf_barnorm
    complex(kind=dp) :: tprod(norb,norb), tprod_bar(norb,norb)
    complex(kind=dp) :: t0(norb,norb), t0_bar(norb,norb)
    complex(kind=dp) :: t1(norb,norb), t1_bar(norb,norb)
    complex(kind=dp) :: tmp1(norb,norb), tmp2(norb,norb), tmp3(norb,norb)

    integer  :: norb2, ierr, info
    integer  :: iorb, jorb, iter
    integer  :: ipiv(norb)
    
    !================================================================!
    !=== Compute t0 and t0_bar :
    !===   t0 =  k00^{-1} * k01^{dag}
    !===   t0_bar =  k00^{-1} * k01
    !================================================================!

    norb2 = norb*norb

    ! compute k00=(e*s00 - h00) --> tmp2
    !         k01=-(e*s01 - h01) --> tmp3
    tmp2(:,:) = energy*cmplx(s_00(:,:),kind=dp) - cmplx(h_00(:,:),kind=dp)
    tmp3(:,:) = cmplx(h_01(:,:),kind=dp) - energy*cmplx(s_01(:,:),kind=dp)

    ! invert k00 --> tmp1
    tmp1 = cmplx_0
    do iorb=1,norb
       tmp1(iorb,iorb) = cmplx_1
    end do
    
    call ZGESV(norb,norb,tmp2,norb,ipiv,tmp1,norb,info)
    if (info.ne.0) then
       write(stderr,*) 'Error in compute_transfer:  '
       write(stderr,*) 'Failure of zgesv with info=', info,' !'
       call comms_abort
    end if
    
    ! compute the t0 and t0_bar matrices
    t0 = cmplx_0
    t0_bar = cmplx_0

    call ZGEMM('N','C',norb,norb,norb,cmplx_1,tmp1,norb,tmp3,norb,cmplx_0,t0,norb)
    call ZGEMM('N','N',norb,norb,norb,cmplx_1,tmp1,norb,tmp3,norb,cmplx_0,t0_bar,norb)

    ! Initialization of the bulk transfer matrices 
    transf(:,:) = t0(:,:)
    tprod(:,:) = t0_bar(:,:)
    transf_bar(:,:) = t0_bar(:,:)
    tprod_bar(:,:) = t0(:,:)


    !================================================================!
    !=== Main iteration loop
    !===   transf(i) = transf(i-1) + tprod(i-1)*t(i) 
    !===   transf_bar(i) = transf_bar(i-1) + tprod_bar(i-1)*t(i) 
    !================================================================!
    
    do iter=1,maxiter
                
       !======================================!
       !=== Compute t1=t(i+1) from t0=t(i) ===!
       !======================================!
       call ZGEMM('N','N',norb,norb,norb,cmplx_1,t0,norb,t0_bar,norb,cmplx_0,tmp1,norb)
       call ZGEMM('N','N',norb,norb,norb,cmplx_1,t0_bar,norb,t0,norb,cmplx_0,tmp2,norb)

       tmp3(:,:) = -tmp1(:,:)-tmp2(:,:)
       do iorb=1,norb
          tmp3(iorb,iorb) = cmplx_1 + tmp3(iorb,iorb)
       end do

       tmp1 = cmplx_0
       do iorb=1,norb
          tmp1(iorb,iorb)=cmplx_1 
       end do

       call ZGESV(norb,norb,tmp3,norb,ipiv,tmp1,norb,info)
       if (info.ne.0) then
          write(stderr,*) 'Error in compute_transfer:  '
          write(stderr,*) 'Failure of zgesv with info=', info,' !'
          call comms_abort
       end if

       call ZGEMM('N','N',norb,norb,norb,cmplx_1,t0,norb,t0,norb,cmplx_0,tmp2,norb)
       call ZGEMM('N','N',norb,norb,norb,cmplx_1,t0_bar,norb,t0_bar,norb,cmplx_0,tmp3,norb)
      
       call ZGEMM('N','N',norb,norb,norb,cmplx_1,tmp1,norb,tmp2,norb,cmplx_0,t1,norb)
       call ZGEMM('N','N',norb,norb,norb,cmplx_1,tmp1,norb,tmp3,norb,cmplx_0,t1_bar,norb)

       !===============================================!
       !=== Compute transf(i+1) and transf_bar(i+1) ===!
       !===============================================!
       call ZGEMM('N','N',norb,norb,norb,cmplx_1,tprod,norb,t1,norb,cmplx_0,tmp1,norb)
       call ZGEMM('N','N',norb,norb,norb,cmplx_1,tprod,norb,t1_bar,norb,cmplx_0,tmp2,norb)
        
       call ZAXPY(norb2,cmplx_1,tmp1,1,transf,1)
       tprod = tmp2

       call ZGEMM('N','N',norb,norb,norb,cmplx_1,tprod_bar,norb,t1_bar,norb,cmplx_0,tmp1,norb)
       call ZGEMM('N','N',norb,norb,norb,cmplx_1,tprod_bar,norb,t1,norb,cmplx_0,tmp2,norb) 
       
       call ZAXPY(norb2,cmplx_1,tmp1,1,transf_bar,1)
       tprod_bar = tmp2

       !=================!
       !=== Update t0 ===!
       !=================!
       t0 = t1
       t0_bar = t1_bar

       !=============================!
       !=== Check the convergence ===!
       !=============================!
       conver = 0.0_dp
       conver2 = 0.0_dp

       do jorb=1,norb
          do iorb=1,norb
              conver=conver+sqrt(real(t1(iorb,jorb),dp)**2+aimag(t1(iorb,jorb))**2)
              conver2=conver2+sqrt(real(t1_bar(iorb,jorb),dp)**2+aimag(t1_bar(iorb,jorb))**2)
          end do
       end do

       if (conver.lt.convtol .and. conver2.lt.convtol) exit

    end do 

    return 

  end subroutine compute_transfer

!====================================================================!
!====================================================================!

  subroutine write_trc_block(hs,block_type,filename)

    use comms, only: pub_on_root, comms_bcast, comms_barrier, &
         pub_total_num_nodes, comms_abort
    use constants, only: dp, stdout, stderr
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit

    implicit none

    ! Arguments
    type(block_hs), intent(in)   :: hs
    character(len=1), intent(in) :: block_type
    character(len=*), intent(in) :: filename

    ! Local variables
    integer    :: io_unit, io, jo, is

    io_unit = utils_unit()
    open(unit=io_unit,file=filename,form='unformatted')

    write(io_unit) pub_cell%num_spins
    write(io_unit) hs%norb
    write(io_unit) block_type
    
    do is = 1, pub_cell%num_spins
       write(io_unit) ((hs%h00(is)%dmtx(io,jo),io=1,hs%norb),jo=1,hs%norb)
    enddo
    if (block_type .ne. 'D') then
       do is = 1, pub_cell%num_spins
          write(io_unit) ((hs%h01(is)%dmtx(io,jo),io=1,hs%norb),jo=1,hs%norb)
       enddo
    endif

    write(io_unit) ((hs%s00%dmtx(io,jo),io=1,hs%norb),jo=1,hs%norb)
    if (block_type .ne. 'D') then
       write(io_unit) ((hs%s01%dmtx(io,jo),io=1,hs%norb),jo=1,hs%norb)
    endif

    close(io_unit)

  end subroutine write_trc_block

!==========================================================================!
  
end module transport
