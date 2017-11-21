! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=================================================================!
!                    Coulomb Cutoff module                        !
!                                                                 !
! This module contains subroutines for initialisation and de-     !
! initialisation of the cutoff Coulomb interaction, and           !
! wrappers for the routines which need to work on different cells !
! or otherwise be treated differently in the presence of cutoff   !
! Coulomb interactions.                                           !
!-----------------------------------------------------------------!
! Written by Nicholas Hine, March 2008                            !
!=================================================================!

module cutoff_coulomb

  use cell_grid, only: GRID_INFO
  use constants, only: DP, PI, stdout, VERBOSE
  use simulation_cell, only: pub_cell, pub_padded_cell

  implicit none

  ! Public routines
  public :: cutoff_coulomb_init
  public :: cutoff_coulomb_exit
  public :: cutoff_coulomb_hartree
  public :: cutoff_coulomb_localpseudo
  public :: cutoff_coulomb_core_density
  public :: cutoff_coulomb_locps_forces
  public :: cutoff_coulomb_nlcc_forces
  public :: cutoff_coulomb_initial_guess
  public :: cutoff_coulomb_struct_fac
  public :: cutoff_coulomb_classical_sfac
  public :: cutoff_coulomb_ii_energy
  public :: cutoff_coulomb_ii_forces

  ! Padded version of the fine grid
  type(GRID_INFO), public :: pub_padded_grid

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cutoff_coulomb_ii_energy(elements,ii_energy,cutoff_type)

    !=====================================================================!
    ! Calculates the direct ion-ion contribution to the energy. Replaces  !
    ! ewald_calculate_energy for finite systems                           !
    !---------------------------------------------------------------------!
    ! Arguments:                                                           !
    ! elements (input): list of all the atoms in the system                !
    ! ii_energy (output): the ion-ion energy                               !
    ! cutoff_type (input): the cutoff type to apply                        !
    !----------------------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                                 !
    ! Modified by Jacek Dziedzic in January 2010 to accept the cutoff type !
    ! in the third argument rather than relying on pub_coulomb_cutoff_type.!
    !======================================================================!

    use comms, only: comms_reduce, pub_on_root, pub_my_node_id
    use geometry, only: POINT, magnitude, OPERATOR(.cross.), OPERATOR(.dot.), &
         OPERATOR(-), OPERATOR(+), OPERATOR(*)
    use ion, only : ELEMENT
    use parallel_strategy, only: pub_first_atom_on_node, &
         pub_num_atoms_on_node
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(ELEMENT),intent(in) :: elements(pub_cell%nat)
    real(dp),intent(out) :: ii_energy
    character(len=*),intent(in) :: cutoff_type

    ! Locals
    integer :: my_first_at, my_last_at  ! range of atoms on this node
    integer :: iatom, jatom             ! atom index counters
    real(dp) :: qi, qj                  ! charge of atoms i and j
    !real(dp) :: minv(3,3)     ! inverse of peridocity-related coordinate system
    type(POINT) :: a1, a2, a3 ! lattice (unit) vectors of periodicity
    type(POINT) :: b1, b2     ! reciprocal lattice vectors of 2D in-plane cell
    type(POINT) :: rij        ! distance between particles i and j
    real(dp) :: self_term     ! Ewald self term for periodic images

    ! 1D Ewald constants
    integer, parameter :: mmax = 6
    real(dp) :: alat, U       ! cell length in 1D systems and (mmax+0.5)*alat

    ! 2D Ewald constants
    real(dp) :: Rmax, Gmax, sigma, Acell
    integer :: n1_real, n2_real, n1_recip, n2_recip

    ii_energy = 0.0d0

    my_first_at = pub_first_atom_on_node(pub_my_node_id)
    my_last_at = my_first_at + pub_num_atoms_on_node(pub_my_node_id) - 1

    select case (cutoff_type)

    ! 0D periodicity
    case('SPHERE','CYLINDER')

       self_term = 0.0d0

       ! Loop over my atoms
       do iatom = my_first_at, my_last_at

          qi = elements(iatom)%ion_charge

          ! Loop over all other atoms
          do jatom = 1, pub_cell%nat

             if (jatom == iatom ) cycle

             ! Find charge and distance of other atom
             qj = elements(jatom)%ion_charge
             rij = elements(iatom)%centre - elements(jatom)%centre

             ! Add point charge contribution to energy
             ii_energy = ii_energy + qi*qj/magnitude(rij)

          enddo
       enddo

       ! remove double counting
       ii_energy = ii_energy * 0.50d0

    ! 1D periodicity
    case('WIRE')

       ! Set up periodicity vector
       a1 = (1.0d0/magnitude(pub_cell%a1))*pub_cell%a1
       alat = magnitude(pub_cell%a1)

       ! Set up constants determined from a1
       U = (real(mmax,DP)+0.5_DP)*alat

       ! Calculate Self Term
       self_term = ewald_self_1d()

       ! Loop over my atoms
       do iatom = my_first_at, my_last_at

          qi = elements(iatom)%ion_charge
          ii_energy = ii_energy + qi*qi*self_term

          ! Loop over other atoms with j>i
          do jatom = iatom + 1, pub_cell%nat

             if (jatom == iatom ) cycle

             ! Find charge and distance of other atom
             qj = elements(jatom)%ion_charge
             rij = elements(iatom)%centre - elements(jatom)%centre

             ! Add point charge contribution to energy
             ii_energy = ii_energy + qi*qj*ewald_1d(rij)

          enddo
       enddo

    ! 2D periodicity
    case('SLAB')

       ! Set up shorthand lattice vectors and '12' plane unit normal
       a1 = pub_cell%a1
       a2 = pub_cell%a2
       a3 = a1.cross.a2
       a3 = (1.0d0/magnitude(a3))*a3

       ! Set up reciprocal lattice vectors in '12' plane
       ! (not the same as pub_cell%b1,b2)
       b1=(2.0_DP*PI/( a1.dot.(a2.cross.a3) ) ) * (a2.cross.a3)
       b2=(2.0_DP*PI/( a1.dot.(a2.cross.a3) ) ) * (a3.cross.a1)

       ! Set up constants
       Acell = magnitude(a1.cross.a2)
       sigma = sqrt(Acell)/2.4d0
       Rmax = sqrt(10.0d0*Acell/2.4d0)
       Gmax = sqrt(40.0d0*PI*PI/Acell)
       n1_real = int(Rmax/magnitude(a1))
       n2_real = int(Rmax/magnitude(a2))
       n1_recip = int(Gmax/magnitude(b1))
       n2_recip = int(Gmax/magnitude(b2))

       ! Calculate Self Term
       self_term = ewald_self_2d()

       ! Loop over my atoms
       do iatom = my_first_at, my_last_at

          qi = elements(iatom)%ion_charge
          ii_energy = ii_energy + 0.5d0*qi*qi*self_term

          ! Loop over other atoms with j>i
          do jatom = iatom + 1, pub_cell%nat

             if (jatom == iatom ) cycle

             ! Find charge and distance of other atom
             qj = elements(jatom)%ion_charge
             rij = elements(iatom)%centre - elements(jatom)%centre

             ! Add point charge contribution to energy
             ii_energy = ii_energy + qi*qj*ewald_2d(rij)

          enddo
       enddo

    ! jd: Safety valve, saves debugging effort
    case default
       call utils_abort('Illegal cutoff type in cutoff_coulomb_ii_energy')

    end select

    call comms_reduce('SUM', ii_energy)

    if (pub_on_root) then
       write(stdout,*)
       if (self_term /= 0.0d0) then
          write(stdout,'(a,f24.14)') 'Self-interaction term : ',self_term
       end if
       write(stdout,'(a,f24.14)')    'Ion-Ion Energy        : ',ii_energy
    end if

  contains

    real(dp) function ewald_1d(rij)
      !-----------------------------------------------------------------------!
      ! Evaluate 1 dimensional Coulomb sum for particle at (0,0,0) interacting!
      ! with particle at d and its periodic images (summed over n different d)!
      ! Algorithm involves accurate approximation based on Euler-Maclaurin    !
      ! summation formula - see e.g. Eq. 4.8 in Comp.Phys.Commun. 84, 156     !
      ! (1994).                                                               !
      ! In the following the two sums are truncated at M = 6ish, and m=2 (Eq. !
      ! 4.5ab) to give an accuracy of roughly 1x10^-7.                        !
      !-----------------------------------------------------------------------!
      implicit none
      type(POINT), intent(in) :: rij

      real(dp), parameter :: E1 = -1.0_DP/24.0_DP, E2 = 7.0_DP/5760_DP
      integer, parameter :: mmax=6

      integer :: m
      real(dp) :: z, alpha, upz, umz, supz2a, sumz2a, upz2a, umz2a

      ! Set up distance in a1 dir and distance squared perp to it
      z = rij.dot.a1
      alpha = rij%X**2 + rij%Y**2 + rij%Z**2 - z**2

      ! 1/r term
      ewald_1d = 1.0d0/sqrt(z**2 + alpha)
      ! 1/(r+Rm) terms from periodic interactions
      do m = 1, mmax
         ewald_1d = ewald_1d + 1.d0/sqrt((z + alat*real(m,dp))**2+alpha)
         ewald_1d = ewald_1d + 1.d0/sqrt((z - alat*real(m,dp))**2+alpha)
      enddo

      upz = U + z
      upz2a = upz*upz + alpha
      supz2a = sqrt(upz2a)
      umz = U - z
      umz2a = umz*umz + alpha
      sumz2a = sqrt(umz2a)

      ! First Euler-MacLaurin term: (H(U-z,alpha) + H(U+z,alpha) - 2 ln(a))/a
      ewald_1d = ewald_1d - log((supz2a + upz)*(sumz2a + umz))/alat

      ! Second Euler-MacLaurin terms: xi(M,r) + xi(M,-r)
      ewald_1d = ewald_1d + (upz/upz2a/supz2a + umz/umz2a/sumz2a)*E1*alat
      ewald_1d = ewald_1d - (9.0_DP - 15.0_DP*upz*upz/upz2a) / &
           (upz2a**2*supz2a)*upz*E2*alat**3
      ewald_1d = ewald_1d - (9.0_DP - 15.0_DP*umz*umz/umz2a) / &
           (umz2a**2*sumz2a)*umz*E2*alat**3

    end function ewald_1d

    real(dp) function ewald_self_1d()

      implicit none

      ! Locals
      real(dp), parameter :: E1 = -1.0_DP/24.0_DP, E2 = 7.0_DP/5760_DP
      integer :: m

      ! 1/(Rm) terms from periodic copies of self
      ewald_self_1d = 0.d0
      do m = 1, mmax
         ewald_self_1d = ewald_self_1d + 2.0d0/(alat*real(m,dp))
      enddo

      ! First Euler-MacLaurin term: (H(U,0) + H(U,0) - 2 ln(a))/a
      ewald_self_1d = ewald_self_1d - 2.0d0*log(2.0d0*U)/alat

      ! Second Euler-MacLaurin terms: 2 xi(M,0)
      ewald_self_1d = ewald_self_1d + (2.0d0/(U**2))*E1*alat
      ewald_self_1d = ewald_self_1d + (12.0d0/(U**4))*E2*alat**3
      ewald_self_1d = ewald_self_1d*0.5d0

    end function ewald_self_1d

    real(dp) function ewald_2d(rij)
      !-----------------------------------------------------------------------!
      ! Evaluate 2 dimensional Coulomb sum for particle at (0,0,0) interacting!
      ! with particle at rij and its periodic images.                         !
      !-----------------------------------------------------------------------!

      use utils, only: utils_erf, utils_erfc

      implicit none

      ! Arguments
      type(POINT), intent(in) :: rij

      ! Locals
      type(POINT) :: Rvec, Gvec, rpar
      real(dp) :: fac1, rperp, zerf, cosGr
      real(dp) :: Rmag, Gmag
      real(dp) :: realsum, recipsum
      integer :: i1, i2, m2

      ! Set up perpendicular and parallel components of rij
      rperp = rij.dot.a3
      rpar = rij - rperp*a3

      ! Real Space Sum
      realsum = 0.0d0
      do i1 = -n1_real, n1_real
         do i2 = -n2_real, n2_real
            Rvec = real(i1,dp)*a1 + real(i2,dp)*a2
            Rmag = magnitude(Rvec)
            if (Rmag < Rmax) then
               Rmag = magnitude(rij - Rvec)
               realsum = realsum + utils_erfc(Rmag/sigma)/Rmag
            end if
         end do
      end do

      ! 2pi/A(z erf(z/sigma)+sigma/sqrt(pi)*exp(-z^2/sigma^2)) term
      zerf = (-2.0d0*PI/Acell)*(rperp*utils_erf(rperp/sigma) + &
           sigma*exp(-rperp*rperp/(sigma*sigma))/sqrt(PI))

      ! Reciprocal Space Sum
      fac1 = 2.0d0*PI/Acell  ! each term includes k and -k, hence the 2
      m2 = 1
      recipsum = 0.0d0
      do i1 = 0, n1_recip
         do i2 = m2, n2_recip
            Gvec = real(i1,dp)*b1 + real(i2,dp)*b2
            Gmag = magnitude(Gvec)
            if (Gmag <= Gmax) then
               cosGr = cos(Gvec.dot.rpar)
               recipsum = recipsum + fac1/Gmag*cosGr*( &
                    exp(-Gmag*rperp)*utils_erfc(sigma*Gmag*0.5d0-rperp/sigma) +&
                    exp(+Gmag*rperp)*utils_erfc(sigma*Gmag*0.5d0+rperp/sigma))
            end if

         end do
         ! Now go from -n2_recip to n2_recip
         m2 = -n2_recip
      end do

      ewald_2d = realsum + zerf + recipsum

    end function ewald_2d

    real(dp) function ewald_self_2d()

      use utils, only: utils_erfc

      implicit none

      ! Locals
      type(POINT) :: Rvec, Gvec
      real(dp) :: Rmag, Gmag, fac1
      real(dp) :: realsum, recipsum, recip0
      integer :: i1, i2, m2

      ! Real Space Sum
      fac1 = 2.0d0
      m2 = 1
      realsum = -2.0d0/(sigma*sqrt(PI))
      do i1 = 0, n1_real
         do i2 = m2, n2_real
            Rvec = real(i1,dp)*a1 + real(i2,dp)*a2
            Rmag = magnitude(Rvec)
            if (Rmag < Rmax) then
               realsum = realsum + fac1*utils_erfc(Rmag/sigma)/Rmag
            end if
         end do
         ! On first '1' step we only go from 1 to n2
         ! subsequently, we need to go from -n2 to n2
         m2 = -n2_real
      end do

      ! Reciprocal Space Sum
      recip0 = -2.0d0*sigma*sqrt(PI)/Acell
      fac1 = 4.0d0*PI/Acell
      m2 = 1
      recipsum = 0.0d0
      do i1 = 0, n1_recip
         do i2 = m2, n2_recip
            Gvec = real(i1,dp)*b1 + real(i2,dp)*b2
            Gmag = magnitude(Gvec)
            if (Gmag <= Gmax) then
               recipsum = recipsum + fac1*utils_erfc(Gmag*sigma*0.5D0)/Gmag
            end if

         end do
         ! Now go from -n2_recip to n2_recip
         m2 = -n2_recip
      end do

      ewald_self_2d = realsum + recipsum + recip0

    end function ewald_self_2d

  end subroutine cutoff_coulomb_ii_energy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cutoff_coulomb_ii_forces(elements,ii_forces)

    !=====================================================================!
    ! Calculates the direct ion-ion contribution to the energy. Replaces  !
    ! ewald_calculate_forces for finite systems                           !
    !---------------------------------------------------------------------!
    ! Arguments:                                                          !
    ! elements (input): list of all the atoms in the system               !
    ! ii_energy (output): all ion-ion forces                              !
    !---------------------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                                !
    !=====================================================================!

    use comms, only : comms_reduce, pub_my_node_id
    use geometry, only : POINT, OPERATOR(.cross.), OPERATOR(.dot.), &
         OPERATOR(+), OPERATOR(-), OPERATOR(*), magnitude
    use ion, only : ELEMENT
    use parallel_strategy, only: pub_first_atom_on_node, &
         pub_num_atoms_on_node
    use rundat, only : pub_coulomb_cutoff_type

    implicit none

    ! Arguments
    type(ELEMENT),intent(in) :: elements(pub_cell%nat)
    real (kind=dp), dimension(1:3,1:pub_cell%nat), intent(out) :: ii_forces

    ! Locals
    integer :: my_first_at, my_last_at  ! range of atoms on this node
    integer :: iatom, jatom             ! atom index counters
    real(dp) :: qi, qj                  ! charge of atoms i and j
    real(dp) :: fij(3)        ! force between particles i and j
    type(POINT) :: rij        ! distance between particles i and j
    real(dp) :: rij3          ! |ri-rj|^3
    type(POINT) :: a1, a2, a3 ! lattice (unit) vectors of periodicity
    type(POINT) :: a1p, a2p   ! unit vectors for 'x' and 'y' directions in 2D
    type(POINT) :: b1, b2     ! reciprocal lattice vectors of 2D in-plane cell
    real(dp) :: minv(3,3)     ! inverse of peridocity-related coordinate system

    ! 1D Ewald Constants
    integer, parameter :: mmax=6
    real(dp) :: alat, U       ! cell length in 1D systems

    ! 2D Ewald Constants
    real(dp) :: Rmax, Gmax, sigma, Acell
    integer :: n1_real, n2_real, n1_recip, m2_recip, n2_recip

    ii_forces(:,:) = 0.0d0

    my_first_at = pub_first_atom_on_node(pub_my_node_id)
    my_last_at = my_first_at + pub_num_atoms_on_node(pub_my_node_id) - 1

    select case (pub_coulomb_cutoff_type)

    ! 0D periodicity
    case('SPHERE','CYLINDER')

       ! Loop over my atoms
       do iatom = my_first_at, my_last_at

          qi = real(elements(iatom)%ion_charge, kind=DP)

          ! Loop over all other atoms j<i
          do jatom = 1, iatom - 1

             ! Find charge and distance of other atom
             qj = real(elements(jatom)%ion_charge,kind=DP)
             rij = elements(iatom)%centre - elements(jatom)%centre
             rij3 = magnitude(rij)**3

             ! Find force between atoms i and j
             fij(1) = (qi*qj/rij3)*rij%X
             fij(2) = (qi*qj/rij3)*rij%Y
             fij(3) = (qi*qj/rij3)*rij%Z

             ! Add contribution of jatom to force on iatom
             ii_forces(:,iatom) = ii_forces(:,iatom) + fij

             ! Add contribution of iatom to force on jatom
             ii_forces(:,jatom) = ii_forces(:,jatom) - fij

          enddo

       enddo

    ! 1D periodicity
    case('WIRE')

       ! Set up periodicity direction vector (becomes 'z')
       alat = magnitude(pub_cell%a1)
       a1 = (1.0d0/alat)*pub_cell%a1

       ! Set up an arbitrary perpendicular unit vector to be 'x'
       a2%X = a1%Y + a1%Z
       a2%Y = -a1%X + a1%Z
       a2%Z = -a1%X - a1%Y
       a2 = (1.0d0/magnitude(a2))*a2

       ! Set up a vector perpendicular to both to be 'y'
       a3 = a1.cross.a2
       a3 = (1.0d0/magnitude(a3))*a3

       ! Set up the inverse of the matrix comprising a1, a2, a3 to enable
       ! us to go back from coordinates relative to the periodicity to the
       ! simulation cell coordinates
       minv(1,1) = a2%X; minv(2,1) = a2%Y; minv(3,1) = a2%Z
       minv(1,2) = a3%X; minv(2,2) = a3%Y; minv(3,2) = a3%Z
       minv(1,3) = a1%X; minv(2,3) = a1%Y; minv(3,3) = a1%Z

       ! Set up constants determined from a1
       U = (real(mmax,DP)+0.5_DP)*alat

       ! Loop over my atoms
       do iatom = my_first_at, my_last_at

          qi = real(elements(iatom)%ion_charge, kind=DP)

          ! Loop over all other atoms j<i
          do jatom = 1, iatom - 1

             ! Find charge and distance of other atom
             qj = real(elements(jatom)%ion_charge,kind=DP)
             rij = elements(iatom)%centre - elements(jatom)%centre

             ! Find force between atoms i and j
             call ewald_force_1d(rij,fij)
             fij = fij*qi*qj

             ! Add contribution of jatom to force on iatom
             ii_forces(:,iatom) = ii_forces(:,iatom) + fij

             ! Add contribution of iatom to force on jatom
             ii_forces(:,jatom) = ii_forces(:,jatom) - fij

          enddo

       enddo

    ! 2D periodicity
    case('SLAB')

       ! Set up shorthand lattice vectors and '12' plane unit normal
       a1 = pub_cell%a1
       a2 = pub_cell%a2
       a3 = a1.cross.a2
       a3 = (1.0d0/magnitude(a3))*a3
       a1p = (1.0d0/magnitude(a1))*a1
       a2p = a3.cross.a1p
       a2p = (1.0d0/magnitude(a2p))*a2p

       ! Set up reciprocal lattice vectors in '12' plane
       ! (not the same as pub_cell%b1,b2)
       b1=(2.0_DP*PI/( a1.dot.(a2.cross.a3) ) ) * (a2.cross.a3)
       b2=(2.0_DP*PI/( a1.dot.(a2.cross.a3) ) ) * (a3.cross.a1)

       ! Set up the inverse of the matrix comprising a1, a2, a3 to enable
       ! us to go back from coordinates relative to the periodicity to the
       ! simulation cell coordinates
       minv(1,1) = a1p%X; minv(2,1) = a1p%Y; minv(3,1) = a1p%Z
       minv(1,2) = a2p%X; minv(2,2) = a2p%Y; minv(3,2) = a2p%Z
       minv(1,3) = a3%X;  minv(2,3) = a3%Y;  minv(3,3) = a3%Z

       ! Set up constants
       Acell = magnitude(a1.cross.a2)
       sigma = sqrt(Acell)/2.4d0
       Rmax = sqrt(10.0d0*Acell/2.4d0)
       Gmax = sqrt(40.0d0*PI*PI/Acell)
       n1_real = int(Rmax/magnitude(a1))
       n2_real = int(Rmax/magnitude(a2))
       n1_recip = int(Gmax/magnitude(b1))
       n2_recip = int(Gmax/magnitude(b2))

       ! Loop over my atoms
       do iatom = my_first_at, my_last_at

          qi = real(elements(iatom)%ion_charge, kind=DP)

          ! Loop over all other atoms j<i
          do jatom = 1, iatom - 1

             ! Find charge and distance of other atom
             qj = real(elements(jatom)%ion_charge,kind=DP)
             rij = elements(iatom)%centre - elements(jatom)%centre

             ! Find force between atoms i and j
             call ewald_force_2d(rij,fij)
             fij(1:3) = fij(1:3)*qi*qj

             ! Add contribution of jatom to force on iatom
             ii_forces(1:3,iatom) = ii_forces(1:3,iatom) + fij(1:3)

             ! Add contribution of iatom to force on jatom
             ii_forces(1:3,jatom) = ii_forces(1:3,jatom) - fij(1:3)

          enddo

       enddo

    end select

    call comms_reduce('SUM', ii_forces)

  contains

    subroutine ewald_force_1d(rij,fij)
      !-----------------------------------------------------------------------!
      ! Evaluates the force on a particle at (0,0,0) due to a particle at rij !
      ! and its 1D periodic images.                                           !
      !-----------------------------------------------------------------------!

      implicit none

      type(POINT), intent(in) :: rij
      real(dp), intent(out) :: fij(3)

      real(dp), parameter :: E1 = -1.0_DP/24.0_DP, E2 = 7.0_DP/5760_DP

      integer :: m
      real(dp) :: rRm3
      real(dp) :: x, y, z, alpha
      real(dp) :: upz, umz, supz2a, sumz2a, upz2a, umz2a
      real(dp) :: f(3), fr(3), fem1(3), fem2(3)

      ! Set up constants determined from a1
      U = (real(mmax,DP)+0.5_DP)*alat

      ! Set up distance in a1 dir and distance squared perp to it
      z = rij.dot.a1
      x = rij.dot.a2
      y = rij.dot.a3
      alpha = x**2 + y**2

      ! 1/r term
      rRm3 = sqrt(z**2 + alpha)**3
      fr(1) = -x/rRm3
      fr(2) = -y/rRm3
      fr(3) = -z/rRm3

      ! 1/(r+Rm) terms from periodic interactions
      do m = 1, mmax
         rRm3 = sqrt((z + alat*real(m,dp))**2+alpha)**3
         fr(1) = fr(1) - x/rRm3
         fr(2) = fr(2) - y/rRm3
         fr(3) = fr(3) - (z + alat*real(m,dp))/rRm3
         rRm3 = sqrt((z - alat*real(m,dp))**2+alpha)**3
         fr(1) = fr(1) - x/rRm3
         fr(2) = fr(2) - y/rRm3
         fr(3) = fr(3) - (z - alat*real(m,dp))/rRm3
      enddo

      upz = U + z
      upz2a = upz*upz + alpha
      supz2a = sqrt(upz2a)
      umz = U - z
      umz2a = umz*umz + alpha
      sumz2a = sqrt(umz2a)

      ! First Euler Maclaurin terms
      fem1(1) = (x/(supz2a*(upz + supz2a)) &
           +  x/(sumz2a*(umz + sumz2a)))/alat
      fem1(2) = (y/(supz2a*(upz + supz2a)) &
           +  y/(sumz2a*(umz + sumz2a)))/alat
      fem1(3) = ((1.0d0 + upz/supz2a)/(upz + supz2a) &
           + (-1.0d0 - umz/sumz2a)/(umz + sumz2a))/alat

      ! Second Euler Maclaurin terms
      fem2(1) = alat*E1*(3.0d0*x*upz/supz2a**5 + 3.0d0*x*umz/sumz2a**5)
      fem2(2) = alat*E1*(3.0d0*y*upz/supz2a**5 + 3.0d0*x*umz/sumz2a**5)
      fem2(3) = -alat*E1*(1.0d0/supz2a**3 - 1.0d0/sumz2a**3 &
           - 3.0d0*upz**2/supz2a**5 + 3.0d0*umz**2/sumz2a**5)

      fem2(1) = -alat**3*E2*(105.d0*x*upz**3/supz2a**9 - 45.d0*x*upz/supz2a**7 &
           + 105.d0*x*umz**3/sumz2a**9 - 45.d0*x*umz/sumz2a**7)

      fem2(2) = -alat**3*E2*(105.d0*y*upz**3/supz2a**9 - 45.d0*y*upz/supz2a**7 &
           + 105.d0*y*umz**3/sumz2a**9 - 45.d0*y*umz/sumz2a**7)

      fem2(3) = -alat**3*E2*(9.0d0/supz2a**5 - 9.0d0/sumz2a**5 &
           - 90.0d0*upz**2/supz2a**7 + 90.0d0*umz**2/sumz2a**7 &
           + 105.0d0*upz**4/supz2a**9 - 105.0d0*umz**4/sumz2a**9)

      ! Above components were defined as d/dx(V), so negate them and add
      f = -fr - fem1 - fem2

      ! Transform back to cartesian coordinates of system
      fij = matmul(minv,f)

    end subroutine ewald_force_1d

    subroutine ewald_force_2d(rij,fij)
      !-----------------------------------------------------------------------!
      ! Evaluates the force on a particle at (0,0,0) due to a particle at rij !
      ! and its 2D periodic images.                                           !
      !-----------------------------------------------------------------------!

      use utils, only: utils_erf, utils_erfc

      implicit none

      ! Arguments
      type(POINT), intent(in) :: rij
      real(dp), intent(out) :: fij(3)

      ! Locals
      type(POINT) :: Rvec, Gvec, rpar ! IN NORMAL CARTESIANS
      real(dp) :: f(3), fr(3), fg(3), fzerf(3)
      real(dp) :: fac1, sinGr, cosGr
      real(dp) :: Rmag, Gmag, expfac1, erffac1, expfac2, erffac2
      real(dp) :: z
      integer :: i1, i2

      ! Set up perpendicular and parallel components of rij

      z = rij.dot.a3
      rpar = rij - z*a3

      fac1 = 2.0d0/(sqrt(PI)*sigma)

      ! Real Space Sum
      fr(1:3) = 0.0d0
      do i1 = -n1_real, n1_real
         do i2 = -n2_real, n2_real
            Rvec = real(i1,dp)*a1 + real(i2,dp)*a2
            Rmag = magnitude(Rvec)
            if (Rmag < Rmax) then
               Rvec = rij - Rvec
               Rmag = magnitude(Rvec)
               expfac1 = fac1*exp(-Rmag**2/sigma**2)/Rmag**2
               erffac1 = utils_erfc(Rmag/sigma)/Rmag**3
               fr(1) = fr(1) - Rvec%X*(expfac1 + erffac1)
               fr(2) = fr(2) - Rvec%Y*(expfac1 + erffac1)
               fr(3) = fr(3) - Rvec%Z*(expfac1 + erffac1)
            end if
         end do
      end do

      ! 2pi/A(z erf(z/sigma)+sigma/sqrt(pi)*exp(-z^2/sigma^2)) term
      fzerf(1) = 0.0d0
      fzerf(2) = 0.0d0
      fzerf(3) = -2.0d0*PI*utils_erf(z/sigma)/Acell

      ! Reciprocal Space Sum
      fac1 = 2.0d0*PI/Acell      ! includes k and -k
      m2_recip = 1
      fg(1:3) = 0.0d0
      do i1 = 0, n1_recip
         do i2 = m2_recip, n2_recip
            Gvec = real(i1,dp)*b1 + real(i2,dp)*b2
            Gmag = magnitude(Gvec)
            if (Gmag <= Gmax) then
               sinGr = sin(Gvec.dot.rpar)
               cosGr = cos(Gvec.dot.rpar)
               expfac1 = exp(+Gmag*z)
               expfac2 = exp(-Gmax*z)
               erffac1 = utils_erfc(sigma*Gmag*0.5d0 - z/sigma)
               erffac2 = utils_erfc(sigma*Gmag*0.5d0 + z/sigma)

               fg(1) = fg(1) - fac1*(Gvec.dot.a1p)/Gmag*( &
                    expfac2*erffac1 + expfac1*erffac2) * sinGr
               fg(2) = fg(2) - fac1*(Gvec.dot.a2p)/Gmag*( &
                    expfac2*erffac1 + expfac1*erffac2) * sinGr
               fg(3) = fg(3) + fac1/Gmag*( 2.0d0/(sqrt(Pi)*sigma)* &
                    (expfac2*exp(-(-z/sigma + Gmag*sigma*0.5d0)**2) &
                    -expfac1*exp(-( z/sigma + Gmag*sigma*0.5d0)**2)) &
                    - Gmag*expfac2*erffac1 + Gmag*expfac1*erffac2) * cosGr

            end if

         end do
         ! Now go from -n2_recip to n2_recip
         m2_recip = -n2_recip
      end do

      ! Translate this back to normal cartesian coordinates
      f(1:3) = matmul(minv,fzerf+fg)

      ! Above components were defined as d/dx(V), so negate them and add
      fij = -fr - f

    end subroutine ewald_force_2d

  end subroutine cutoff_coulomb_ii_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cutoff_coulomb_initial_guess(density_orig, elements, struct_fac)

    !=========================================================================!
    ! Wrapper for density_initial_guess_recip to switch to the padded cell,   !
    ! initialise the density on the padded grid and extract the relevant      !
    ! sections of the density to the standard cell.                           !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                                    !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use density, only : density_initial_guess_recip
    use ion, only : ELEMENT
    use utils, only : utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: density_orig(pub_fine_grid%ld1, &
         pub_fine_grid%ld2,pub_fine_grid%max_slabs12,pub_cell%num_spins)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    complex(kind=DP), intent(in) :: struct_fac(pub_padded_cell%num_pspecies, &
         pub_padded_grid%ld3, pub_padded_grid%ld2, &
         pub_padded_grid%max_slabs23)

    ! Locals
    real(kind=DP), allocatable :: density_pad(:,:,:,:) ! Padded density
    integer :: ierr          ! Error flag
    integer :: is            ! Spin counter

    ! jd: Takes care of padding between n1 and ld1 etc.
    density_orig = 0.0_DP

    ! Allocate padded density
    allocate(density_pad(pub_padded_grid%ld1, &
         pub_padded_grid%ld2,pub_padded_grid%max_slabs12, &
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('cutoff_coulomb_initial_guess', &
         'density_pad', ierr)

    ! Calculate the local pseudopotential on the padded fine grid
    call density_initial_guess_recip(density_pad, &
         elements, struct_fac, pub_padded_grid)

    ! Extract relevant sections of density from padded grid
    do is=1,pub_cell%num_spins
       call cc_padded_cell_extract(density_orig(:,:,:,is), &
            density_pad(:,:,:,is))
    end do

    ! Deallocate padded potential
    deallocate(density_pad,stat=ierr)
    call utils_dealloc_check('cutoff_coulomb_initial_guess',&
         'density_pad', ierr)

  end subroutine cutoff_coulomb_initial_guess


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cutoff_coulomb_struct_fac(struct_fac, elements)

    !========================================================================!
    ! Wrapper for pseudo_make_structure_factor to switch to the padded cell, !
    ! calculate the structure factor and switch back                         !
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                                   !
    !========================================================================!

    use ion, only : ELEMENT
    use pseudopotentials, only : pseudo_make_structure_factor

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: struct_fac(pub_padded_cell%num_pspecies, &
         pub_padded_grid%ld3, pub_padded_grid%ld2, pub_padded_grid%max_slabs23)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! Calculate structure factor for padded cell
    call pseudo_make_structure_factor(struct_fac, elements, pub_padded_grid)

  end subroutine cutoff_coulomb_struct_fac


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cutoff_coulomb_classical_sfac(struct_fac_classical)

    !========================================================================!
    ! Wrapper for pseudo_make_structure_factor to switch to the padded cell, !
    ! calculate the structure factor and switch back                         !
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                                   !
    !========================================================================!

    use ion, only : ELEMENT
    use classical_pot, only : classical_pot_struct_fac

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: struct_fac_classical( &
         pub_padded_grid%ld3, pub_padded_grid%ld2, &
         pub_padded_grid%max_slabs23)

    ! Calculate structure factor for padded cell
    call classical_pot_struct_fac(struct_fac_classical,pub_padded_grid)

  end subroutine cutoff_coulomb_classical_sfac


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cutoff_coulomb_localpseudo(localpseudo, struct_fac, &
       struct_fac_classical, elements)

    !=========================================================================!
    ! Wrapper for pseudopotentials_local_on_fine to switch to the padded      !
    ! cell, calculate the local pseudopotential on the padded grid and extract!
    ! the relevant sections to the original grid                              !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                                    !
    ! Modified by Jacek Dziedzic, 13/05/2010 to pass 'elements' to            !
    ! pseudopotentials_local_on_fine.                                         !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use ion, only: ELEMENT
    use paw, only: paw_tcore_hartree_on_grid
    use pseudopotentials, only : pseudopotentials_local_on_grid
    use rundat, only: pub_paw
    use utils, only : utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: localpseudo(pub_fine_grid%ld1, &
         pub_fine_grid%ld2,pub_fine_grid%max_slabs12)
    complex(kind=DP), intent(in) :: struct_fac(pub_padded_grid%ld3,&
         pub_padded_grid%ld2, pub_padded_grid%max_slabs23, &
         pub_padded_cell%num_pspecies)
    complex(kind=DP), intent(in) :: struct_fac_classical( &
         pub_padded_grid%ld3,pub_padded_grid%ld2,pub_padded_grid%max_slabs23)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! Locals
    ! Padded potential
    real(kind=DP), allocatable :: localpseudo_pad(:,:,:)
    integer :: ierr          ! Error flag

    ! Allocate padded potential
    allocate(localpseudo_pad(pub_padded_grid%ld1, &
         pub_padded_grid%ld2,pub_padded_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('cutoff_coulomb_localpseudo', &
         'localpseudo_pad', ierr)

    ! Calculate the local pseudopotential on the padded fine grid
    if (.not.pub_paw) then
       call pseudopotentials_local_on_grid(localpseudo_pad, &
            struct_fac,struct_fac_classical,pub_padded_grid,elements)
    else
       call paw_tcore_hartree_on_grid(localpseudo_pad, struct_fac, &
            struct_fac_classical,pub_padded_grid)
    end if

    ! Extract relevant sections of localpseudo from padded grid
    call cc_padded_cell_extract(localpseudo, localpseudo_pad)

    ! Deallocate padded potential
    deallocate(localpseudo_pad,stat=ierr)
    call utils_dealloc_check('cutoff_coulomb_localpseudo',&
         'localpseudo_pad', ierr)

  end subroutine cutoff_coulomb_localpseudo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cutoff_coulomb_core_density(core_density_orig,struct_fac)

    !=========================================================================!
    ! Wrapper for pseudopotentials_core_density to switch to the padded       !
    ! cell, calculate the core density on the padded grid and extract         !
    ! the relevant sections to the original grid                              !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, February 2009                                 !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use parallel_strategy, only: pub_elements_on_node
    use paw, only: paw_tcore_density
    use pseudopotentials, only : pseudopotentials_core_density
    use rundat, only: pub_paw
    use utils, only : utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: core_density_orig(pub_fine_grid%ld1, &
         pub_fine_grid%ld2,pub_fine_grid%max_slabs12)
    complex(kind=DP), intent(in) :: struct_fac(pub_padded_grid%ld3,&
         pub_padded_grid%ld2, pub_padded_grid%max_slabs23, &
         pub_padded_cell%num_pspecies)

    ! Locals
    ! Padded core density
    real(kind=DP), allocatable :: core_density_pad(:,:,:)
    integer :: ierr          ! Error flag

    ! Allocate padded potential
    allocate(core_density_pad(pub_padded_grid%ld1, &
         pub_padded_grid%ld2,pub_padded_grid%max_slabs12), &
         stat=ierr)
    call utils_alloc_check('cutoff_coulomb_core_density', &
         'core_density_pad', ierr)

    ! Calculate the core density on the padded fine grid
    if (.not.pub_paw) then
       call pseudopotentials_core_density(core_density_pad, &
            struct_fac,pub_padded_grid)
    else
       call paw_tcore_density(core_density_pad, &
            struct_fac,pub_padded_grid)
    end if

    ! Extract relevant sections of core density from padded grid
    call cc_padded_cell_extract(core_density_orig, core_density_pad)

    ! Deallocate padded core density
    deallocate(core_density_pad,stat=ierr)
    call utils_dealloc_check('cutoff_coulomb_core_density',&
         'core_density_pad', ierr)

  end subroutine cutoff_coulomb_core_density


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cutoff_coulomb_locps_forces(density_fine,elements,locps_forces)

    !=======================================================================!
    ! Wrapper for pseudo_local_calculate_forces to switch to the padded     !
    ! cell, then calculate the local pseudopotential forces in the padded   !
    ! cell and switch back                                                  !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                                  !
    !=======================================================================!

    use cell_grid, only: pub_fine_grid
    use ion, only : ELEMENT
    use paw, only: paw_tcore_hartree_calc_forces
    use pseudopotentials, only : pseudo_local_calculate_forces
    use rundat, only: pub_paw
    use utils, only : utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=dp), intent(inout) :: density_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2,pub_fine_grid%max_slabs12,pub_cell%num_spins)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    real(kind=dp), intent(out) :: locps_forces(1:3,pub_cell%nat)

    ! Locals
    real(kind=DP), allocatable :: density_fine_pad(:,:,:,:) ! Padded density
    integer :: ierr          ! Error flag
    integer :: is            ! Spin counter

    ! Allocate padded density
    allocate(density_fine_pad(pub_padded_grid%ld1,pub_padded_grid%ld2, &
         pub_padded_grid%max_slabs12,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('cutoff_coulomb_locps_forces','density_fine_pad', &
         ierr)

    ! Initialise padded density
    density_fine_pad = 0.0d0

    ! Deposit relevant sections of density to padded grid
    do is=1,pub_cell%num_spins
       call cc_padded_cell_deposit(density_fine_pad(:,:,:,is), &
            density_fine(:,:,:,is))
    end do

    ! Calculate local pseudopotential forces
    if (.not.pub_paw) then
       call pseudo_local_calculate_forces(density_fine_pad,pub_padded_grid, &
            elements,locps_forces)
    else
       call paw_tcore_hartree_calc_forces(density_fine_pad,pub_padded_grid, &
            elements,locps_forces)
    end if

    ! Deallocate padded density
    deallocate(density_fine_pad, stat=ierr)
    call utils_dealloc_check('cutoff_coulomb_locps_forces',&
         'density_fine_pad', ierr)

  end subroutine cutoff_coulomb_locps_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cutoff_coulomb_nlcc_forces(density_fine,core_density_fine, &
       elements,nlcc_forces,nhat_den_grad)

    !=======================================================================!
    ! Wrapper for pseudo_nlcc_calculate_forces to switch to the padded      !
    ! cell, calculate the nlcc forces in the padded cell and switch back    !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                                  !
    !=======================================================================!

    use cell_grid, only: pub_fine_grid
    use ion, only : ELEMENT
    use paw, only: paw_nlcc_calculate_forces
    use pseudopotentials, only : pseudo_nlcc_calculate_forces
    use rundat, only: pub_aug, pub_aug_den_dim
    use utils, only : utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=dp), intent(inout) :: density_fine(pub_fine_grid%ld1,&
        pub_fine_grid%ld2,pub_fine_grid%max_slabs12,pub_cell%num_spins)
    real(kind=DP), intent(in) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    real(kind=dp), intent(out) :: nlcc_forces(1:3,pub_cell%nat)
    real(kind=dp), optional, intent(in) :: nhat_den_grad(pub_fine_grid%ld1, &
         pub_fine_grid%ld2,pub_fine_grid%max_slabs12,0:pub_aug_den_dim, &
         pub_cell%num_spins)

    ! Locals
    real(kind=DP), allocatable :: density_fine_pad(:,:,:,:) ! Padded density
    real(kind=DP), allocatable :: core_density_fine_pad(:,:,:) ! Padded core
                                                               ! density
    real(kind=DP), allocatable :: nhat_den_grad_fine_pad(:,:,:,:,:) ! Padded 
                                                        ! compensation density
    integer :: ierr          ! Error flag
    integer :: is            ! Spin counter
    integer :: iaug          ! augmentation density components

    ! Allocate padded densities
    allocate(density_fine_pad(pub_padded_grid%ld1, &
         pub_padded_grid%ld2,pub_padded_grid%max_slabs12,pub_cell%num_spins), &
         stat=ierr)
    call utils_alloc_check('cutoff_coulomb_nlcc_forces','density_fine_pad', &
         ierr)
    allocate(core_density_fine_pad(pub_padded_grid%ld1, &
         pub_padded_grid%ld2,pub_padded_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('cutoff_coulomb_nlcc_forces', &
         'core_density_fine_pad',ierr)
    if (pub_aug) then
       allocate(nhat_den_grad_fine_pad(pub_padded_grid%ld1, &
            pub_padded_grid%ld2,pub_padded_grid%max_slabs12, &
            0:pub_aug_den_dim,pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('cutoff_coulomb_nlcc_forces', &
            'nhat_den_grad_fine_pad',ierr)
    end if

    ! Initialise padded density
    density_fine_pad = 0.0d0
    core_density_fine_pad = 0.0d0

    ! Deposit relevant sections of density (and compensation density) to
    ! padded grid
    do is=1,pub_cell%num_spins
       call cc_padded_cell_deposit(density_fine_pad(:,:,:,is), &
            density_fine(:,:,:,is))
       if (pub_aug) then
          do iaug=0,pub_aug_den_dim
             call cc_padded_cell_deposit(nhat_den_grad_fine_pad(:,:,:,iaug,is),&
                  nhat_den_grad(:,:,:,iaug,is))
          end do
       end if
    end do

    ! Deposit relevant sections of core density to padded grid
    call cc_padded_cell_deposit(core_density_fine_pad(:,:,:), &
         core_density_fine(:,:,:))

    ! Calculate nlcc pseudopotential forces
    if (.not.pub_aug) then
       call pseudo_nlcc_calculate_forces(density_fine_pad, &
            core_density_fine_pad,pub_padded_grid,elements,nlcc_forces)
    else
       call paw_nlcc_calculate_forces(density_fine_pad, &
            core_density_fine_pad,nhat_den_grad_fine_pad,pub_padded_grid, &
            elements,nlcc_forces)
    end if

    ! Deallocate padded densities
    if (pub_aug) then
       deallocate(nhat_den_grad_fine_pad,stat=ierr)
       call utils_dealloc_check('cutoff_coulomb_nlcc_forces',&
            'nhat_den_grad_fine_pad',ierr)
    end if
    deallocate(core_density_fine_pad,stat=ierr)
    call utils_dealloc_check('cutoff_coulomb_nlcc_forces',&
         'core_density_fine_pad',ierr)
    deallocate(density_fine_pad,stat=ierr)
    call utils_dealloc_check('cutoff_coulomb_nlcc_forces',&
         'density_fine_pad',ierr)

  end subroutine cutoff_coulomb_nlcc_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cutoff_coulomb_hartree(hartree_potential,density)

    !=======================================================================!
    ! Wrapper for hartree_on_grid to switch to the padded cell, then        !
    ! calculate the Hartree potential and extract the relevant sections for !
    ! original grid                                                         !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                                  !
    !=======================================================================!

    use cell_grid, only: pub_fine_grid
    use hartree, only : hartree_on_grid
    use utils, only : utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    ! Density in real space:
    real(kind=DP), intent(inout) :: density(pub_fine_grid%ld1, &
         pub_fine_grid%ld2,pub_fine_grid%max_slabs12,pub_cell%num_spins)
    ! Potential in real space:
    real(kind=DP), intent(out) :: hartree_potential(pub_fine_grid%ld1, &
         pub_fine_grid%ld2,pub_fine_grid%max_slabs12,pub_cell%num_spins)

    ! Locals
    real(kind=DP), allocatable :: density_pad(:,:,:,:) ! Padded density
    real(kind=DP), allocatable :: hartree_pad(:,:,:,:) ! Padded potential
    integer :: ierr          ! Error flag
    integer :: is            ! Spin counter

    ! Allocate padded density
    allocate(density_pad(pub_padded_grid%ld1, &
         pub_padded_grid%ld2,pub_padded_grid%max_slabs12,pub_cell%num_spins), &
         stat=ierr)
    call utils_alloc_check('cutoff_coulomb_hartree','density_pad',ierr)

    ! Deposit relevant sections of density to padded grid
    ! (density_pad is zeroed in cc_padded_cell_deposit)
    do is=1,pub_cell%num_spins
       call cc_padded_cell_deposit(density_pad(:,:,:,is), &
            density(:,:,:,is))
    end do

    ! Allocate padded Hartree potential
    allocate(hartree_pad(pub_padded_grid%ld1, &
         pub_padded_grid%ld2,pub_padded_grid%max_slabs12,pub_cell%num_spins), &
         stat=ierr)
    call utils_alloc_check('cutoff_coulomb_hartree','hartree_pad', &
         ierr)

    ! Calculate Cutoff Hartree potential
    call hartree_on_grid(hartree_pad,density_pad,pub_padded_grid)

    ! Extract relevant sections of hartree potential from padded grid
    do is=1,pub_cell%num_spins
       call cc_padded_cell_extract(hartree_potential(:,:,:,is), &
            hartree_pad(:,:,:,is))
    end do

    ! Deallocate padded Hartree potential
    deallocate(hartree_pad,stat=ierr)
    call utils_dealloc_check('cutoff_coulomb_hartree','hartree_pad', ierr)

    ! Deallocate padded density
    deallocate(density_pad, stat=ierr)
    call utils_dealloc_check('cutoff_coulomb_hartree','density_pad', ierr)

  end subroutine cutoff_coulomb_hartree


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cutoff_coulomb_init(elements)

    !=======================================================================!
    ! This subroutine initialises the padded grid information and calls     !
    ! parallel_strategy_determine for both the padded grid and the original !
    ! grid, and stores the results of both to allow the two grids to be     !
    ! exchanged quickly.                                                    !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                                  !
    !=======================================================================!

    use cell_grid, only: pub_fine_grid, pub_dbl_grid, cell_grid_distribute
    use comms, only : comms_abort, comms_reduce, pub_my_node_id, pub_on_root, &
         pub_total_num_nodes
    use constants, only: NORMAL
    use fourier, only : fourier_init_cell
    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(.cross.), &
         OPERATOR(-), OPERATOR(*), magnitude
    use ion, only : ELEMENT
    use parallel_strategy, only: pub_first_atom_on_node, pub_num_atoms_on_node
    use rundat, only : pub_coulomb_cutoff_type, pub_coulomb_radius, &
         pub_coulomb_length, pub_coulomb_cutoff_write_int, pub_output_detail, &
         pub_fine_grid_scale, pub_fine_is_dbl
    use simulation_cell, only: CELL_INFO
    use utils, only : utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! Copy in pub_cell parameters that weren't previously initialised
    pub_padded_cell%num_projectors = pub_cell%num_projectors
    pub_padded_cell%num_spins = pub_cell%num_spins
    pub_padded_cell%spin_fac = pub_cell%spin_fac

    ! Check padded_lattice_cart is bigger than lattice_cart
    call internal_check_padded_cell

    ! Print psinc grid info for padded cell
    !call cc_print_padded_cell

    ! Initialise parallel strategy information for padded cell
    call cell_grid_distribute(pub_padded_grid,pub_padded_cell, &
         pub_fine_grid_scale,'Padded Grid',.true.)

    if (.not.pub_fine_is_dbl) then

       ! ndmh: pub_fine_grid and pub_padded grid may have ended up with
       ! ndmh: different spacings
       if ((magnitude(pub_fine_grid%da1)/=magnitude(pub_padded_grid%da1)).or. &
            (magnitude(pub_fine_grid%da2)/=magnitude(pub_padded_grid%da2)).or. &
            (magnitude(pub_fine_grid%da3)/=magnitude(pub_padded_grid%da3))) then
          call utils_abort('Error in cutoff_coulomb_init: Padded and fine grids&
               & have different spacings')
       end if

    end if

    ! Initialise whole-cell FFT routines for padded cell
    call fourier_init_cell(pub_padded_grid)

    ! Check atoms are sufficiently far from cell boundaries that density does
    ! not spill out of home cell
    call internal_check_boundaries

    ! Check periodic images of simulation cell do not impinge on home cell
    call internal_check_images

    ! Set up transform of Coulomb interaction
    if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(a)',advance='no') &
         'Calculating transform of cutoff Coulomb interaction ...'

    call cc_interaction_on_recip(pub_padded_grid,pub_padded_grid%coulomb_recip,&
         pub_coulomb_cutoff_type,pub_coulomb_radius,pub_coulomb_length)

    if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(a)',advance='yes') ' done'

    ! Write out the real-space interaction if required
    if (pub_coulomb_cutoff_write_int) then
       call internal_write_int(pub_padded_grid%coulomb_recip,elements)
    end if

contains

    subroutine internal_write_int(coulomb_recip,elements)

      !=======================================================================!
      ! This subroutine writes out a .cube or .grd file of the interaction    !
      ! being used, by taking the Fourier transform of the terms stored in    !
      ! coulomb_grid.                                                         !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! coulomb_recip (input): interaction coefficients in recip space        !
      ! elements (input): list of all the atoms in the system                 !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2008                                  !
      !=======================================================================!

      use fourier, only : fourier_apply_cell_backward
      use visual, only : visual_scalarfield

      implicit none

      ! Arguments
      real(kind=DP), intent(in) :: coulomb_recip(pub_padded_grid%ld3,&
           pub_padded_grid%ld2,pub_padded_grid%max_slabs23)
      type(ELEMENT), intent(in) :: elements(pub_cell%nat)

      ! Local variables
      integer :: ierr              ! Error flag
      complex(kind=DP), allocatable :: v_recip_pad(:,:,:)
      real(kind=DP), allocatable :: v_real_pad(:,:,:)
      real(kind=DP), allocatable :: v_real(:,:,:)

      allocate(v_real_pad(pub_padded_grid%ld1, pub_padded_grid%ld2, &
           pub_padded_grid%max_slabs12),stat=ierr)
      call utils_alloc_check('internal_write_int (cutoff_coulomb_init)', &
           'v_real_pad',ierr)
      allocate(v_recip_pad(pub_padded_grid%ld3,pub_padded_grid%ld2, &
           pub_padded_grid%max_slabs23),stat=ierr)
      call utils_alloc_check('internal_write_int (cutoff_coulomb_init)', &
           'v_recip_pad',ierr)
      allocate(v_real(pub_fine_grid%ld1, pub_fine_grid%ld2, &
           pub_fine_grid%max_slabs12),stat=ierr)
      call utils_alloc_check('internal_write_int (cutoff_coulomb_init)', &
           'v_real',ierr)

      v_recip_pad(:,:,:) = cmplx(coulomb_recip(:,:,:), 0.0d0, kind=DP)

      ! FFT the local ionic potential from reciprocal to real space
      call fourier_apply_cell_backward(v_real_pad,v_recip_pad,pub_padded_grid)

      ! Write out the interaction on the normal cell
      !call visual_scalarfield(v_real_pad(:,:,:),pub_padded_grid, &
      !     'Coulomb interaction in real space:','_padcoulomb',elements)

      ! Extract relevant sections of v_real from padded grid
      call cc_padded_cell_extract(v_real,v_real_pad)

      ! Write out the interaction on the normal cell
      call visual_scalarfield(v_real(:,:,:),pub_fine_grid, &
           'Coulomb interaction in real space:','_coulomb',elements)

      ! Create periodic Coulomb interaction
      !v_recip_pad(:,:,:) = cmplx(1/recip_grid(4,:,:,:)**2, 0.0d0, kind=DP)
      !if (pub_on_root) v_recip_pad(1,1,1) = cmplx( 0.0d0, 0.0d0, kind=DP)

      ! FFT the local ionic potential from reciprocal to real space
      !call fourier_apply_cell_backward(v_real_pad,v_recip_pad)

      ! Extract relevant sections of v_real from padded grid
      !call cc_padded_cell_extract(v_real, v_real_pad)

      ! Write out the interaction on the normal cell
      !call visual_scalarfield(v_real(:,:,:), &
      !     'Coulomb interaction in real space:', elements, &
      !     '_percoulomb')

      ! Deallocate Arrays
      deallocate(v_real,stat=ierr)
      call utils_dealloc_check('internal_write_int (cutoff_coulomb_init)', &
           'v_real',ierr)
      deallocate(v_recip_pad,stat=ierr)
      call utils_dealloc_check('internal_write_int (cutoff_coulomb_init)', &
           'v_recip_pad',ierr)
      deallocate(v_real_pad,stat=ierr)
      call utils_dealloc_check('internal_write_int (cutoff_coulomb_init)', &
           'v_real_pad',ierr)

    end subroutine internal_write_int

    subroutine internal_check_padded_cell

      !=======================================================================!
      ! Checks that various quantities match between the padded and orginal   !
      ! cells in along all three axes                                         !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, December 2008                               !
      !=======================================================================!

      implicit none

      ! Locals
      real(kind=DP), parameter :: tol=0.00001_DP

      if (abs(pub_cell%d1-pub_padded_cell%d1)>tol) then
         if (pub_on_root) write(stdout,'(a/a)') 'Error in &
              &internal_check_padded_cell (cutoff_coulomb_init):', &
              &'d1 in padded cell does not match d1 in original cell'
         call comms_abort
      end if

      if (abs(pub_cell%d2-pub_padded_cell%d2)>tol) then
         if (pub_on_root) write(stdout,'(a/a)') 'Error in &
              &internal_check_padded_cell (cutoff_coulomb_init):', &
              &'d2 in padded cell does not match d2 in original cell'
         call comms_abort
      end if

      if (abs(pub_cell%d3-pub_padded_cell%d3)>tol) then
         if (pub_on_root) write(stdout,'(a/a)') 'Error in &
              &internal_check_padded_cell (cutoff_coulomb_init):', &
              &'d3 in padded cell does not match d3 in original cell'
         call comms_abort
      end if

      if (magnitude(pub_cell%a1) > magnitude(pub_padded_cell%a1)) then
         if (pub_on_root) write(stdout,'(a/a)') 'Error in &
              &internal_check_padded_cell (cutoff_coulomb_init):', &
              &'Padded cell a1 is smaller than original cell a1'
         call comms_abort
      end if

      if (magnitude(pub_cell%a2) > magnitude(pub_padded_cell%a2)) then
         if (pub_on_root) write(stdout,'(a/a)') 'Error in &
              &internal_check_padded_cell (cutoff_coulomb_init):', &
              &'Padded cell a2 is smaller than original cell a2'
         call comms_abort
      end if

      if (magnitude(pub_cell%a3) > magnitude(pub_padded_cell%a3)) then
         if (pub_on_root) write(stdout,'(a/a)') 'Error in &
              &internal_check_padded_cell (cutoff_coulomb_init):', &
              &'Padded cell a3 is smaller than original cell a3'
         call comms_abort
      end if

      if ( (pub_padded_cell%a1.dot.pub_cell%a1) < &
           magnitude(pub_cell%a1)*magnitude(pub_padded_cell%a1) ) then
         if (pub_on_root) write(stdout,'(a/a)') 'Error in &
              &internal_check_padded_cell (cutoff_coulomb_init):', &
              &'Padded cell a1 is not parallel to original cell a1'
         call comms_abort
      end if

      if ( (pub_padded_cell%a2.dot.pub_cell%a2) < &
           magnitude(pub_cell%a2)*magnitude(pub_padded_cell%a2) ) then
         if (pub_on_root) write(stdout,'(a/a)') 'Error in &
              &internal_check_padded_cell (cutoff_coulomb_init):', &
              &'Padded cell a2 is not parallel to original cell a2'
         call comms_abort
      end if

      if ( (pub_padded_cell%a3.dot.pub_cell%a3) < &
           magnitude(pub_cell%a3)*magnitude(pub_padded_cell%a3) ) then
         if (pub_on_root) write(stdout,'(a/a)') 'Error in &
              &internal_check_padded_cell (cutoff_coulomb_init):', &
              &'Padded cell a3 is not parallel to original cell a3'
         call comms_abort
      end if

    end subroutine internal_check_padded_cell

    subroutine internal_check_boundaries

      !=======================================================================!
      ! Cycles through the atoms in the elements array and checks that none   !
      ! of their NGWF radii are such that any of their electron density falls !
      ! outside of the home unit cell, since in cutoff coulomb such density   !
      ! will be in the wrong place in the padded cell.                        !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, April 2008                                  !
      !=======================================================================!

      implicit none

      ! Locals
      integer :: my_first_at     ! First atom on this node
      integer :: my_last_at      ! Last atom on this node
      integer :: nfail           ! Number of atoms spilling outside cell
      integer :: iatom           ! Atom index
      real(dp) :: rngwf          ! NGWF radius of current atom
      type(POINT) :: oo          ! Cell origin
      type(POINT) :: a1, a2, a3  ! Shorthand variables for lattice vectors
      type(POINT) :: rat         ! Position of current atom

      oo%X = 0.0d0
      oo%Y = 0.0d0
      oo%Z = 0.0d0
      a1 = pub_cell%a1
      a2 = pub_cell%a2
      a3 = pub_cell%a3

      my_first_at = pub_first_atom_on_node(pub_my_node_id)
      my_last_at = my_first_at + pub_num_atoms_on_node(pub_my_node_id) - 1

      nfail = 0

      if (pub_on_root.and.(pub_output_detail>=VERBOSE)) &
           write(stdout,'(a)',advance='no') &
           'Checking all ngwfs are fully within unit cell ...'

      ! Loop over my atoms
      do iatom = my_first_at, my_last_at

         rngwf = elements(iatom)%radius
         rat = elements(iatom)%centre

         select case (pub_coulomb_cutoff_type)
         case('NONE')
            return

         case('SPHERE','CYLINDER')

            if ((pointtoplanedist(rat,a1,a2,oo) < rngwf).or. &
                 (pointtoplanedist(rat,a2,a1,a3) < rngwf).or. &
                 (pointtoplanedist(rat,a2,a3,oo) < rngwf).or. &
                 (pointtoplanedist(rat,a3,a2,a1) < rngwf).or. &
                 (pointtoplanedist(rat,a3,a1,oo) < rngwf).or. &
                 (pointtoplanedist(rat,a1,a3,a2) < rngwf)) then
               nfail = nfail + 1
            end if

         case('WIRE')

            if ((pointtoplanedist(rat,a1,a2,oo) < rngwf).or. &
                 (pointtoplanedist(rat,a2,a1,a3) < rngwf).or. &
                 (pointtoplanedist(rat,a3,a1,oo) < rngwf).or. &
                 (pointtoplanedist(rat,a1,a3,a2) < rngwf)) then
               nfail = nfail + 1
            end if

         case('SLAB')

            if ((pointtoplanedist(rat,a1,a2,oo) < rngwf).or. &
                 (pointtoplanedist(rat,a2,a1,a3) < rngwf)) then
               nfail = nfail + 1
            end if

         end select

      enddo

      call comms_reduce('SUM',nfail)

      if (pub_on_root) then
         if (nfail > 0) then
            ! ndmh: Give error message (and what task we were doing if not
            ! ndmh: already mentioned).
            if (pub_output_detail<VERBOSE) write(stdout,'(a)',advance='no') &
                 'Checking all ngwfs are fully within unit cell ...'
            write(stdout,'(a)') ' failed'
            write(stdout,'(i5,a)') nfail, ' atom(s) have NGWFs overlapping at &
                 & least one cell boundary in'
            write(stdout,'(a)') ' a non-periodic direction'
            call comms_abort
         else
            if (pub_output_detail>=VERBOSE) write(stdout,'(a)') ' done'
         end if
      end if

    end subroutine internal_check_boundaries

    subroutine internal_check_images

      !=======================================================================!
      ! Checks the periodic image of the Coulomb interaction arising from     !
      ! each of the nearest neighbour simulation cells to see if it impinges  !
      ! upon the home unit cell.                                              !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2009                                  !
      !=======================================================================!

      use simulation_cell, only: cutoff_coulomb_tol

      implicit none

      type(POINT) :: image(3)        ! origins of periodic images of interaction
      type(POINT) :: a1,a2,a3        ! original cell lattice vectors
      character(4) :: string(3)      ! strings for lattice vectors for warnings
      character(80) :: warning(3)    ! warning string
      real(kind=DP) :: dmin          ! minimum allowable point-to-plane distance

      dmin = -1.0_DP ! qoh: Initialise to prevent compiler warning

      ! Find the vectors to each of the origins of the nearest neighbour cells
      ! in the +ve directions (no need to check -ve ones too by symmetry)
      image(1) = pub_padded_cell%a1
      image(2) = pub_padded_cell%a2
      image(3) = pub_padded_cell%a3

      string(1) = 'a1'
      string(2) = 'a2'
      string(3) = 'a3'

      a1 = pub_cell%a1
      a2 = pub_cell%a2
      a3 = pub_cell%a3

      warning(1) = 'WARNING in internal_check_images (cutoff_coulomb_init): '
      warning(2) = '     Periodic images of Coulomb interaction from &
           &padded_cell%'
      warning(3) = '     may impinge upon home unit cell. Check %block &
           &padded_lattice_cart.'

      if (pub_on_root) then

         ! Find cutoff along a1 axis
         select case (pub_coulomb_cutoff_type)
         case('CYLINDER')
            dmin = pub_coulomb_length + cutoff_coulomb_tol
         case('SPHERE')
            dmin = pub_coulomb_radius + cutoff_coulomb_tol
         case('SLAB','WIRE')
            dmin = 0.0_DP
         end select

         ! Check distance from padded cell a1 point to '23'-plane of orig cell
         if (abs(pointtoplanedist(image(1),a2,a3,a1))<dmin) then
            write(stdout,'(a)') trim(warning(1))
            write(stdout,'(a,a)') trim(warning(2)),trim(adjustl(string(1)))
            write(stdout,'(a)') trim(warning(3))
         end if

         ! Find cutoff along a2 axis
         select case (pub_coulomb_cutoff_type)
         case('CYLINDER','SPHERE','WIRE')
            dmin = pub_coulomb_radius + cutoff_coulomb_tol
         case('SLAB')
            dmin = 0.0_DP
         end select

         ! Check distance from padded cell a2 point to '13'-plane of orig cell
         if (abs(pointtoplanedist(image(2),a3,a1,a2))<dmin) then
            write(stdout,'(a)') trim(warning(1))
            write(stdout,'(a,a)') trim(warning(2)),trim(adjustl(string(2)))
            write(stdout,'(a)') trim(warning(3))
         end if

         ! Find cutoff along a3 axis
         select case (pub_coulomb_cutoff_type)
         case('CYLINDER','SPHERE','WIRE')
            dmin = pub_coulomb_radius + cutoff_coulomb_tol
         case('SLAB')
            dmin = pub_coulomb_length + cutoff_coulomb_tol
         end select

         ! Check distance from padded cell a3 point to '12'-plane of orig cell
         if (abs(pointtoplanedist(image(3),a1,a2,a3))<dmin) then
            write(stdout,'(a)') trim(warning(1))
            write(stdout,'(a,a)') trim(warning(2)),trim(adjustl(string(3)))
            write(stdout,'(a)') trim(warning(3))
         end if
      end if

    end subroutine internal_check_images

    real(dp) function pointtoplanedist(p,v1,v2,v0)

      !======================================================================!
      ! Finds the distance from point p to a plane defined by vectors v1, v2 !
      ! and origin v0                                                        !
      !----------------------------------------------------------------------!
      ! Written by Nicholas Hine, April 2008                                 !
      !======================================================================!

      implicit none

      ! Arguments
      type(POINT), intent(in) :: p       ! point to which to find distance
      type(POINT), intent(in) :: v1, v2  ! two vectors lying in the plane
      type(POINT), intent(in) :: v0      ! a point on the plane

      ! Locals
      type(POINT) :: n  ! plane normal

      n = v1.cross.v2
      n = (1.0d0/magnitude(n))*n
      pointtoplanedist = n.dot.(p - v0)

    end function pointtoplanedist

  end subroutine cutoff_coulomb_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cutoff_coulomb_exit

    !=======================================================================!
    ! This subroutine deallocates memory allocated by cutoff_coulomb_init   !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine, April 2008                                  !
    !=======================================================================!

    use cell_grid, only: cell_grid_exit, pub_fine_grid
    use utils, only : utils_dealloc_check

    implicit none

    ! Deallocate storage arrays for padded cell
    call cell_grid_exit(pub_padded_grid)

  end subroutine cutoff_coulomb_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !=========================================================================!
  ! INTERNAL SUBROUTINES                                                    !
  !=========================================================================!

  subroutine cc_padded_cell_extract(denpot,denpot_pad)

    !==========================================================!
    ! Takes a parallelised array of data for the padded cell   !
    ! and extracts from it the relevant slices of needed for a !
    ! parallelised array of data for the original cell         !
    !----------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                     !
    !==========================================================!

    use cell_grid, only: cell_grid_extract_box, pub_fine_grid
    use comms, only : pub_total_num_nodes, pub_my_node_id
    use utils, only : utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: denpot(pub_fine_grid%ld1, &
         pub_fine_grid%ld2,pub_fine_grid%max_slabs12)
    real(kind=DP), intent(in) :: denpot_pad(pub_padded_grid%ld1, &
         pub_padded_grid%ld2,pub_padded_grid%max_slabs12)

    ! Locals
    real(kind=DP), allocatable :: denpot_buffer(:,:,:)
    integer :: n_slabs12     ! Number of slabs of density_fine on node
    integer :: ierr          ! Error flag

    ! Allocate denpot_buffer
    allocate(denpot_buffer(pub_fine_grid%ld1, pub_fine_grid%ld2, &
         pub_padded_grid%max_group_slabs12), stat=ierr)
    call utils_alloc_check('cc_padded_cell_extract','denpot_buffer', &
         ierr)

    ! jd: This takes care of clearing the padding between n1 and ld1 etc.
    denpot = 0.0_DP
    denpot_buffer = 0.0_DP

    ! Extract relevant section of denpot from padded grid
    n_slabs12 = pub_fine_grid%last_slab12(pub_my_node_id) - &
         pub_fine_grid%first_slab12(pub_my_node_id) + 1
    call cell_grid_extract_box(denpot, denpot_buffer, &
         denpot_pad, pub_padded_grid, pub_fine_grid%n1, pub_fine_grid%n2, &
         n_slabs12, pub_fine_grid%ld1, pub_fine_grid%ld2, &
         1, 1, pub_fine_grid%first_slab12(pub_my_node_id), .true., .false.)

    ! Deallocate denpot_buffer
    deallocate(denpot_buffer, stat=ierr)
    call utils_dealloc_check('cc_padded_cell_extract','denpot_buffer', &
         ierr)

  end subroutine cc_padded_cell_extract


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cc_padded_cell_deposit(denpot_pad, denpot)

    !==========================================================!
    ! Takes a parallelised array of data for the original cell !
    ! and deposits it to a parallelised array of data for the  !
    ! padded cell                                              !
    !----------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                     !
    !==========================================================!

    use cell_grid, only: cell_grid_deposit_box, pub_fine_grid
    use comms, only : pub_total_num_nodes, pub_my_node_id
    use utils, only : utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: denpot_pad(pub_padded_grid%ld1, &
         pub_padded_grid%ld2,pub_padded_grid%max_slabs12,1)
    real(kind=DP), intent(in) :: denpot(pub_fine_grid%ld1, &
         pub_fine_grid%ld2,pub_fine_grid%max_slabs12,1)

    ! Locals
    real(kind=DP), allocatable :: denpot_buffer(:,:,:,:)
    integer :: n_slabs12     ! Number of slabs of density_fine on node
    integer :: ierr          ! Error flag

    ! Allocate denpot_buffer
    allocate(denpot_buffer(pub_fine_grid%ld1, pub_fine_grid%ld2, &
         pub_padded_grid%max_group_slabs12, 1), stat=ierr)
    call utils_alloc_check('cc_padded_cell_deposit','denpot_buffer', &
         ierr)

    ! Deposit relevant section of density to padded grid
    n_slabs12 = pub_fine_grid%last_slab12(pub_my_node_id) - &
         pub_fine_grid%first_slab12(pub_my_node_id) + 1


    ! jd: Need to clear these, as cell_grid_deposit_box relies on this
    denpot_pad = 0.0_DP
    denpot_buffer = 0.0_DP

    call cell_grid_deposit_box(denpot_pad, denpot, &
         denpot_buffer, pub_padded_grid, pub_fine_grid%n1, &
         pub_fine_grid%n2, n_slabs12, pub_fine_grid%ld1, &
         pub_fine_grid%ld2, 1, 1, &
         pub_fine_grid%first_slab12(pub_my_node_id), .true., .false.)

    ! Deallocate denpot_buffer
    deallocate(denpot_buffer, stat=ierr)
    call utils_dealloc_check('cc_padded_cell_deposit','denpot_buffer', &
         ierr)

  end subroutine cc_padded_cell_deposit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cc_print_padded_cell()

    use comms, only : pub_on_root

    implicit none

    ! print padded grids related info
    if (pub_on_root) then

       write(stdout,'(a)') ''
       write(stdout,'(a)') '============================== PSINC grid sizes ================================'
       write(stdout,'(3(a,i4))') '                          Padded cell: ',&
            pub_padded_cell%total_pt1 ,' x',pub_padded_cell%total_pt2 ,' x',pub_padded_cell%total_pt3

       write(stdout,'(3(a,i4))') '                                  PPD: ',&
            pub_padded_cell%n_pt1 ,' x',pub_padded_cell%n_pt2 ,' x',pub_padded_cell%n_pt3

       write(stdout,'(a)') '================================================================================'

    end if

  end subroutine cc_print_padded_cell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cc_interaction_on_recip(grid,coulomb_recip,cutoff_type,&
       coulomb_radius,coulomb_length)

    !==========================================================================!
    ! This subroutine calculates the Fourier components of the Coulomb         !
    ! potential for cutoff interactions on a sphere, cylinder, wire (infinite  !
    ! cylinder) or slab.                                                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine, March 2008                                     !
    !==========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only : comms_abort, comms_barrier, pub_on_root, pub_my_node_id, &
         pub_total_num_nodes
    use services, only : services_1d_interpolation
    use simulation_cell, only: CELL_INFO
    use timer, only : timer_clock
    use utils, only : utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(inout) :: grid
    real(kind=DP), intent(inout) :: coulomb_recip(grid%ld3,grid%ld2, &
         grid%max_slabs23)
    ! Coulomb cutoffs
    real(kind=DP), intent(in) :: coulomb_radius,coulomb_length
    ! Type of cutoff
    character(len=80) :: cutoff_type

    ! Locals
    integer :: i3,i2,islab23     ! Reciprocal grid loop counters
    integer :: my_first_slab23   ! First slab23 on this node
    integer :: n3,n2,nslab23     ! Reciprocal grid size shorthands
    integer :: nx,nr,nHx,nHxmax  ! Integration grid numbers of points
    integer :: ix,ir,iHx         ! Integration loop counters
    integer :: ierr              ! Allocation error flag

    real(kind=DP) :: Gvec(3), Glen
    real(kind=DP) :: Gx, Grho    ! Axial and Radial G-vectors for cylinder/wire
    real(kind=DP) :: Gz, Gpar    ! Out-of-Plane and In-Plane G-vectors for slab
    real(kind=DP) :: Rc, Lc      ! Coulomb cutoff lengths
    real(kind=DP) :: agrd, bgrd  ! Logarithmic grid variables
    real(kind=dp) :: rL

    ! Integration grid variables for Gx convolution
    real(kind=DP) :: Hxmax, dx, drho
    real(kind=DP),allocatable :: Hx(:),dHx(:),fHx(:,:)
    real(kind=DP),allocatable :: sinc(:,:),bessk0_aHxR(:),aHxR_bessk1_aHxR(:)
    real(kind=DP),allocatable :: bessj0_aGrR(:,:),aGrR_bessj1_aGrR(:,:)
    ! Integration grid variables for axial integral for cylinder
    real(kind=DP),allocatable :: x(:),fx(:)
    ! Integration grid variables for radial integral for cylinder/wire
    real(kind=DP),allocatable :: rho(:),fr(:),logr(:)

    call timer_clock('cc_interaction_on_recip',1)

    ! Initialise shorthand variables
    n3 = grid%n3
    n2 = grid%n2
    nslab23 = grid%num_slabs23
    my_first_slab23 = grid%first_slab23(pub_my_node_id)

    select case (cutoff_type)
    case('NONE')
       ! No cutoff - nothing to do
       return

    case('SPHERE')

       Rc = coulomb_radius

       ! Loop over reciprocal space grid on this node
       do islab23=1,nslab23       ! along b1
          do i2=1,n2/2+1          ! along b2
             do i3=1,n3/2+1       ! along b3
                call cell_grid_recip_pt(Gvec,islab23 + my_first_slab23 - 1, &
                     i2,i3,grid)
                Glen = sqrt(sum(Gvec(1:3)**2))

                coulomb_recip(i3,i2,islab23) = coulomb_recip(i3,i2,islab23) * &
                     (1.0_DP - cos(Glen * Rc))

             end do  ! b3

             do i3=n3/2+2,n3 ! Copy in negative G along b3
                coulomb_recip(i3,i2,islab23)= coulomb_recip(n3+2-i3,i2,islab23)
             end do ! b3

          end do ! b2

          do i2=n2/2+2,n2 ! Copy in negative G along b2
             coulomb_recip(1:n3,i2,islab23)= coulomb_recip(1:n3,n2+2-i2,islab23)
          end do ! b2

       end do ! b1

       ! Override for G=0 term
       if (pub_my_node_id==grid%node_slab23(1)) then
          coulomb_recip(1,1,1) = 0.5d0 * Rc**2
       end if

    case('CYLINDER')

       Rc = coulomb_radius
       Lc = coulomb_length

       ! Allocate arrays. ndmh_2011_09_08: larger arrays for all integrals
       ! to ensure more accurate convergence

       ! Prepare x arrays
       nx = 20001
       allocate(x(nx),fx(nx),stat=ierr)
       dx=Lc/real(nx-1,kind=DP)
       do ix = 1, nx
          x(ix) = real(ix-1,kind=DP)*dx
       enddo

       ! Prepare rho arrays
       nr = 20001
       allocate(rho(nr),fr(nr),logr(nr),stat=ierr)
       drho = Rc/real(nr-1,dp)
       rho(1) = 0.0_DP
       logr(1) = 0.0_DP
       do ir = 2, nr
          rho(ir) = real(ir-1,dp)*drho
          rL = sqrt(rho(ir)**2+Lc**2)
          logr(ir) = (log(rL+Lc)-log((rL-Lc)))*rho(ir)
       enddo

       ! Prepare Hx arrays
       nHxmax = 20001
       nHx = nHxmax
       allocate(Hx(nHx),dHx(nHx),fHx(nHx,2),stat=ierr)

       ! Prepare precalculated sinc and bessel function arrays
       allocate(sinc(nHx,2),stat=ierr)
       allocate(bessk0_aHxR(nHx),stat=ierr)
       allocate(aHxR_bessk1_aHxR(nHx),stat=ierr)
       allocate(bessj0_aGrR(n2/2+1,n3/2+1),stat=ierr)
       allocate(aGrR_bessj1_aGrR(n2/2+1,n3/2+1),stat=ierr)

       do i2=1,n2/2+1,1 ! positive G along b2
          do i3=1,n3/2+1,1 ! positive G along b3
             call cell_grid_recip_pt(Gvec,1,i2,i3,grid)
             Grho = sqrt(Gvec(2)**2 + Gvec(3)**2)
             bessj0_aGrR(i2,i3) = bessj0(abs(Grho)*Rc)
             aGrR_bessj1_aGrR(i2,i3) = abs(Grho)*Rc*bessj1(abs(Grho)*Rc)
          end do
       end do

       ! Loop over reciprocal space grid on this node
       do islab23=1,nslab23,1 !,-1 ! all G along b1
          call cell_grid_recip_pt(Gvec,islab23 + my_first_slab23 - 1,1,1, &
               grid)
          Gx = abs(Gvec(1))
          call cylinder_CC_general_newGx()
          if (pub_on_root) write(stdout,*) islab23

          do i2=1,n2/2+1,1 ! positive G along b2
             do i3=1,n3/2+1,1 ! positive G along b3

                call cell_grid_recip_pt(Gvec,islab23 + my_first_slab23 - 1, &
                     i2,i3,grid)
                Grho = sqrt(Gvec(2)**2 + Gvec(3)**2)

                if(Gx==0.0d0)then
                   if(Grho==0.0d0)then
                      coulomb_recip(i3,i2,islab23) = cylinder_CC_Gx0Gr0()
                   else
                      coulomb_recip(i3,i2,islab23) = cylinder_CC_Gx0()
                   end if
                else
                   if(Grho==0.0d0)then
                      coulomb_recip(i3,i2,islab23) = cylinder_CC_Gr0()
                   else
                      coulomb_recip(i3,i2,islab23) = cylinder_CC_general()
                   end if
                end if

             end do  ! b3

             do i3=n3/2+2,n3 ! Copy in negative G along b3
                coulomb_recip(i3,i2,islab23) = coulomb_recip(n3+2-i3,i2,islab23)
             end do ! b3

          end do ! b2

          do i2=n2/2+2,n2 ! Copy in negative G along b2
             coulomb_recip(1:n3,i2,islab23) = coulomb_recip(1:n3,n2+2-i2,islab23)
          end do ! b2

          i2 = n2/2; i3 = n3/2

       end do ! b1

       deallocate(aGrR_bessj1_aGrR,stat=ierr)
       deallocate(bessj0_aGrR,stat=ierr)
       deallocate(aHxR_bessk1_aHxR,stat=ierr)
       deallocate(bessk0_aHxR,stat=ierr)
       deallocate(sinc,stat=ierr)

       deallocate(Hx,dHx,fHx,stat=ierr)
       deallocate(rho,fr,logr,stat=ierr)
       deallocate(x,fx,stat=ierr)

    case('WIRE')

       Rc = coulomb_radius

       ! Prepare rho arrays
       nr = 10001
       allocate(rho(nr),stat=ierr)
       allocate(fr(nr),stat=ierr)
       drho = Rc/real(nr-1,dp)
       do ir = 1, nr
          rho(ir) = real(ir-1,dp)*drho
       enddo

       ! Loop over reciprocal space grid on this node
       do islab23=1,nslab23 ! all G along b1
          do i2=1,n2/2+1 ! positive G along b2
             do i3=1,n3/2+1 ! positive G along b3
             
                call cell_grid_recip_pt(Gvec,islab23 + my_first_slab23 - 1, &
                     i2,i3,grid)

                Gx = abs(Gvec(1))
                Grho = sqrt(Gvec(2)**2 + Gvec(3)**2)

                if(Gx==0.0d0)then
                   if(Grho==0.0d0)then
                      coulomb_recip(i3,i2,islab23) = wire_CC_Gx0Gr0()
                   else
                      coulomb_recip(i3,i2,islab23) = wire_CC_Gx0()
                   end if
                else
                   coulomb_recip(i3,i2,islab23) = wire_CC_general()
                end if

             end do  ! b3

             do i3=n3/2+2,n3 ! Copy in negative G along b3
                coulomb_recip(i3,i2,islab23) = coulomb_recip(n3+2-i3,i2,islab23)
             end do ! b3

          end do ! b2

          do i2=n2/2+2,n2 ! Copy in negative G along b2
             coulomb_recip(1:n3,i2,islab23) = coulomb_recip(1:n3,n2+2-i2,islab23)
          end do ! b2

       end do ! b1

       !call internal_print_terms()

       deallocate(fr,stat=ierr)
       deallocate(rho,stat=ierr)

    case('SLAB')

       Rc = coulomb_length

       ! Loop over reciprocal space grid on this node
       do islab23=1,nslab23 ! all G along b1

          do i2=1,n2/2+1 ! positive G along b2
             do i3=1,n3/2+1 ! positive G along b3
             
                call cell_grid_recip_pt(Gvec,islab23 + my_first_slab23 - 1, &
                     i2,i3,grid)

                Gz = abs(Gvec(3))
                Gpar = sqrt(Gvec(1)**2 + Gvec(2)**2)

                if(Gpar==0.0d0)then
                   if(Gz==0.0d0)then
                      coulomb_recip(i3,i2,islab23) = slab_CC_Gpar0Gz0()
                   else
                      coulomb_recip(i3,i2,islab23) = slab_CC_Gpar0()
                   end if
                else
                   coulomb_recip(i3,i2,islab23) = slab_CC_general()
                end if

             end do  ! b3

             do i3=n3/2+2,n3 ! Copy in negative G along b3
                coulomb_recip(i3,i2,islab23) = coulomb_recip(n3+2-i3,i2,islab23)
             end do ! b3

          end do ! b2

          do i2=n2/2+2,n2 ! Copy in negative G along b2
             coulomb_recip(1:n3,i2,islab23) = coulomb_recip(1:n3,n2+2-i2,islab23)
          end do ! b2

       end do ! b1

    end select

    call timer_clock('cc_interaction_on_recip',2)

  contains

    !=========================================================================!
    ! INTERNAL FUNCTIONS FOR SLAB CUTOFFS                                     !
    !=========================================================================!

    subroutine internal_print_terms

      !=======================================================================!
      ! Prints the terms of the Fourier transform of the interaction to stdout!
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2008                                  !
      !=======================================================================!

      integer :: node

      do node=0,pub_total_num_nodes-1

         if(node==pub_my_node_id)then
            do islab23=1,nslab23         ! all G along b1
               do i2=1,n2/2              ! positive G along b2
                  do i3=1,n3/2           ! positive G along b3
                     write(stdout,'(5f18.12)') coulomb_recip(i3,i2,islab23)
                  end do   ! b3
                  write(stdout,*)
               end do      ! b2
            end do         ! b1
         end if

         call comms_barrier()

      end do

    end subroutine internal_print_terms

    real(kind=dp) function slab_CC_Gpar0Gz0()
      !=======================================================================!
      ! This function returns the Fourier component of the Coulomb            !
      ! interaction cutoff on a slab for Gpar=0, Gz=0                         !
      !-----------------------------------------------------------------------!
      ! Takes parameter Rc from parent routine                                !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2008                                  !
      !=======================================================================!

      slab_CC_Gpar0Gz0 = -0.5d0*Rc**2

    end function slab_CC_Gpar0Gz0

    real(kind=dp) function slab_CC_Gpar0()
      !=======================================================================!
      ! This function returns the Fourier component of the Coulomb interaction!
      ! cutoff on a slab for Gpar=0, Gz>0                                     !
      !-----------------------------------------------------------------------!
      ! Takes parameters Gz, Rc from parent routine                           !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2008                                  !
      !=======================================================================!

      real(kind=dp) :: invGsq

      invGsq = 1.0d0/Gz**2

      slab_CC_Gpar0 = invGsq*(1 - (Gz*Rc)*sin(Gz*Rc) - cos(Gz*Rc))

    end function slab_CC_Gpar0

    real(kind=dp) function slab_CC_general()
      !=======================================================================!
      ! This function returns the Fourier component of the Coulomb interaction!
      ! cutoff on an infinite cylinder for Gx=0, Gr=0                         !
      !-----------------------------------------------------------------------!
      ! Takes parameters Gpar, Gz, Rc from parent routine                     !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2008                                  !
      !=======================================================================!

      real(kind=dp) :: invGsq

      invGsq = 1.0d0/(Gpar**2 + Gz**2)

      slab_CC_general = invGsq*(1 + &
           exp(-Gpar*Rc)*((Gz/Gpar)*sin(Gz*Rc) - cos(Gz*Rc)))

    end function slab_CC_general

    !=========================================================================!
    ! INTERNAL FUNCTIONS FOR WIRE CUTOFFS                                     !
    !=========================================================================!

    real(kind=dp) function wire_CC_Gx0Gr0()
      !=======================================================================!
      ! This function returns the Fourier component of the Coulomb interaction!
      ! cutoff on an infinite cylinder for Gx=0, Gr=0                         !
      !-----------------------------------------------------------------------!
      ! Uses Rc from parent routine                                           !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2008                                  !
      !=======================================================================!

      wire_CC_Gx0Gr0 = -0.25d0*Rc**2*(2.0d0*log(Rc) - 1.0d0)

    end function wire_CC_Gx0Gr0

    real(kind=dp) function wire_CC_Gx0()
      !=======================================================================!
      ! This function returns the Fourier component of the Coulomb interaction!
      ! cutoff on an infinite cylinder for Gx=0, Gr>0                         !
      !-----------------------------------------------------------------------!
      ! Takes parameters Grho, Rc from parent routine                         !
      ! Uses rho, drho and fr arrays from parent routine                      !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2008                                  !
      !=======================================================================!

      fr(1) = 0.0d0
      do ir = 2, nr
         fr(ir) = -log(rho(ir))*bessj0(Grho*rho(ir))*rho(ir)
      enddo

      wire_CC_Gx0 = internal_intreg(nr,drho,fr)

    end function wire_CC_Gx0

    real(kind=dp) function wire_CC_general()
      !=======================================================================!
      ! This function returns the Fourier component of the Coulomb interaction!
      ! cutoff on an infinite cylinder for Gx>0, Gr>=0                        !
      !-----------------------------------------------------------------------!
      ! Takes parameters Gx, Grho, Rc from parent routine                     !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2008                                  !
      !=======================================================================!

      real(kind=dp) :: invGsq,aGxR,aGrR

      invGsq = 1.0d0/(Grho**2 + Gx**2)
      aGxR = abs(Gx)*Rc
      aGrR = abs(Grho)*Rc
      wire_CC_general = invGsq*(1.0d0 + aGrR*bessj1(aGrR)*bessk0(aGxR) &
                                    & - aGxR*bessj0(aGrR)*bessk1(aGxR))

    end function wire_CC_general

    !=========================================================================!
    ! INTERNAL FUNCTIONS FOR CYLINDER CUTOFFS                                 !
    !=========================================================================!

    real(kind=dp) function cylinder_CC_Gx0Gr0()
      !=======================================================================!
      ! This function returns the Fourier component of the Coulomb interaction!
      ! cutoff on a cylinder for Gx=0, Gr=0                                   !
      !-----------------------------------------------------------------------!
      ! Uses Rc, Lc from parent routine                                       !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2008                                  !
      !=======================================================================!

      real(kind=dp) :: RLc

      RLc = sqrt(Rc**2+Lc**2)

      cylinder_CC_Gx0Gr0 = 0.5d0*(Lc*(RLc-Lc) + Rc**2*log((Lc+RLc)/Rc))

    end function cylinder_CC_Gx0Gr0

    real(kind=dp) function cylinder_CC_Gr0()
      !=======================================================================!
      ! This function returns the Fourier component of the Coulomb interaction!
      ! cutoff on a cylinder for Gx>0, Gr=0                                   !
      !-----------------------------------------------------------------------!
      ! Takes parameters Gx, Rc, Lc from parent routine                       !
      ! Uses x, dx and fx arrays from parent routine                          !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2008                                  !
      !=======================================================================!

      do ix = 1, nx
         fx(ix) = (sqrt(Rc**2 + x(ix)**2) - x(ix))*cos(Gx*x(ix))
      enddo

      cylinder_CC_Gr0 = internal_intreg(nx,dx,fx)

    end function cylinder_CC_Gr0

    real(kind=dp) function cylinder_CC_Gx0()
      !=======================================================================!
      ! This function returns the Fourier component of the Coulomb interaction!
      ! cutoff on a cylinder for Gx=0, Gr>0                                   !
      !-----------------------------------------------------------------------!
      ! Takes parameters Grho, Rc, Lc from parent routine                     !
      ! Uses rho, drho and fr arrays from parent routine                      !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2008                                  !
      !=======================================================================!

      fr(1) = 0.0d0
      do ir=2,nr
         fr(ir) = logr(ir)*bessj0(Grho*rho(ir))
      enddo

      cylinder_CC_Gx0 = 0.5d0*internal_intreg(nr,drho,fr)

    end function cylinder_CC_Gx0

    subroutine cylinder_CC_general_newGx()
      !=======================================================================!
      ! This subroutine prepares the arrays dGHx, sinc and aHxR which are     !
      ! used to speed up the evaluation of cylinder_CC_general                !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2009                                  !
      !=======================================================================!

      real(kind=dp) :: dGHx1,dGHx2,aHxR

      ! Prepare logarithmic grid for this value of Gx
      Hxmax = 100.0_DP/Rc
      bgrd = 14.0_DP/real(nHx,dp)
      agrd = real(Hxmax,dp) / (exp(bgrd*real(nHx-1,dp))-1.0_DP)
      do iHx = 1, nHx
         Hx(iHx) = agrd*(exp(bgrd*real(iHx,dp))-1.0_DP)
         dHx(iHx) = bgrd*(Hx(iHx)+agrd)
         aHxR = abs(Hx(iHx))*Rc
         bessk0_aHxR(iHx) = bessk0(aHxR)
         aHxR_bessk1_aHxR(iHx) = aHxR*bessk1(aHxR)
      enddo

      ! Prepare pre-calculated array for sinc function
      do iHx = 1, nHx
         dGHx1 = Gx - Hx(iHx)
         dGHx2 = Gx + Hx(iHx)
         sinc(iHx,1) = sin(dGHx1*Lc)/dGHx1
         sinc(iHx,2) = sin(dGHx2*Lc)/dGHx2
      enddo

    end subroutine cylinder_CC_general_newGx

    real(kind=dp) function cylinder_CC_general()
      !=======================================================================!
      ! This function returns the Fourier component of the Coulomb interaction!
      ! cutoff on a cylinder for Gx>0, Gr>0                                   !
      !-----------------------------------------------------------------------!
      ! Takes parameters Gx, Grho, Rc, Lc from parent routine                 !
      ! Uses Hx, dHx and fHx arrays from parent routine                       !
      !-----------------------------------------------------------------------!
      ! Written by Nicholas Hine, March 2008                                  !
      !=======================================================================!

      implicit none

      real(kind=DP) :: invGsq,bF,int1,int2,int3

      ! The part of the integral with no Bessel functions is analytic
      int1 = (1.0_DP + exp(-Grho*Lc)*( (Gx/Grho)*sin(Gx*Lc) - cos(Gx*Lc))) &
           & /(Gx**2 + Grho**2)

      ! The rest of the integral is split into negative and positive parts,
      ! calculated in int1 and int2 respectively. Because of the shared parts
      ! of both we evaluate them together

      do iHx = 1, nHx

         bF = (aGrR_bessj1_aGrR(i2,i3)*bessk0_aHxR(iHx) - &
              bessj0_aGrR(i2,i3)*aHxR_bessk1_aHxR(iHx))
         invGsq = 1.0_DP/(Grho**2 + Hx(iHx)**2)
         fHx(iHx,1) = sinc(iHx,1)*invGsq*bF
         fHx(iHx,2) = sinc(iHx,2)*invGsq*bF

      enddo

      int2 = internal_intlog(nHx,dHx,fHx(:,1))/PI
      int3 = internal_intlog(nHx,dHx,fHx(:,2))/PI

      cylinder_CC_general = int1 + int2 + int3

    end function cylinder_CC_general

    !=========================================================================!
    ! INTEGRATION ROUTINES                                                    !
    !=========================================================================!

    real(kind=dp) function internal_intlog(npts,dr,func)

      !=======================================================================!
      ! Simpson's rule integrator for function stored on a logarithmic grid.  !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! 1) npts  : input  : number of points in logarithmic grid (odd)        !
      ! 2) dr    : input  : the grid spacing                                  !
      ! 3) func  : input  : the function to be integrated                     !
      !=======================================================================!

      implicit none

      ! Arguments
      integer, intent(in) :: npts
      real(kind=dp), intent(in) :: dr(npts)
      real(kind=dp), intent(in) :: func(npts)

      ! Locals
      integer :: i
      real(kind=dp) :: r12,f1,f2,f3

      if(npts==2*(npts/2))then
         if(pub_on_root)then
            write(stdout,'(a)') 'Error: internal_intlog (cutoff_coulomb_mod): &
                 & npts should be odd'
            call comms_abort()
         end if
      endif

      internal_intlog = 0.0_dp
      r12 = 1.0_dp/12.0_dp
      f3 = func(1)*dr(1)*r12
      do i = 2,npts-1,2
         f1 = f3
         f2 = func(i)*dr(i)*r12
         f3 = func(i+1)*dr(i+1)*r12
         internal_intlog = internal_intlog + 5.0_dp*f1+8.0_dp*f2-1.0_dp*f3
         internal_intlog = internal_intlog - 1.0_dp*f1+8.0_dp*f2+5.0_dp*f3
      end do

      return

    end function internal_intlog

    real(kind=DP) function internal_intreg(npts,dr,func)

      !=======================================================================!
      ! Simpson's rule integrator for function stored on the regular grid.    !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! 1) npts  : input  : number of points in logarithmic grid (odd)        !
      ! 2) dr    : input  : the grid spacing                                  !
      ! 3) func  : input  : the function to be integrated                     !
      !=======================================================================!

      implicit none

      ! Arguments
      integer, intent(in) :: npts
      real(kind=dp), intent(in) :: dr
      real(kind=dp), intent(in) :: func(npts)

      ! Locals
      real(kind=dp) :: a(4)

      if(npts==2*(npts/2))then
         if(pub_on_root)then
            write(stdout,'(a)') 'Error: internal_intreg (cutoff_coulomb_mod): &
                 & npts should be odd'
            call comms_abort()
         end if
      endif

      a(1)=17.0_dp/48.0_dp
      a(2)=59.0_dp/48.0_dp
      a(3)=43.0_dp/48.0_dp
      a(4)=49.0_dp/48.0_dp

      internal_intreg = a(1)*(func(1)+func(npts))+a(2)*(func(2)+func(npts-1)) &
           + a(3)*(func(3)+func(npts-2))+a(4)*(func(4)+func(npts-3)) &
           + sum(func(5:npts-4))
      internal_intreg=internal_intreg*dr

      return

    end function internal_intreg


    !=========================================================================!
    ! BESSEL FUNCTION ROUTINES                                                !
    !=========================================================================!

      !============================================!
      ! Bessel J_0(x) function in double precision !
      !--------------------------------------------!
      ! Written by Nicholas Hine, 2008             !
      !============================================!

      real(kind=DP) function bessj0(x)

         implicit none

         ! Argument
         real(kind=DP), intent(in) :: x

         ! Locals
         real(kind=DP),parameter :: pi4 = PI*0.25_DP
         real(kind=DP) :: a(0:7), b(0:64), c(0:69), d(0:51)
         real(kind=DP) :: t,y,v,w,theta
         integer :: i,k

         data (a(i), i = 0, 7) /                                &
        &    -0.0000000000023655394d0, 0.0000000004708898680d0, &
        &    -0.0000000678167892231d0, 0.0000067816840038636d0, &
        &    -0.0004340277777716935d0, 0.0156249999999992397d0, &
        &    -0.2499999999999999638d0, 0.9999999999999999997d0 /
         data (b(i), i = 0, 12) /                               &
        &    0.0000000000626681117d0, -0.0000000022270614428d0, &
        &    0.0000000662981656302d0, -0.0000016268486502196d0, &
        &    0.0000321978384111685d0, -0.0005005237733315830d0, &
        &    0.0059060313537449816d0, -0.0505265323740109701d0, &
        &    0.2936432097610503985d0, -1.0482565081091638637d0, &
        &    1.9181123286040428113d0, -1.1319199475221700100d0, &
        &    -0.1965480952704682000d0 /
         data (b(i), i = 13, 25) /                              &
        &    0.0000000000457457332d0, -0.0000000015814772025d0, &
        &    0.0000000455487446311d0, -0.0000010735201286233d0, &
        &    0.0000202015179970014d0, -0.0002942392368203808d0, &
        &    0.0031801987726150648d0, -0.0239875209742846362d0, &
        &    0.1141447698973777641d0, -0.2766726722823530233d0, &
        &    0.1088620480970941648d0, 0.5136514645381999197d0,  &
        &    -0.2100594022073706033d0 /
         data (b(i), i = 26, 38) /                              &
        &    0.0000000000331366618d0, -0.0000000011119090229d0, &
        &    0.0000000308823040363d0, -0.0000006956602653104d0, &
        &    0.0000123499947481762d0, -0.0001662951945396180d0, &
        &    0.0016048663165678412d0, -0.0100785479932760966d0, &
        &    0.0328996815223415274d0, -0.0056168761733860688d0, &
        &    -0.2341096400274429386d0, 0.2551729256776404262d0, &
        &    0.2288438186148935667d0 /
         data (b(i), i = 39, 51) /                              &
        &    0.0000000000238007203d0, -0.0000000007731046439d0, &
        &    0.0000000206237001152d0, -0.0000004412291442285d0, &
        &    0.0000073107766249655d0, -0.0000891749801028666d0, &
        &    0.0007341654513841350d0, -0.0033303085445352071d0, &
        &    0.0015425853045205717d0, 0.0521100583113136379d0,  &
        &    -0.1334447768979217815d0, -0.1401330292364750968d0,&
        &    0.2685616168804818919d0 /
         data (b(i), i = 52, 64) /                              &
        &    0.0000000000169355950d0, -0.0000000005308092192d0, &
        &    0.0000000135323005576d0, -0.0000002726650587978d0, &
        &    0.0000041513240141760d0, -0.0000443353052220157d0, &
        &    0.0002815740758993879d0, -0.0004393235121629007d0, &
        &    -0.0067573531105799347d0, 0.0369141914660130814d0, &
        &    0.0081673361942996237d0, -0.2573381285898881860d0, &
        &    0.0459580257102978932d0 /
         data (c(i), i = 0, 13) /                                  &
        &    -0.00000000003009451757d0, -0.00000000014958003844d0, &
        &    0.00000000506854544776d0, 0.00000001863564222012d0,   &
        &    -0.00000060304249068078d0, -0.00000147686259937403d0, &
        &    0.00004714331342682714d0, 0.00006286305481740818d0,   &
        &    -0.00214137170594124344d0, -0.00089157336676889788d0, &
        &    0.04508258728666024989d0, -0.00490362805828762224d0,  &
        &    -0.27312196367405374426d0, 0.04193925184293450356d0 /
         data (c(i), i = 14, 27) /                                 &
        &    -0.00000000000712453560d0, -0.00000000041170814825d0, &
        &    0.00000000138012624364d0, 0.00000005704447670683d0,   &
        &    -0.00000019026363528842d0, -0.00000533925032409729d0, &
        &    0.00001736064885538091d0, 0.00030692619152608375d0,   &
        &    -0.00092598938200644367d0, -0.00917934265960017663d0, &
        &    0.02287952522866389076d0, 0.10545197546252853195d0,   &
        &    -0.16126443075752985095d0, -0.19392874768742235538d0 /
         data (c(i), i = 28, 41) /                                 &
        &    0.00000000002128344556d0, -0.00000000031053910272d0,  &
        &    -0.00000000334979293158d0, 0.00000004507232895050d0,  &
        &    0.00000036437959146427d0, -0.00000446421436266678d0,  &
        &    -0.00002523429344576552d0, 0.00027519882931758163d0,  &
        &    0.00097185076358599358d0, -0.00898326746345390692d0,  &
        &    -0.01665959196063987584d0, 0.11456933464891967814d0,  &
        &    0.07885001422733148815d0, -0.23664819446234712621d0 /
         data (c(i), i = 42, 55) /                                 &
        &    0.00000000003035295055d0, 0.00000000005486066835d0,   &
        &    -0.00000000501026824811d0, -0.00000000501246847860d0, &
        &    0.00000058012340163034d0, 0.00000016788922416169d0,   &
        &    -0.00004373270270147275d0, 0.00001183898532719802d0,  &
        &    0.00189863342862291449d0, -0.00113759249561636130d0,  &
        &    -0.03846797195329871681d0, 0.02389746880951420335d0,  &
        &    0.22837862066532347461d0, -0.06765394811166522844d0 /
         data (c(i), i = 56, 69) /                                 &
        &    0.00000000001279875977d0, 0.00000000035925958103d0,   &
        &    -0.00000000228037105967d0, -0.00000004852770517176d0, &
        &    0.00000028696428000189d0, 0.00000440131125178642d0,   &
        &    -0.00002366617753349105d0, -0.00024412456252884129d0, &
        &    0.00113028178539430542d0, 0.00708470513919789080d0,   &
        &    -0.02526914792327618386d0, -0.08006137953480093426d0, &
        &    0.16548380461475971846d0, 0.14688405470042110229d0 /
         data (d(i), i = 0, 12) /                                  &
        &    1.059601355592185731d-14, -2.71150591218550377d-13,   &
        &    8.6514809056201638d-12, -4.6264028554286627d-10,      &
        &    5.0815403835647104d-8, -1.76722552048141208d-5,       &
        &    0.16286750396763997378d0, 2.949651820598278873d-13,   &
        &    -8.818215611676125741d-12, 3.571119876162253451d-10,  &
        &    -2.631924120993717060d-8, 4.709502795656698909d-6,    &
        &    -5.208333333333283282d-3 /
         data (d(i), i = 13, 25) /                                 &
        &    7.18344107717531977d-15, -2.51623725588410308d-13,    &
        &    8.6017784918920604d-12, -4.6256876614290359d-10,      &
        &    5.0815343220437937d-8, -1.76722551764941970d-5,       &
        &    0.16286750396763433767d0, 2.2327570859680094777d-13,  &
        &    -8.464594853517051292d-12, 3.563766464349055183d-10,  &
        &    -2.631843986737892965d-8, 4.709502342288659410d-6,    &
        &    -5.2083333332278466225d-3 /
         data (d(i), i = 26, 38) /                                 &
        &    5.15413392842889366d-15, -2.27740238380640162d-13,    &
        &    8.4827767197609014d-12, -4.6224753682737618d-10,      &
        &    5.0814848128929134d-8, -1.76722547638767480d-5,       &
        &    0.16286750396748926663d0, 1.7316195320192170887d-13,  &
        &    -7.971122772293919646d-12, 3.544039469911895749d-10,  &
        &    -2.631443902081701081d-8, 4.709498228695400603d-6,    &
        &    -5.2083333315143653610d-3 /
         data (d(i), i = 39, 51) /                                 &
        &    3.84653681453798517d-15, -2.04464520778789011d-13,    &
        &    8.3089298605177838d-12, -4.6155016158412096d-10,      &
        &    5.0813263696466650d-8, -1.76722528311426167d-5,       &
        &    0.16286750396650065930d0, 1.3797879972460878797d-13,  &
        &    -7.448089381011684812d-12, 3.512733797106959780d-10,  &
        &    -2.630500895563592722d-8, 4.709483934775839193d-6,    &
        &    -5.2083333227940760113d-3 /
         w = abs(x)
         if (w .lt. 1.0d0) then
             t = w * w
             y = ((((((a(0) * t + a(1)) * t + &
        &        a(2)) * t + a(3)) * t + a(4)) * t + &
        &        a(5)) * t + a(6)) * t + a(7)
         else if (w .lt. 8.5d0) then
             t = w * w * 0.0625d0
             k = int(t)
             t = t - (k + 0.5d0)
             k = k * 13
             y = (((((((((((b(k) * t + b(k + 1)) * t + &
        &        b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + &
        &        b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + &
        &        b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + &
        &        b(k + 11)) * t + b(k + 12)
         else if (w .lt. 12.5d0) then
             k = int(w)
             t = w - (k + 0.5d0)
             k = 14 * (k - 8)
             y = ((((((((((((c(k) * t + c(k + 1)) * t + &
        &        c(k + 2)) * t + c(k + 3)) * t + c(k + 4)) * t + &
        &        c(k + 5)) * t + c(k + 6)) * t + c(k + 7)) * t + &
        &        c(k + 8)) * t + c(k + 9)) * t + c(k + 10)) * t + &
        &        c(k + 11)) * t + c(k + 12)) * t + c(k + 13)
         else
             v = 24.0d0 / w
             t = v * v
             k = 13 * (int(t))
             y = ((((((d(k) * t + d(k + 1)) * t + &
        &        d(k + 2)) * t + d(k + 3)) * t + d(k + 4)) * t + &
        &        d(k + 5)) * t + d(k + 6)) * sqrt(v)
             theta = (((((d(k + 7) * t + d(k + 8)) * t + &
        &        d(k + 9)) * t + d(k + 10)) * t + d(k + 11)) * t + &
        &        d(k + 12)) * v - pi4
             y = y * cos(w + theta)
         end if

         bessj0 = y

      end function bessj0


      !============================================!
      ! Bessel J_1(x) function in double precision !
      !--------------------------------------------!
      ! Written by Nicholas Hine, 2008             !
      !============================================!

      real(kind=DP) function bessj1(x)

         implicit none

         ! Argument
         real(kind=DP), intent(in) :: x

         ! Locals
         real(kind=DP),parameter :: pi4 = PI*0.25_DP
         real(kind=DP) :: a(0:7), b(0:64), c(0:69), d(0:51)
         real(kind=DP) :: t,y,v,w,theta
         integer :: i,k

         data (a(i), i = 0, 7) /                                  &
        &    -0.00000000000014810349d0, 0.00000000003363594618d0, &
        &    -0.00000000565140051697d0, 0.00000067816840144764d0, &
        &    -0.00005425347222188379d0, 0.00260416666666662438d0, &
        &    -0.06249999999999999799d0, 0.49999999999999999998d0 /
         data (b(i), i = 0, 12) /                                 &
        &    0.00000000000243721316d0, -0.00000000009400554763d0, &
        &    0.00000000306053389980d0, -0.00000008287270492518d0, &
        &    0.00000183020515991344d0, -0.00003219783841164382d0, &
        &    0.00043795830161515318d0, -0.00442952351530868999d0, &
        &    0.03157908273375945955d0, -0.14682160488052520107d0, &
        &    0.39309619054093640008d0, -0.47952808215101070280d0, &
        &    0.14148999344027125140d0 /
         data (b(i), i = 13, 25) /                                &
        &    0.00000000000182119257d0, -0.00000000006862117678d0, &
        &    0.00000000217327908360d0, -0.00000005693592917820d0, &
        &    0.00000120771046483277d0, -0.00002020151799736374d0, &
        &    0.00025745933218048448d0, -0.00238514907946126334d0, &
        &    0.01499220060892984289d0, -0.05707238494868888345d0, &
        &    0.10375225210588234727d0, -0.02721551202427354117d0, &
        &    -0.06420643306727498985d0 /
         data (b(i), i = 26, 38) /                                  &
        &    0.000000000001352611196d0, -0.000000000049706947875d0, &
        &    0.000000001527944986332d0, -0.000000038602878823401d0, &
        &    0.000000782618036237845d0, -0.000012349994748451100d0, &
        &    0.000145508295194426686d0, -0.001203649737425854162d0, &
        &    0.006299092495799005109d0, -0.016449840761170764763d0, &
        &    0.002106328565019748701d0, 0.058527410006860734650d0,  &
        &    -0.031896615709705053191d0 /
         data (b(i), i = 39, 51) /                                  &
        &    0.000000000000997982124d0, -0.000000000035702556073d0, &
        &    0.000000001062332772617d0, -0.000000025779624221725d0, &
        &    0.000000496382962683556d0, -0.000007310776625173004d0, &
        &    0.000078028107569541842d0, -0.000550624088538081113d0, &
        &    0.002081442840335570371d0, -0.000771292652260286633d0, &
        &    -0.019541271866742634199d0, 0.033361194224480445382d0, &
        &    0.017516628654559387164d0 /
         data (b(i), i = 52, 64) /                                  &
        &    0.000000000000731050661d0, -0.000000000025404499912d0, &
        &    0.000000000729360079088d0, -0.000000016915375004937d0, &
        &    0.000000306748319652546d0, -0.000004151324014331739d0, &
        &    0.000038793392054271497d0, -0.000211180556924525773d0, &
        &    0.000274577195102593786d0, 0.003378676555289966782d0,  &
        &    -0.013842821799754920148d0, -0.002041834048574905921d0,&
        &    0.032167266073736023299d0 /
         data (c(i), i = 0, 13) /                                   &
        &    -0.00000000001185964494d0, 0.00000000039110295657d0,   &
        &    0.00000000180385519493d0, -0.00000005575391345723d0,   &
        &    -0.00000018635897017174d0, 0.00000542738239401869d0,   &
        &    0.00001181490114244279d0, -0.00033000319398521070d0,   &
        &    -0.00037717832892725053d0, 0.01070685852970608288d0,   &
        &    0.00356629346707622489d0, -0.13524776185998074716d0,   &
        &    0.00980725611657523952d0, 0.27312196367405374425d0 /
         data (c(i), i = 14, 27) /                                  &
        &    -0.00000000003029591097d0, 0.00000000009259293559d0,   &
        &    0.00000000496321971223d0, -0.00000001518137078639d0,   &
        &    -0.00000057045127595547d0, 0.00000171237271302072d0,   &
        &    0.00004271400348035384d0, -0.00012152454198713258d0,   &
        &    -0.00184155714921474963d0, 0.00462994691003219055d0,   &
        &    0.03671737063840232452d0, -0.06863857568599167175d0,   &
        &    -0.21090395092505707655d0, 0.16126443075752985095d0 /
         data (c(i), i = 28, 41) /                                  &
        &    -0.00000000002197602080d0, -0.00000000027659100729d0,  &
        &    0.00000000374295124827d0, 0.00000003684765777023d0,    &
        &    -0.00000045072801091574d0, -0.00000327941630669276d0,  &
        &    0.00003571371554516300d0, 0.00017664005411843533d0,    &
        &    -0.00165119297594774104d0, -0.00485925381792986774d0,  &
        &    0.03593306985381680131d0, 0.04997877588191962563d0,    &
        &    -0.22913866929783936544d0, -0.07885001422733148814d0 /
         data (c(i), i = 42, 55) /                                  &
        &    0.00000000000516292316d0, -0.00000000039445956763d0,   &
        &    -0.00000000066220021263d0, 0.00000005511286218639d0,   &
        &    0.00000005012579400780d0, -0.00000522111059203425d0,   &
        &    -0.00000134311394455105d0, 0.00030612891890766805d0,   &
        &    -0.00007103391195326182d0, -0.00949316714311443491d0,  &
        &    0.00455036998246516948d0, 0.11540391585989614784d0,    &
        &    -0.04779493761902840455d0, -0.22837862066532347460d0 /
         data (c(i), i = 56, 69) /                                  &
        &    0.00000000002697817493d0, -0.00000000016633326949d0,   &
        &    -0.00000000433134860350d0, 0.00000002508404686362d0,   &
        &    0.00000048528284780984d0, -0.00000258267851112118d0,   &
        &    -0.00003521049080466759d0, 0.00016566324273339952d0,   &
        &    0.00146474737522491617d0, -0.00565140892697147306d0,   &
        &    -0.02833882055679300400d0, 0.07580744376982855057d0,   &
        &    0.16012275906960187978d0, -0.16548380461475971845d0 /
         data (d(i), i = 0, 12) /                                   &
        &    -1.272346002224188092d-14, 3.370464692346669075d-13,   &
        &    -1.144940314335484869d-11, 6.863141561083429745d-10,   &
        &    -9.491933932960924159d-8, 5.301676561445687562d-5,     &
        &    0.1628675039676399740d0, -3.652982212914147794d-13,    &
        &    1.151126750560028914d-11, -5.165585095674343486d-10,   &
        &    4.657991250060549892d-8, -1.186794704692706504d-5,     &
        &    1.562499999999994026d-2 /
         data (d(i), i = 13, 25) /                                  &
        &    -8.713069680903981555d-15, 3.140780373478474935d-13,   &
        &    -1.139089186076256597d-11, 6.862299023338785566d-10,   &
        &    -9.491926788274594674d-8, 5.301676558106268323d-5,     &
        &    0.1628675039676466220d0, -2.792555727162752006d-13,    &
        &    1.108650207651756807d-11, -5.156745588549830981d-10,   &
        &    4.657894859077370979d-8, -1.186794650130550256d-5,     &
        &    1.562499999987299901d-2 /
         data (d(i), i = 26, 38) /                                  &
        &    -6.304859171204770696d-15, 2.857249044208791652d-13,   &
        &    -1.124956921556753188d-11, 6.858482894906716661d-10,   &
        &    -9.491867953516898460d-8, 5.301676509057781574d-5,     &
        &    0.1628675039678191167d0, -2.185193490132496053d-13,    &
        &    1.048820673697426074d-11, -5.132819367467680132d-10,   &
        &    4.657409437372994220d-8, -1.186794150862988921d-5,     &
        &    1.562499999779270706d-2 /
         data (d(i), i = 39, 51) /                                  &
        &    -4.740417209792009850d-15, 2.578715253644144182d-13,   &
        &    -1.104148898414138857d-11, 6.850134201626289183d-10,   &
        &    -9.491678234174919640d-8, 5.301676277588728159d-5,     &
        &    0.1628675039690033136d0, -1.755122057493842290d-13,    &
        &    9.848723331445182397d-12, -5.094535425482245697d-10,   &
        &    4.656255982268609304d-8, -1.186792402114394891d-5,     &
        &    1.562499998712198636d-2 /
         w = abs(x)
         if (w .lt. 1) then
             t = w * w
             y = (((((((a(0) * t + a(1)) * t + &
        &        a(2)) * t + a(3)) * t + a(4)) * t + &
        &        a(5)) * t + a(6)) * t + a(7)) * w
         else if (w .lt. 8.5d0) then
             t = w * w * 0.0625d0
             k = int(t)
             t = t - (k + 0.5d0)
             k = k * 13
             y = ((((((((((((b(k) * t + b(k + 1)) * t + &
        &        b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + &
        &        b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + &
        &        b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + &
        &        b(k + 11)) * t + b(k + 12)) * w
         else if (w .lt. 12.5d0) then
             k = int(w)
             t = w - (k + 0.5d0)
             k = 14 * (k - 8)
             y = ((((((((((((c(k) * t + c(k + 1)) * t + &
        &        c(k + 2)) * t + c(k + 3)) * t + c(k + 4)) * t + &
        &        c(k + 5)) * t + c(k + 6)) * t + c(k + 7)) * t + &
        &        c(k + 8)) * t + c(k + 9)) * t + c(k + 10)) * t + &
        &        c(k + 11)) * t + c(k + 12)) * t + c(k + 13)
         else
             v = 24.0d0 / w
             t = v * v
             k = 13 * (int(t))
             y = ((((((d(k) * t + d(k + 1)) * t + &
        &        d(k + 2)) * t + d(k + 3)) * t + d(k + 4)) * t + &
        &        d(k + 5)) * t + d(k + 6)) * sqrt(v)
             theta = (((((d(k + 7) * t + d(k + 8)) * t + &
        &        d(k + 9)) * t + d(k + 10)) * t + d(k + 11)) * t + &
        &        d(k + 12)) * v - pi4
             y = y * sin(w + theta)
         end if

         if (x .lt. 0) y = -y

         bessj1 = y

      end function bessj1

      !============================================!
      ! Bessel K_0(x) function in double precision !
      !--------------------------------------------!
      ! Written by Nicholas Hine, 2008             !
      !============================================!

      real(kind=DP) function bessk0(x)

         implicit none

         ! Argument
         real(kind=DP), intent(in) :: x

         ! Locals
         real(kind=DP) :: a(0:15), b(0:111), c(0:134), d(0:39)
         real(kind=DP) :: t,y
         integer :: i,k

         data (a(i), i = 0, 15) /                                   &
        &    2.4307270476772195953d-12, 4.7091666363785304370d-10,  &
        &    6.7816861334344265568d-8, 6.7816840204737508252d-6,    &
        &    4.3402777777915334676d-4, 1.5624999999999872796d-2,    &
        &    2.5000000000000000448d-1, 9.9999999999999999997d-1,    &
        &    6.5878327432224993071d-12, 1.2083308769932888218d-9,   &
        &    1.6271062073716412046d-7, 1.4914719278555277887d-5,    &
        &    8.4603509071212245667d-4, 2.5248929932162333910d-2,    &
        &    2.7898287891460312491d-1, 1.1593151565841244874d-1 /
         data (b(i), i = 0, 13) /                                   &
        &    -4.6430702971053162197d-13, 1.0377936059563728230d-11, &
        &    -1.0298475936392057807d-10, 5.3632747492333959219d-10, &
        &    -2.1674628861036068105d-10, -2.3316071545820437669d-8, &
        &    2.2557819578691704059d-7, -9.2325694638587080009d-7,   &
        &    -3.3569097781613661759d-6, 8.7355061305812582974d-5,   &
        &    -6.8021202111645760475d-4, 2.7434654781323362319d-4,   &
        &    1.0031787169953909561d-1, 4.2102443824070833334d-1 /
         data (b(i), i = 14, 27) /                                  &
        &    4.1447451117883103686d-12, -3.4026589638576604315d-11, &
        &    9.3398790624638977468d-12, 1.5184181750799852630d-9,   &
        &    -1.1364911665083029464d-8, 2.0619457602095915719d-8,   &
        &    3.0431018037572243630d-7, -2.9749736264474555510d-6,   &
        &    8.0143661611467038568d-6, 8.0937525149549218398d-5,    &
        &    -1.0356346549612699886d-3, 2.8534806627578638795d-3,   &
        &    9.7369634474060441807d-2, 3.2175066577856452683d-1 /
         data (b(i), i = 28, 41) /                                  &
        &    1.1170882570740727520d-13, -8.2865909408297066068d-11, &
        &    9.4656678749191182763d-10, -3.5832019841847883380d-9,  &
        &    -9.5017955656904252761d-9, 1.5200595674883329093d-7,   &
        &    -3.8663262571356059980d-7, -3.3350340828235103499d-6,  &
        &    2.9359886663960844231d-5, -1.1266401822556801563d-5,   &
        &    -1.2113572742435576205d-3, 6.3158973673701376253d-3,   &
        &    8.8291790250128171341d-2, 2.2833982383240512262d-1 /
         data (b(i), i = 42, 55) /                                  &
        &    -3.2880638807053948433d-11, 4.3194884830465283512d-10, &
        &    -1.7455089683104033093d-9, -3.2437330799994764516d-9,  &
        &    4.7393655539139519778d-8, -1.1929265603456272466d-8,   &
        &    -1.3177845881013419388d-6, 3.3873375636197969526d-6,   &
        &    3.2729835880668256625d-5, -1.8367283883002494561d-4,   &
        &    -8.2830996454188084408d-4, 9.5512732229514251931d-3,   &
        &    7.2233832113719266702d-2, 1.4753187103603405298d-1 /
         data (b(i), i = 56, 69) /                                  &
        &    7.9998492614150860098d-11, -7.0257346702686139490d-10, &
        &    7.8898821627084586270d-10, 1.1294796399671507085d-8,   &
        &    -1.1360539648638059137d-8, -3.0346309115270564487d-7,  &
        &    3.2235585426189451721d-7, 8.3575612102298214948d-6,    &
        &    -8.5169628089198208211d-6, -2.5740175232173357342d-4,  &
        &    1.2462734014689152770d-4, 1.0683232869192203450d-2,    &
        &    5.1515690033637395779d-2, 8.5465862953544883657d-2 /
         data (b(i), i = 70, 83) /                                  &
        &    -8.6111506537356531608d-11, 5.1862926131024597823d-10, &
        &    7.5884324949371110022d-10, -6.4011975813006767417d-9,  &
        &    -4.1966181325111763156d-8, 9.1306285446881485314d-8,   &
        &    1.3573638315827954034d-6, 4.8683213252735694701d-7,    &
        &    -3.8805424608710197066d-5, -1.1838986468688980610d-4,  &
        &    9.2796213947750964945d-4, 8.9611057737319027776d-3,    &
        &    3.1464453915862785606d-2, 4.4267648087536630780d-2 /
         data (b(i), i = 84, 97) /                                  &
        &    4.4400123834164610288d-11, -1.1411233140911074336d-10, &
        &    -8.8200670702467059830d-10, -1.9686735373323381456d-9, &
        &    1.9921120728941773855d-8, 1.4543974418584834740d-7,    &
        &    1.8238418041265854754d-8, -4.5363700392899066037d-6,   &
        &    -2.1688068222527688542d-5, 4.5496062166687034700d-5,   &
        &    1.0435238076080528284d-3, 5.8374528996419979931d-3,    &
        &    1.6611210710425455850d-2, 2.0756008367065750538d-2 /
         data (b(i), i = 98, 111) /                                 &
        &    -6.5166519951106397214d-12, -5.8572182858788539580d-11,&
        &    1.5550375065815375404d-10, 1.9526509484993563229d-9,   &
        &    9.2637123346818426594d-9, -1.4136471501812055943d-8,   &
        &    -4.3024895710889717172d-7, -2.3235612243330592076d-6,  &
        &    4.0380616133862188804d-7, 9.2783767992909743602d-5,    &
        &    7.2964887597817095035d-4, 3.1316245282223273413d-3,    &
        &    7.8028233022066428316d-3, 9.0014807263791058095d-3 /
         data (c(i), i = 0, 14) /                                   &
        &    4.5161032649342790231d-11, -4.2774336988557091369d-11, &
        &    6.0998467173896677777d-10, 1.9845167242599996944d-9,   &
        &    1.3097678767280215271d-8, 7.4505822268382641286d-8,    &
        &    4.2893920879106814989d-7, 2.3900851955655303104d-6,    &
        &    1.2533473009382380357d-5, 5.9693359063879871983d-5,    &
        &    2.4775070661087304580d-4, 8.5106703131389516508d-4,    &
        &    2.2500105115665788755d-3, 4.0446134454521634600d-3,    &
        &    3.6910983340425942762d-3 /
         data (c(i), i = 15, 29) /                                  &
        &    3.5732826433251464989d-12, -3.2906649482312266258d-12, &
        &    7.0873811190464760555d-11, 2.9551320580484177120d-10,  &
        &    2.2776940472505079894d-9, 1.5175463612815010036d-8,    &
        &    9.9462487812170164133d-8, 6.1448757797853901100d-7,    &
        &    3.4869531882907360750d-6, 1.7615836644757657443d-5,    &
        &    7.6373536037879531886d-5, 2.7098571871205999668d-4,    &
        &    7.3399047381788927036d-4, 1.3439197177355085297d-3,    &
        &    1.2439943280131230863d-3 /
         data (c(i), i = 30, 44) /                                  &
        &    3.6343547836242523646d-13, 9.7997961751276137602d-14,  &
        &    1.0184692699811569047d-11, 6.1495184828957652064d-11,  &
        &    5.0238328349302602543d-10, 3.7498626376004337661d-9,   &
        &    2.6689445483857236307d-8, 1.7591899737346368084d-7,    &
        &    1.0486448307010701679d-6, 5.4986458466257148573d-6,    &
        &    2.4521456351751345323d-5, 8.8900942259143832228d-5,    &
        &    2.4483947714068300190d-4, 4.5418248688489693045d-4,    &
        &    4.2479574186923180694d-4 /
         data (c(i), i = 45, 59) /                                  &
        &    5.2460389348163395857d-14, 7.4802063026503503540d-14,  &
        &    2.0012201610651998417d-12, 1.4887306044735163359d-11,  &
        &    1.2946705414232940350d-10, 1.0391628915892803144d-9,   &
        &    7.8091180499677328456d-9, 5.3694223626907660084d-8,    &
        &    3.3063914804658509029d-7, 1.7776972424421486506d-6,    &
        &    8.0833148098458320202d-6, 2.9755556304448817780d-5,    &
        &    8.2945928349220642178d-5, 1.5536921180500112883d-4,    &
        &    1.4647070522281538711d-4 /
         data (c(i), i = 60, 74) /                                  &
        &    9.7531436733955514559d-15, 2.4084291220447154982d-14,  &
        &    4.7654956400897494468d-13, 4.0200949504810597783d-12,  &
        &    3.6726577109162191533d-11, 3.0939005665422637601d-10,  &
        &    2.4122848979784500179d-9, 1.7071884462645525505d-8,    &
        &    1.0752238955654933405d-7, 5.8844190041189462347d-7,    &
        &    2.7136083303224014597d-6, 1.0102477728604441135d-5,    &
        &    2.8420490721532571809d-5, 5.3637016379451944413d-5,    &
        &    5.0881312956459247572d-5 /
         data (c(i), i = 75, 89) /                                  &
        &    2.1732049868189377260d-15, 7.2720052142815590531d-15,  &
        &    1.2803083795536820100d-13, 1.1696825543787717167d-12,  &
        &    1.1083298191597132094d-11, 9.6536661252658773139d-11,  &
        &    7.7242553835198536397d-10, 5.5798366267110575620d-9,   &
        &    3.5721345296543414370d-8, 1.9806931547193682466d-7,    &
        &    9.2312964655319555313d-7, 3.4666258590861079959d-6,    &
        &    9.8224698307751177077d-6, 1.8648773453825584428d-5,    &
        &    1.7780062316167651812d-5 /
         data (c(i), i = 90, 104) /                                 &
        &    5.5012463763851934112d-16, 2.2254763392767319419d-15,  &
        &    3.7187669817701214965d-14, 3.5819585377733489628d-13,  &
        &    3.4866061263191556694d-12, 3.1101633450629652910d-11,  &
        &    2.5358235662235617663d-10, 1.8597629779492599046d-9,   &
        &    1.2052654739462999992d-8, 6.7501417351172136833d-8,    &
        &    3.1720052198654584574d-7, 1.1993651363602981832d-6,    &
        &    3.4179130317623363474d-6, 6.5208606745808860158d-6,    &
        &    6.2430205476536771454d-6 /
         data (c(i), i = 105, 119) /                                &
        &    1.5225407517829491689d-16, 6.9834820025664405161d-16,  &
        &    1.1380182837138781431d-14, 1.1369488761077196511d-13,  &
        &    1.1291168681618466716d-12, 1.0250757630526871007d-11,  &
        &    8.4765287317253141514d-11, 6.2886627779402596211d-10,  &
        &    4.1142865598366029316d-9, 2.3223773435632014408d-8,    &
        &    1.0985095234166396934d-7, 4.1766260951820336228d-7,    &
        &    1.1958609263543792991d-6, 2.2907574647671878055d-6,    &
        &    2.2008253973114914005d-6 /
         data (c(i), i = 120, 134) /                                &
        &    4.4863058691420695911d-17, 2.2437356594371819978d-16,  &
        &    3.6107964803015652759d-15, 3.7031193629853392081d-14,  &
        &    3.7341552790439784371d-13, 3.4355950129497564468d-12,  &
        &    2.8719942600171304499d-11, 2.1499646844509516453d-10,  &
        &    1.4171810843455227171d-9, 8.0501118772875784153d-9,    &
        &    3.8281889106330295876d-8, 1.4621673458431979989d-7,    &
        &    4.2029868696411098586d-7, 8.0785884122023473025d-7,    &
        &    7.7845438614204963209d-7 /
         data (d(i), i = 0, 7) /                                    &
        &    -7.9737703860537066166d-14, 1.9543834380466766627d-12, &
        &    -4.7230794431646733538d-11, 1.4001773785771252004d-9,  &
        &    -5.4864553020583098585d-8, 3.1601984250143742772d-6,   &
        &    -3.3708783204090252161d-4, 1.6180215937964160437d-1 /
         data (d(i), i = 8, 15) /                                   &
        &    -5.2593898374798632343d-14, 1.7725913926973236457d-12, &
        &    -4.6672234858122387294d-11, 1.3991653503828889207d-9,  &
        &    -5.4863400156413929639d-8, 3.1601976099900075541d-6,   &
        &    -3.3708783171335864627d-4, 1.6180215937958433760d-1 /
         data (d(i), i = 16, 23) /                                  &
        &    -3.6135496189875398132d-14, 1.5466239429618130284d-12, &
        &    -4.5320259146602122624d-11, 1.3945974109459385552d-9,  &
        &    -5.4853994841172088787d-8, 3.1601858228022739196d-6,   &
        &    -3.3708782339998302320d-4, 1.6180215937704286491d-1 /
         data (d(i), i = 24, 31) /                                  &
        &    -2.5640663123518180635d-14, 1.3288079339404032671d-12, &
        &    -4.3368537955908371563d-11, 1.3848103653102203186d-9,  &
        &    -5.4824335664256344123d-8, 3.1601315173126153586d-6,   &
        &    -3.3708776779035695640d-4, 1.6180215935248373474d-1 /
         data (d(i), i = 32, 39) /                                  &
        &    -1.8678321325292127767d-14, 1.1354310934105733311d-12, &
        &    -4.1057197297998608931d-11, 1.3693990961296350970d-9,  &
        &    -5.4762428935047089835d-8, 3.1599817092775027963d-6,   &
        &    -3.3708756559715893599d-4, 1.6180215923508144240d-1 /
         if (x .lt. 0.86d0) then
             t = x * x
             y = ((((((a(0) * t + a(1)) * t + &
        &        a(2)) * t + a(3)) * t + a(4)) * t + &
        &        a(5)) * t + a(6)) * t + a(7)
             y = ((((((a(8) * t + a(9)) * t + &
        &        a(10)) * t + a(11)) * t + a(12)) * t + &
        &        a(13)) * t + a(14)) * t + a(15) - y * log(x)
         else if (x .lt. 4.15d0) then
             t = x - 5 / x
             k = int(t + 5)
             t = (k - 4) - t
             k = k * 14
             y = ((((((((((((b(k) * t + b(k + 1)) * t + &
        &        b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + &
        &        b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + &
        &        b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + &
        &        b(k + 11)) * t + b(k + 12)) * t + b(k + 13)
         else if (x .lt. 12.5d0) then
             k = int(x)
             t = (k + 1) - x
             k = 15 * (k - 4)
             y = (((((((((((((c(k) * t + c(k + 1)) * t + &
        &        c(k + 2)) * t + c(k + 3)) * t + c(k + 4)) * t + &
        &        c(k + 5)) * t + c(k + 6)) * t + c(k + 7)) * t + &
        &        c(k + 8)) * t + c(k + 9)) * t + c(k + 10)) * t + &
        &        c(k + 11)) * t + c(k + 12)) * t + c(k + 13)) * t + &
        &        c(k + 14)
         else
             t = 60 / x
             k = 8 * (int(t))
             y = (((((((d(k) * t + d(k + 1)) * t + &
        &        d(k + 2)) * t + d(k + 3)) * t + d(k + 4)) * t + &
        &        d(k + 5)) * t + d(k + 6)) * t + d(k + 7)) * &
        &        sqrt(t) * exp(-x)
         end if

         bessk0 = y

      end function bessk0

      !============================================!
      ! Bessel K_1(x) function in double precision !
      !--------------------------------------------!
      ! Written by Nicholas Hine, 2008             !
      !============================================!

      real(kind=DP) function bessk1(x)

         implicit none

         ! Argument
         real(kind=DP), intent(in) :: x

         ! Locals
         real(kind=DP) :: a(0:15), b(0:119), c(0:134), d(0:39)
         real(kind=DP) :: t,y,v
         integer :: i,k

         data (a(i), i = 0, 15) /                                   &
        &    1.5151605362537935201d-13, 3.3637909513536510350d-11,  &
        &    5.6514041131016827202d-9, 6.7816840255069534052d-7,    &
        &    5.4253472222259226487d-5, 2.6041666666666637057d-3,    &
        &    6.2500000000000000090d-2, 5.0000000000000000000d-1,    &
        &    -8.9790303384748696588d-11, -1.4029047449249185771d-8, &
        &    -1.5592893752540998113d-6, -1.1253607018469017569d-4,  &
        &    -4.6421827665011579173d-3, -8.5370719728648409609d-2,  &
        &    -3.0796575782920629660d-1, 1.0000000000000000004d0 /
         data (b(i), i = 0, 14) /                                   &
        &    -9.4055461896630579928d-12, 3.1307934665844664773d-11, &
        &    4.2005295001519243251d-10, -4.1636196779679820012d-9,  &
        &    1.4483026181700966164d-8, 1.1661000205428816914d-8,    &
        &    -3.5023996724943046209d-7, 1.4404279316339005012d-6,   &
        &    5.3581564157158242080d-7, -3.5249754038612334639d-5,   &
        &    1.7150324075631641453d-4, -4.1276362611239191024d-5,   &
        &    -4.6943110979636602591d-3, 3.5085369853392357659d-2,   &
        &    2.0063574339907819159d-1 /
         data (b(i), i = 15, 29) /                                  &
        &    3.3998989888944034586d-11, 7.1558979072937373055d-11,  &
        &    -2.9226856932927698732d-9, 1.4591620256525610213d-8,   &
        &    -6.6141635609854161666d-9, -1.9991101838984472332d-7,  &
        &    5.9185836628873530572d-7, 1.9562880347358085687d-6,    &
        &    -1.5814366450418102764d-5, 7.6791682910944612028d-6,   &
        &    2.8354678948323983936d-4, -1.0217932669418690641d-3,   &
        &    -3.2205661865210048433d-3, 4.3497494842354644077d-2,   &
        &    1.6110284302315089935d-1 /
         data (b(i), i = 30, 44) /                                  &
        &    -2.0933987679665541827d-10, 7.9503322090520447316d-10, &
        &    3.8000150948242575774d-9, -2.3076136195585571309d-8,   &
        &    -2.3902880302550799653d-8, 3.1914500937804377478d-7,   &
        &    3.2639909831082417694d-7, -5.3166994792995439449d-6,   &
        &    -3.1109524694269240094d-6, 9.2575906966353273247d-5,   &
        &    7.5063709094147644746d-7, -1.7416491592625765379d-3,   &
        &    1.2138560335171676007d-3, 4.5879687144659643175d-2,    &
        &    1.1566544716132846709d-1 /
         data (b(i), i = 45, 59) /                                  &
        &    3.1582384905164908749d-10, -1.9959561818098999516d-9,  &
        &    8.6959328920030927557d-10, 1.1642778282445577109d-8,   &
        &    4.3552264337818440471d-8, -1.5057982160481803238d-7,   &
        &    -1.0101729117980989857d-6, 7.7002002510177612013d-7,   &
        &    1.9580574235590194233d-5, 1.9358461980242834361d-5,    &
        &    -3.3932339942485532728d-4, -9.3416673584325090073d-4,  &
        &    5.5800080455912847227d-3, 3.8668683062477179235d-2,    &
        &    7.2651643500517000658d-2 /
         data (b(i), i = 60, 74) /                                  &
        &    -1.1554749629758510059d-10, 8.2270678758893273006d-10, &
        &    -5.0211156951551538591d-10, -1.4929179050554858361d-9, &
        &    -2.7107940791526366702d-8, -4.2204764086705349384d-8,  &
        &    3.7253098167927628867d-7, 2.4374697215363361156d-6,    &
        &    1.4141942006909768370d-6, -4.8766389019473918231d-5,   &
        &    -2.1681387247526720978d-4, 2.9325729929653405236d-4,   &
        &    6.4087534504827239815d-3, 2.6054628289709454356d-2,    &
        &    4.0156431128194184336d-2 /
         data (b(i), i = 75, 89) /                                  &
        &    2.5506555170746221691d-11, -1.3521164018407978152d-10, &
        &    -8.3281235274106699399d-11, -9.7764849575562351891d-10,&
        &    3.4661828409940354542d-9, 3.9760633711791357544d-8,    &
        &    1.5902906645504529930d-7, -1.4919441249454941275d-7,   &
        &    -5.3779684992094351263d-6, -2.7513862296246223142d-5,  &
        &    -9.7880089725297162007d-6, 7.0787668964515789714d-4,   &
        &    4.6968199862345387583d-3, 1.4745740181663320127d-2,    &
        &    2.0048622219583455723d-2 /
         data (b(i), i = 90, 104) /                                 &
        &    -3.4824483072529265585d-12, 1.5157161810563380451d-12, &
        &    8.5303859696700686144d-12, 3.3455414203743741076d-10,  &
        &    2.0226016353844285376d-9, 5.3128154003266334990d-9,    &
        &    -3.0799322316418042137d-8, -4.4455408272954712128d-7,  &
        &    -2.4293274626893384034d-6, -3.2129079340119038354d-6,  &
        &    5.9225403683075388850d-5, 5.6822962576781683532d-4,    &
        &    2.7152446516406682732d-3, 7.4075873691848838485d-3,    &
        &    9.3044450815739269849d-3 /
         data (b(i), i = 105, 119) /                                &
        &    -2.7683216166276377232d-13, 3.1986676777610155465d-12, &
        &    9.4142986954031445666d-12, 6.7934609179456399334d-11,  &
        &    3.4851529411470029330d-11, -2.5785248508896551557d-9,  &
        &    -2.8310220027112571258d-8, -1.6384131113072271115d-7,  &
        &    -3.2521663350596379097d-7, 4.0381388757622307160d-6,   &
        &    5.1917606978077281001d-5, 3.3420027947470126154d-4,    &
        &    1.3699550623118247094d-3, 3.4405619148342271096d-3,    &
        &    4.1042919106665762794d-3 /
         data (c(i), i = 0, 14) /                                   &
        &    4.5281968025889407937d-12, 1.0806749918195271176d-11,  &
        &    9.6200972728717669027d-11, 5.7214227063625263650d-10,  &
        &    3.6077804282954825099d-9, 2.2465236858536681852d-8,    &
        &    1.3676961264308735230d-7, 7.9561767489531997361d-7,    &
        &    4.3014380065615550573d-6, 2.0921713905550285590d-5,    &
        &    8.8079183950590176926d-5, 3.0549414408830252064d-4,    &
        &    8.1295715613927890473d-4, 1.4679809476357079195d-3,    &
        &    1.3439197177355090057d-3 /
         data (c(i), i = 15, 29) /                                  &
        &    7.6019964430402432637d-13, -2.2616198599158271190d-13, &
        &    1.7904450823779000744d-11, 9.1467054855312232717d-11,  &
        &    7.1378582044879519122d-10, 4.9925255415445769102d-9,   &
        &    3.3767315471315546644d-8, 2.1350774539167751457d-7,    &
        &    1.2314353082655232903d-6, 6.2918685053670619181d-6,    &
        &    2.7493229298777000013d-5, 9.8085825401369821771d-5,    &
        &    2.6670282677770444935d-4, 4.8967895428135985381d-4,    &
        &    4.5418248688489697144d-4 /
         data (c(i), i = 30, 44) /                                  &
        &    9.4180115230375147213d-14, 7.5943117003734061145d-14,  &
        &    3.0335730243874287654d-12, 2.0202796115462268051d-11,  &
        &    1.6839020189186971198d-10, 1.2907875663127201526d-9,   &
        &    9.3547676125865798920d-9, 6.2471974110281880722d-8,    &
        &    3.7585985422997380441d-7, 1.9838348288114906484d-6,    &
        &    8.8884862203671982034d-6, 3.2333259238682810218d-5,    &
        &    8.9266668913380400243d-5, 1.6589185669844051903d-4,    &
        &    1.5536921180500113394d-4 /
         data (c(i), i = 45, 59) /                                  &
        &    1.5425475332301107271d-14, 2.8674534590132490434d-14,  &
        &    6.5078462279160216936d-13, 5.0939757793961391211d-12,  &
        &    4.4979837460748975520d-11, 3.6662925847520171711d-10,  &
        &    2.7848878755089582413d-9, 1.9298120059339477820d-8,    &
        &    1.1950323861976892013d-7, 6.4513432758147478287d-7,    &
        &    2.9422095033982461936d-6, 1.0854433321174584937d-5,    &
        &    3.0307433185818899481d-5, 5.6840981443065017850d-5,    &
        &    5.3637016379451945253d-5 /
         data (c(i), i = 60, 74) /                                  &
        &    3.1077953698439839352d-15, 8.6899496170729520378d-15,  &
        &    1.6258562067326054104d-13, 1.4104842571366761537d-12,  &
        &    1.3019455544084110747d-11, 1.1070466372863950239d-10,  &
        &    8.6890603844230597917d-10, 6.1793722175049967488d-9,   &
        &    3.9058865943755615801d-8, 2.1432806981070368523d-7,    &
        &    9.9034657762983230155d-7, 3.6925185861895664251d-6,    &
        &    1.0399877577259449786d-5, 1.9644939661550210015d-5,    &
        &    1.8648773453825584597d-5 /
         data (c(i), i = 75, 89) /                                  &
        &    7.2831555285336286457d-16, 2.6077534095895783532d-15,  &
        &    4.4881202059263153495d-14, 4.1674329383944385626d-13,  &
        &    3.9760100480223728037d-12, 3.4835976355351183010d-11,  &
        &    2.7993254212770249700d-10, 2.0286513276830758107d-9,   &
        &    1.3018343087118439152d-8, 7.2315927974997999365d-8,    &
        &    3.3750708681924201599d-7, 1.2688020879407355571d-6,    &
        &    3.5980954090811587848d-6, 6.8358260635246667316d-6,    &
        &    6.5208606745808860557d-6 /
         data (c(i), i = 90, 104) /                                 &
        &    1.9026412343503745875d-16, 8.0073765508732553766d-16,  &
        &    1.3245754278055523992d-14, 1.2885201653055058502d-13,  &
        &    1.2600129301230402587d-12, 1.1283306843147549277d-11,  &
        &    9.2261481309646814329d-11, 6.7812033168299846818d-10,  &
        &    4.4020645304595102132d-9, 2.4685719238301517679d-8,    &
        &    1.1611886719473112509d-7, 4.3940380936523135466d-7,    &
        &    1.2529878285546791905d-6, 2.3917218527087570384d-6,    &
        &    2.2907574647671878160d-6 /
         data (c(i), i = 105, 119) /                                &
        &    5.3709522135744366512d-17, 2.5239221050372845433d-16,  &
        &    4.0933147145899083360d-15, 4.1152784247617592367d-14,  &
        &    4.0998840572769381012d-13, 3.7319354625807158852d-12,  &
        &    3.0921671702920868014d-11, 2.2975898538634445343d-10,  &
        &    1.5049754445782364328d-9, 8.5030864719789148982d-9,    &
        &    4.0250559391118423810d-8, 1.5312755642491878591d-7,    &
        &    4.3865020375297892208d-7, 8.4059737392822153101d-7,    &
        &    8.0785884122023473319d-7 /
         data (d(i), i = 0, 7) /                                    &
        &    9.2371554649979581914d-14, -2.3111336195788410887d-12, &
        &    5.7728710326649832559d-11, -1.8002298130091372598d-9,  &
        &    7.6810375010517145638d-8, -5.2669973752193823306d-6,   &
        &    1.0112634961227401357d-3, 1.6180215937964160466d-1 /
         data (d(i), i = 8, 15) /                                   &
        &    6.1381146507252683381d-14, -2.1034499679806301862d-12, &
        &    5.7090233460448415278d-11, -1.7990724350642330817d-9,  &
        &    7.6809056078388019946d-8, -5.2669964425290062357d-6,   &
        &    1.0112634957478283390d-3, 1.6180215937970716383d-1 /
         data (d(i), i = 16, 23) /                                  &
        &    4.2458150578401296419d-14, -1.8435733128339016981d-12, &
        &    5.5534955081564656595d-11, -1.7938162188526358466d-9,  &
        &    7.6798230945934117807d-8, -5.2669828728791921259d-6,   &
        &    1.0112634861753356559d-3, 1.6180215938263409582d-1 /
         data (d(i), i = 24, 31) /                                  &
        &    3.0314798962267007518d-14, -1.5915009905364214455d-12, &
        &    5.3275907427402047438d-11, -1.7824862013841369751d-9,  &
        &    7.6763890356447075810d-8, -5.2669199860465945909d-6,   &
        &    1.0112634217687349189d-3, 1.6180215941108227283d-1 /
         data (d(i), i = 32, 39) /                                  &
        &    2.2211515002229271212d-14, -1.3664088221521734796d-12, &
        &    5.0585177270502341602d-11, -1.7645432205894533462d-9,  &
        &    7.6691805594577373698d-8, -5.2667455286976269634d-6,   &
        &    1.0112631862810974580d-3, 1.6180215954783127877d-1 /
         if (x .lt. 0.8d0) then
             t = x * x
             y = (((((((a(0) * t + a(1)) * t + &
        &        a(2)) * t + a(3)) * t + a(4)) * t + &
        &        a(5)) * t + a(6)) * t + a(7)) * x
             y = (((((((a(8) * t + a(9)) * t + &
        &        a(10)) * t + a(11)) * t + a(12)) * t + &
        &        a(13)) * t + a(14)) * t + a(15)) / x + &
        &        y * log(x)
         else if (x .lt. 5.5d0) then
             v = 3.0d0 / x
             t = x - v
             k = int(t + 3)
             t = (k - 2) - t
             k = k * 15
             y = ((((((((((((((b(k) * t + b(k + 1)) * t + &
        &        b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + &
        &        b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + &
        &        b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + &
        &        b(k + 11)) * t + b(k + 12)) * t + b(k + 13)) * t + &
        &        b(k + 14)) * v
         else if (x .lt. 12.5d0) then
             k = int(x)
             t = (k + 1) - x
             k = 15 * (k - 5)
             y = (((((((((((((c(k) * t + c(k + 1)) * t + &
        &        c(k + 2)) * t + c(k + 3)) * t + c(k + 4)) * t + &
        &        c(k + 5)) * t + c(k + 6)) * t + c(k + 7)) * t + &
        &        c(k + 8)) * t + c(k + 9)) * t + c(k + 10)) * t + &
        &        c(k + 11)) * t + c(k + 12)) * t + c(k + 13)) * t + &
        &        c(k + 14)
         else
             t = 60 / x
             k = 8 * (int(t))
             y = (((((((d(k) * t + d(k + 1)) * t + &
        &        d(k + 2)) * t + d(k + 3)) * t + d(k + 4)) * t + &
        &        d(k + 5)) * t + d(k + 6)) * t + d(k + 7)) * &
        &        sqrt(t) * exp(-x)
         end if
         bessk1 = y
         end function bessk1

  end subroutine cc_interaction_on_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module cutoff_coulomb
