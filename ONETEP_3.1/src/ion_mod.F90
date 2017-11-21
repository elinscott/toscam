! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by Chris-Kriton Skylaris
!
!   January 2001
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module ion


  use constants, only: dp
  use geometry, only: POINT

  implicit none

  private

  ! cks: this structure defines quantities which are all associated with
  !      a particular atom.
  type ELEMENT
     ! cks: cartesian coordinates of atom in a.u.
     type(POINT)      :: centre
     ! cks: chemical symbol of atom
     character(len=2) :: symbol
     ! cks: name of pseudopotential (file) on atom
     character(len=64) :: pseudo_name
     ! ndmh: names of NGWF sets for valence NGWFs, conduction NGWFs, and
     ! ndmh: auxiliary basis
     character(len=64) :: ngwf_set
     character(len=64) :: cond_ngwf_set
     character(len=64) :: aux_ngwf_set
     ! ndmh: name of core wavefunction file for this atom
     character(len=64) :: core_wf_name
     ! cks: atomic number and number of NGWF functions on atom
     integer          :: atomic_number,nfunctions
     ! ndmh: number of conduction NGWFs, radius of conduction NGWFs
     integer :: nfunctions_cond
     real(kind=DP) :: radius_cond
     ! ndmh: number of auxiliary NGWFs
     integer :: nfunctions_aux
     ! qoh: "species number" to which element belongs (in %block species)
     integer          :: species_number
     ! qoh: "pseudopotential species number" of element
     integer          :: pspecies_number
     ! ndmh: species id string
     character(len=4) :: species_id
     ! cks: charge of the ionic core of a particular element.
     ! pa: changed from integer to real
     real(kind=DP)          :: ion_charge
     ! cks: number of non-local projectors on atom
     integer :: nprojectors
     ! ndmh: number of PAW partial waves on atom
     integer :: npawpws
     ! ndmh: number of PAW core orbitals on atom
     integer :: ncorewfs
     ! cks: NGWF region radius for atom in a.u.
     ! cks: Note that when a halo is defined this is the
     ! cks: sum of the NGWF radius and the halo length
     real(kind=DP) :: radius
     ! cks: max core radius for any projector of this atom in a.u.
     real(kind=DP) :: max_core_radius
     ! ndmh: max core wavefunction radius for any core wf of this atom in a.u.
     real(kind=DP) :: max_core_wf_radius
     ! aam: ionic constraint information
     character(len=5) :: ion_constraint_type
     real(kind=dp)    :: ion_constraint(3)
     ! cks: selective NGWF plotting information
     logical          :: ngwf_plot
     ! smmd: ionic velocity
     real(kind=DP)    :: ion_velocity(3)
     ! smmd: group of atoms to which element belongs
     integer          :: group_id

  end type ELEMENT

  ! cks: public type definitions
  public :: ELEMENT


end module ion
