! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!===================================================================!
!                  Parallel multigrid solver module                 !
!-------------------------------------------------------------------!

module multigrid

  ! jd: Please do not remove these dependencies
  use comms
  use constants
  use timer
  use utils

  implicit none
  private

  !public :: multigrid_solver
  public :: dummy

contains

  subroutine dummy
  end subroutine dummy

end module multigrid
