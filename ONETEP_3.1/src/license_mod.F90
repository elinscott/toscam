!=============================================================================!
!                          L I C E N S E                                      !
!=============================================================================!
!                                                                             !
!-----------------------------------------------------------------------------!
! Module for checking licenses:                                               !
!                                                                             !
!-----------------------------------------------------------------------------!
! Written by Victor Milman, v0.1, 15 November 2004                            !
!=============================================================================!

module license
  use comms
  use constants
  implicit none      
#ifdef ACCELRYS
#ifndef PAPERLICENSE
#define LIC_BAD_EXIT_STATUS 131
#define LICENSE_F_INTERFACE
#include "MS_ModelingLicensePolicy.h"
#endif
#endif

  private       
  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: license_checkout                 
  public :: license_checkin
  public :: license_setup_traps
  public :: license_sendHeartbeat
  interface
     function signal_handler(sig)
        integer :: signal_handler
        integer, intent(in) :: sig
     end function signal_handler
  end interface

  !---------------------------------------------------------------------------!
  !                       P u b l i c    V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  integer, public, parameter :: LIC_SUCCESS=+1 ! Value returned on success
  integer, public, parameter :: LIC_FAILURE=-1 ! Value returned on failure

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  integer, save :: lic_handle
  integer, save :: VA_ID

contains
  subroutine license_checkout(num_licenses,status)
    !=========================================================================!
    ! Checks out required number of licenses                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  num_licenses, in, required number of licenses                          !
    ! status, out                                                             !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !  lic_handle: stored to be used at the end for checking in licenses      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !  comms                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    ! positive number of licenses                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v1.0, 15.11.2004                              !
    !=========================================================================!
    implicit none

    integer, intent(in) :: num_licenses
    integer, intent(out) :: status

    integer :: mystatus

    status = LIC_SUCCESS
#ifdef ACCELRYS
    if (num_licenses<=0) then
       write (*,'(a)') 'Error in license_checkout: non-positive number of licenses'
       call comms_abort
    end if

    ! create a list of PIDs - needed to kill from Perl script later
    call reportpid()

    ! True if we are on the root node
    if (.not.pub_on_root) return
#ifndef PAPERLICENSE
    ! Only check licences and do io when on root node

    VA_ID = GET_VA_ID_F()
    mystatus = SET_VA_F(VA_ID,"System.Heartbeats", "LS_LICENSE_MANUAL")
    if (mystatus /= 1) status = LIC_FAILURE

    lic_handle = VALID_LICENSE_MP_VA_F(K_onetep,num_licenses,VA_ID)
    if ( lic_handle .eq. 0 ) then
       write(stdout,'(a)') 'Licensing Error !'
       write(stdout,'(a)')                                         &
            &       'Error: Manual heartbeat setup for MS_onetep license failed'
       mystatus = LIC_BAD_EXIT_STATUS 
       status = LIC_FAILURE
    else
       write(stdout,'(a)')                                         &
            &       'License checkout of MS_onetep successful'
       write(stdout,*)
    end if
#endif
#endif

    return

  end subroutine license_checkout



  subroutine license_checkin(status)
    !=========================================================================!
    ! Checks in licenses                                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! status, out                                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v1.0, 15.11.2004                              !
    !=========================================================================!

    implicit none

    integer, intent(out) :: status

    integer :: mystatus

    status = LIC_SUCCESS

    ! True if we are on the root node
    if (.not.pub_on_root) return

    ! Only checkin licenses whan it's an Accelrys build with electronic licensing
#ifdef ACCELRYS
#ifndef PAPERLICENSE
    mystatus = RELEASE_VA_ID_F(VA_ID)
    if (mystatus /= 1) status = LIC_FAILURE

    mystatus = RELEASE_VALID_LICENSE_F(lic_handle)
    if (mystatus /= 0) status = LIC_FAILURE
#endif
#endif

    return

  end subroutine license_checkin

  subroutine license_setup_traps()
    !=========================================================================!
    ! Sets up traps for signals (SIGTERM, etc.)                               !
    ! Needed to exit gracefully when killed                                   !
    ! Required to cope with SGI: manual kill of a process also kills the      !
    ! parent, and so takes down Materials Studio Gateway                      !
    !-------------------------------------------------------------------------!
    ! System function SIGNAL is used; On WIN32 (CVF compiler) it comes from   !
    ! DFPORT library, on UNIX it is declared as external.                     !
    !-------------------------------------------------------------------------!
    ! Arguments: none                                                         !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v1.0, 15.11.2004                              !
    !=========================================================================!

! cks: 2005_8_5: Comment out everything in order to compile on Opteron SUN V40z 
#ifdef ACCELRYS
#ifdef WIN32
    use DFPORT, only : signal
#endif

    implicit none

#ifndef WIN32
    integer,external :: signal
#endif

    integer :: handler

    ! Ideally we would like access to the values of SIGINT and SIGTERM,
    ! however there does not seem to be OS-independent Fortran access to
    ! to those values.
    handler = signal(2,signal_handler,-1)   ! SIGINT
    handler = signal(15,signal_handler,-1)  ! SIGTERM
#endif

    return
  end subroutine license_setup_traps

  subroutine license_sendHeartbeat(status)
    !=========================================================================!
    ! Send manual heartbeats                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! status, out                                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v1.0, 15.10.2004                              !
    !=========================================================================!
    implicit none
    integer,intent(out) :: status

    integer :: idum

    status = LIC_SUCCESS

    ! True if we are on the root node
    if (.not.pub_on_root) return

#ifdef ACCELRYS
#ifndef PAPERLICENSE
    idum = SEND_HEARTBEAT_F( lic_handle )
    if (idum /= 0) then
        write(stdout,'(a)') 'Licensing Error !'
        write(stdout,'(a)') 'Send heartbeat for MS_onetep failed'
        status = LIC_BAD_EXIT_STATUS 
    endif
#endif
#endif

    
    return
  end subroutine license_sendHeartbeat

end module license

integer function signal_handler(sig)

    !=========================================================================!
    ! Sets up traps for signals (SIGTERM, etc.)                               !
    ! Needed to exit gracefully when killed                                   !
    ! Required to cope with SGI: manual kill of a process also kills the      !
    ! parent, and so takes down Materials Studio Gateway                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! sig - dummy argument                                                    !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v1.0, 11.09.2003                              !
    ! Updated by Keith Refson,  v1.1, 24.04.2004                              !
    !=========================================================================!
    implicit none
    integer,intent(in) :: sig

    write(*,'(a)') 'Trapped SIGINT or SIGTERM. Exiting...'

    stop
    
    signal_handler = 1

    return
end function signal_handler
