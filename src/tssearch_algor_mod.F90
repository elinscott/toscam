! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                    T S S E A R C H _ A L G O R                              !
!=============================================================================!
!                                                                             !
! $Id: tssearch_algor.F90                                                     !
!                                                                             !
!-----------------------------------------------------------------------------!
! This module contains the transition state search algorithms                 !
!                                                                             !
!-----------------------------------------------------------------------------!
! Written from  "Module Specification for New Plane Wave Code",  M. Segall,   !
!    P. Lindan, M. Probert, C. Pickard, P. Hasnip, S. Clark and M. Payne      !
!                           Copyright 1999/2000                               !
!-----------------------------------------------------------------------------!
! Written by Niri Govind, v0.1, 27 Nov 2001                                   !
!-----------------------------------------------------------------------------!
! modification information                                                    !
!=============================================================================!

module tssearch_algor

  !use tssearch_utils

  implicit none                                 !Impose strong typing

  private                                       !Everything is private ...

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: tssearch_algor_LST
  public :: tssearch_algor_QST
  public :: tssearch_algor_LST_Test
  public :: tssearch_algor_QST_Test
  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  ! Define the kind parameter for double precison
  integer, parameter :: dprec=kind(1.0d0)

  real(kind=dprec), parameter :: rcut = 20.0_dprec

  !tolerance information
  real(kind=dprec), parameter :: eps_tol = 1e-4_dprec
  real(kind=dprec), parameter :: eps_lstqstbracket = eps_tol
  real(kind=dprec), parameter :: eps_lstqstbrent = eps_tol
  real(kind=dprec), parameter :: eps_linemin = eps_tol
  real(kind=dprec), parameter :: eps_projmin = 1e-3_dprec

  !optimizer information
  real(kind=dprec), parameter :: stpmax_varmet = 5.0_dprec
  real(kind=dprec), parameter :: stpmax_linemin = 3.0e-1_dprec
  real(kind=dprec), parameter :: stpmin_linemin = 1.0e-10_dprec
  real(kind=dprec), parameter :: fallmax = 5.0e-2_dprec
  real(kind=dprec), parameter :: brac_size = 2.0e-1_dprec
  real(kind=dprec), parameter :: stpscal = 1.0_dprec
  real(kind=dprec), parameter :: zero = 1.0e-8_dprec

  !declare pi locally
  !the other constants are not used
  real(kind=dprec), parameter :: pi= 3.141592653589793238462643383279502884197_dprec

  !monitoring information
  character(len=3), public,save :: tsstat
  integer,save :: maxloop
  real (kind=dprec), allocatable, dimension(:),save  :: gwrite
  real (kind=dprec), allocatable, dimension(:),save  :: gwrite_reac
  real (kind=dprec), allocatable, dimension(:),save  :: gwrite_prod
  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

contains

  subroutine tssearch_algor_LST(mtsmethod,nat,ndim,period,impose,&
                            x0,g0,ereactant,x1,g1,eproduct,&
                            xm,gm,elst,xclean1,xclean2,nfix,idfix,&
                            tsfound,neg,output_file,&
                            tssearch_cg_max_iter,tssearch_force_tol,&
                            tssearch_disp_tol)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, Dec 4 2001                                !
    !=========================================================================!
    use tssearch_utils, only: energy_conv, energy_label_ts, force_conv, &
         root_write, ts_stdout, force_label_ts, tssearch_utils_io_abort, &
         tssearch_utils_cart_to_frac, tssearch_utils_energy_gradient, &
         tssearch_utils_write, tssearch_utils_write_forces, &
         tssearch_utils_write_summary, tssearch_utils_get_images

    implicit none

    ! Arguments
    integer       :: mtsmethod,nat,ndim,neg,nfix
    integer       :: idfix(ndim)
    real(kind=dprec) :: x0(ndim),g0(ndim),ereactant
    real(kind=dprec) :: x1(ndim),g1(ndim),eproduct
    real(kind=dprec) :: xclean1(ndim),xclean2(ndim)
    real(kind=dprec) :: xm(ndim),gm(ndim),elst
    logical       :: tsfound,impose,period
    character(len=*),  intent(in)    :: output_file
    integer :: tssearch_cg_max_iter
    real(kind=dprec) :: tssearch_force_tol,tssearch_disp_tol
    
    ! Local variables
    integer       :: ii
    integer       :: inloop,maxinloop
    integer       :: numextr1,numextr0
    integer       :: ierr,i,istatus

    real(kind=dprec) :: tolsrh
    real(kind=dprec) :: s0dotg2,disp_rms,s02
    real(kind=dprec) :: p(ndim),xopt(ndim),gopt(ndim)
    real(kind=dprec) :: s0dotg,p_frac(ndim)
    real(kind=dprec) :: xmin1(ndim),xmin2(ndim)
    real(kind=dprec) :: fmin1,fmin2,fm
    real(kind=dprec) :: lower_bound
    real(kind=dprec) :: beta_cg,s0(ndim),tang(ndim),h0(ndim),sg,dgdf(ndim),hg
    real(kind=dprec) :: g_tang,angle,f_old,eta,dedf
    real(kind=dprec) :: pm,g_rms,g2,f_move,g2_old,x_move
    real(kind=dprec) :: x_old(ndim),g_old(ndim),angmax
    real(kind=dprec) :: grdmin
    real(kind=dprec) :: extremum_sep,rmsval
    real(kind=dprec) :: a1old,fa10ld,tol_brent,flopt0,fopt0,fopt1
    real(kind=dprec) :: xmax1(ndim),xmax2(ndim),xx(ndim),g(ndim)
    real(kind=dprec) :: gmax1(ndim),gmax2(ndim)
    real(kind=dprec) :: s(ndim),flopt1,flopt
    real(kind=dprec) :: f0, f1, fopt
    real(kind=dprec) :: tol, pnorm
    real(kind=dprec) :: a, b, c, fa, fb, fc
    real(kind=dprec) :: a0, b0, c0, fa0, fb0, fc0
    real(kind=dprec) :: a1, b1, c1, fa1, fb1, fc1
    real(kind=dprec) :: bracket_size_f,bracket_size_b
    real(kind=dprec) :: gdum(ndim)
    real(kind=dprec) :: Emax_old,out_of_bounds
    real(kind=dprec) :: fvalue,fvm
    
    logical       :: qstcycle,converged,polak,ldouble

    !====================================================!
    !set necessary variables                             !
    !====================================================!
    maxinloop           = tssearch_cg_max_iter
    grdmin              = tssearch_force_tol
    tolsrh              = tssearch_disp_tol
    ldouble             = .FALSE.
    out_of_bounds       = 1.0e+20_dprec
    tol                 = tolsrh
    Emax_old            = out_of_bounds

    !lst does not require pm
    pm                  = 0.0_dprec      
    eta                 = fallmax
    converged           = .FALSE.
    bracket_size_f      = brac_size
    bracket_size_b      = 0.5_dprec*brac_size
    extremum_sep        = 0.2_dprec
    qstcycle            = .FALSE.
    polak               = .TRUE.
    tsfound             = .FALSE.

    tsstat              = 'LST'
    maxloop             = 0
    lower_bound         = 0.01_dprec
    fopt                = 0.0_dprec
    !======================================================!

    !Write method header
    if (root_write) then
     write(ts_stdout,'(a)') '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(ts_stdout,'(a)') '+++++++++++++++++++++++ Linear synchronous transit method ++++++++++++++++++++++'
     write(ts_stdout,'(a)') '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(ts_stdout,*)
    end if

    !bootstrap superposition and periodicity
    if (period) then
      impose = .FALSE.
    end if

    !allocate gwrite,gwrite_reac,gwrite_prod
    if (mtsmethod /= 3) then
    allocate(gwrite(ndim),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating gwrite in tssearch_algor_LST')
    allocate(gwrite_reac(ndim),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating gwrite_reac in tssearch_algor_LST')
    allocate(gwrite_prod(ndim),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating gwrite_prod in tssearch_algor_LST')
    end if 

    !fix atom and superpose incompatibility
    if (.not.period .and. impose .and. nfix.gt.0) then
       ierr = 1
       if (root_write) then
         write(ts_stdout,'(a)') 'Error: Molecular superposition and fix atoms are incompatible'
       end if
       if (ierr/=0) call tssearch_utils_io_abort('Error detected in tssearch_algor_LST')
    end if

    !calculate p vector, check structures and setup tolerance for maximum search
    pnorm = 0.0_dprec
    do ii = 1,ndim
      p(ii) = x1(ii)-x0(ii)
      pnorm = pnorm + p(ii)*p(ii)
    end do
    pnorm=dsqrt(pnorm)

    !print out difference vector
    !first convert to fractional
    call tssearch_utils_cart_to_frac(p,p_frac)
    if (root_write) then
      write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(ts_stdout,'(a,a)') ' Difference vector: Fractional components '
      write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(ts_stdout,666)(p_frac(i),i=1,ndim)
      write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(ts_stdout,*)
    end if

    !check search tolerance
    call tssearch_algor_dflmin(ndim,p,tolsrh,tol)

    !set image information
    if (period) then
      call tssearch_utils_get_images(rcut)
      if (root_write) then
        write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(ts_stdout,'(a)') ' Periodic cell setup complete'
        write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(ts_stdout,*)
      end if
    end if

    !scan test
    !ngeom = 5
    !pm = 0.0_dprec
    !call tssearch_algor_scan_test(nat,ndim,period,ngeom,x0,x1,xm,pm, &
    !                     iopt,idfix,xclean1,xclean2,neg,output_file)
    !return

    !do we need to superpose molecular structures ?
    if (.not.period) then
      if (impose) then
         call tssearch_algor_rmscalc(nat,ndim,x0,x1,rmsval)
         if (root_write) then
            write(ts_stdout,'(a)')       ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,'(a,f18.5)') ' RMS deviation of structures before superposition :',rmsval
         end if
         call tssearch_algor_supmol(nat,ndim,x0,x1,rmsval)
         if (root_write) then
           write(ts_stdout,'(a,f18.5)')  ' RMS deviation of structures after superposition  :',rmsval
           write(ts_stdout,'(a)')        ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,*)
         end if
      else
         call tssearch_algor_rmscalc(nat,ndim,x0,x1,rmsval)
         if (root_write) then
           write(ts_stdout,'(a)')       ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,'(a,f18.5)') ' RMS deviation of structures :',rmsval
           write(ts_stdout,'(a)')       ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,*)
         end if
      end if
    end if

    !=================================!
    !reactant structure information   !
    !=================================!

    !calculate the path-coordinate for the reactant
    call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x0,fvalue)

    !calculate the energy-gradient of the reactant
    call tssearch_utils_energy_gradient(ndim,x0,g0,f0)

    !update the energy counter
    neg = neg + 1

    !zero out gradients on fixed atoms
    do i = 1,ndim
       gwrite_reac(i) = g0(i)       !keep a copy of all the gradients to write out
       if (idfix(i).eq.1) then
            g0(i) = 0.0_dprec 
       end if
    end do

    !store reactant information
    ereactant = f0
    call tssearch_utils_write(x0,gwrite_reac,ereactant,neg,output_file,&
                       tsstat,maxloop,fvalue)

    !write path-coordinate and energy of structure to output file
    if (root_write) then
     call tssearch_utils_write_summary(5,ts_stdout,f0,0.0_dprec,0.0_dprec,fvalue)
     call tssearch_utils_write_forces(ndim,ts_stdout,x0,gwrite_reac)
    end if

    !=================================!
    !product structure information    !
    !=================================!

    !calculate the path-coordinate for the product
    call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x1,fvalue)
  
    !calculate the energy-gradient of the product
    call tssearch_utils_energy_gradient(ndim,x1,g1,f1)

    !update the energy counter
    neg = neg + 1

    !zero out gradients on fixed atoms
    do i = 1,ndim
       gwrite_prod(i) = g1(i)
       if (idfix(i).eq.1) then
            g1(i) = 0.0_dprec 
       end if
    end do

    !store product information
    eproduct = f1
    call tssearch_utils_write(x1,gwrite_prod,eproduct,neg,output_file,&
                       tsstat,maxloop,fvalue)

    !write path-coordinate and energy of structure to output file
    if (root_write) then
     call tssearch_utils_write_summary(6,ts_stdout,f1,0.0_dprec,0.0_dprec,fvalue)
     call tssearch_utils_write_forces(ndim,ts_stdout,x1,gwrite_prod)
    end if

    !zero out gradients on fixed atoms
    do i=1,ndim
         if (idfix(i) .eq. 1) then
      g0(i)   = 0.0_dprec
      g1(i)   = 0.0_dprec
            gm(i)   = 0.0_dprec
            gopt(i) = 0.0_dprec
            g(i)    = 0.0_dprec
         endif
    enddo

    !prepare temporary arrays
    fmin1 = f0
    fmin2 = f1
    fm = fopt
    do i=1,ndim
         xmin1(i) = x0(i)
         xmin2(i) = x1(i)
         p(i) = xmin2(i) - xmin1(i)
    end do

    if (root_write) then
      write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(ts_stdout,'(a)') ' LST maximum bracketing from reactant'
      write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(ts_stdout,*)
    end if

    !update status flag

    !start bracketing maximum from reactant
    !write(ts_stdout,'(a)') ' Commencing bracketing from reactant'
    !write(ts_stdout,*)
    !write(ts_stdout,'(a)')       ' Bracketing maximum from reactant'
    !write(ts_stdout,'(a,f18.4)') ' Bracket size      : ',bracket_size_f
    !write(ts_stdout,'(a,f18.4)') ' Bracket tolerance : ',tol

    !prepare for maximum bracketing
    numextr0= -1

    !call bracketing routine
    call tssearch_algor_lstqstbracket(nat, ndim, idfix, xmin1, xmin2, fmin1, &
                 a0, b0, c0, 1, xm, pm, 1.0_dprec, tol, 1.0_dprec, bracket_size_f, &
                 fa0, fb0, fc0, numextr0, period, gdum, 0, neg, &
                 xclean1, xclean2, output_file)

    if (numextr0.eq.0) then
           ldouble = .true.
           if (root_write) then
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,'(a)') ' LST maximum bracketing from reactant failed'
             !write(ts_stdout,'(a)') ' Initiating bracketing from product'
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,*)
           end if
    else
           ldouble = .false.
           if (root_write) then
             write(ts_stdout,'(a)')     ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,'(a)')     ' LST maximum bracketing from reactant successful'
             write(ts_stdout,'(a,i6)')  ' Cumulative number of energy/gradient calls :',neg
             write(ts_stdout,'(a)')     ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,*)
           end if

             !initialize a1,b1,c1,fa1,fb1,fc1 if ldouble = .false.
       a1=a0
       b1=b0
       c1=c0
       fa1=fa0
       fb1=fb0
       fc1=fc0

             !bypass bracket from product if bracket from reactant is successful
             goto 1212
    end if

    !bracket from product if first pass failed
    if (ldouble) then
         if (root_write) then
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,'(a)') ' LST maximum bracketing from product'
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,*)
         end if
          !write(ts_stdout,'(a)') ' Commencing bracketing from product'
          !write(ts_stdout,*) 
          !write(ts_stdout,'(a)')       ' Bracketing from product'
          !write(ts_stdout,'(a,f18.4)') ' Bracket size      : ',bracket_size_b
          !write(ts_stdout,'(a,f18.4)') ' Bracket tolerance : ',tol
          !write(ts_stdout,*)

          !prepare for maximum bracketing
          numextr1= -1

          !flip sign of p vector to bracket in the opposite direction
          do ii = 1,ndim
              p(ii) = -p(ii)
          end do

          !call bracketing routine
          call tssearch_algor_lstqstbracket(nat, ndim, idfix, xmin2, xmin1, fmin2, &
             a1, b1, c1, 1, xm, pm, 1.0_dprec, tol, 1.0_dprec, bracket_size_b, &
             fa1, fb1, fc1, numextr1, period, gdum, 0, neg, &
             xclean1, xclean2, output_file)

          !flip back sign of p vector
          do ii = 1,ndim
              p(ii) = -p(ii)
          end do

          !change reference back to x0
          a1old=a1
          fa10ld=fa1
          a1=1.0_dprec-c1
          fa1=fc1
          b1=1.0_dprec-b1
          c1=1.0_dprec-a1old
          fc1=fa10ld

          !failure to bracket from both directions
          if (numextr0.eq.0.and.numextr1.eq.0) then
            if (root_write) then
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,'(a)')         ' Unable to bracket maximum from both endpoints'
             write(ts_stdout,'(a)')         ' Please check endpoints'
             write(ts_stdout,'(a,i6)')      ' Cumulative number of energy/gradient calls: ',neg
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,*)
            end if
             ierr = 1
             if (ierr/=0) call tssearch_utils_io_abort('Error detected in tssearch_algor_LST')
          end if

          !if reactant is successful and product is not
          if (numextr1.eq.0) then
             if (root_write) then
              write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              write(ts_stdout,'(a)') ' Bracketing from product failed'
              write(ts_stdout,'(a)') ' Using bracket from reactant'
              write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              write(ts_stdout,*)
             end if
       a1=a0
       b1=b0
       c1=c0
       fa1=fa0
       fb1=fb0
       fc1=fc0
           else
             if (root_write) then
              write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              write(ts_stdout,'(a)')      ' LST maximum bracketing from product successful'
              write(ts_stdout,'(a,i6)')   ' Cumulative number of energy/gradient calls: ',neg
              write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              write(ts_stdout,*)
             end if
             !initialize a1,b1,c1,fa1,fb1,fc1 in case ldouble = .false.
       a0=a1
       b0=b1
       c0=c1
       fa0=fa1
       fb0=fb1
       fc0=fc1
             goto 1212
           end if
    endif

 1212  continue

    !============================!
    !start maximization stage    !
    !============================!
    !maximize along lst optimized path
    if (root_write) then
     write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(ts_stdout,'(a)') ' Commencing maximum search'
     write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(ts_stdout,*)
    end if
  
    !update status


    !set search tolerance for brent
    tol_brent=dmin1(tol,lower_bound)

    !prepare for brent maximization
    if ((a0.le.a1.and.a1.lt.b0).or.(a1.le.a0.and.a0.lt.b1)) then

    !take the smaller of (a0, a1)
    if (a0.lt.a1)then
       a=a0
       fa=fa0
    else
       a=a1
       fa=fa1
    endif

    !take the larger of (c0, c1)
    if (c0.gt.c1)then
       c=c0
       fc=fc0
    else
       c=c1
       fc=fc1
    endif

    !take the higher energy of (b0, b1)
    if (fb0.gt.fb1)then
       b=b0
       fb=fb0
    else
       b=b1
       fb=fb1
    endif


         !write(ts_stdout,*) 
         !write(ts_stdout,'(a,f18.4)') ' Brent Tolerance : ', tol_brent

         !call maximization routine
         call tssearch_algor_lstqstbrent(nat,ndim,idfix,xmin1,xmin2,1,pm,xm, &
                b,a,c,fb,fa,fc,numextr0,tol_brent,flopt0, &
                fopt0,period,xopt,gopt,neg,xclean1,xclean2, &
                emax_old,output_file)

         !update variables
         flopt1=flopt0
         fopt1=fopt0
         emax_old = fopt0
         elst = fopt0

         !find the rms gradient at the maximum
         g2  = 0.0_dprec
         do i = 1,ndim
          gwrite(i) = gopt(i)
          if (idfix(i).eq.1) then
            gopt(i) = 0.0_dprec
          end if
          g2  = g2 + gopt(i)*gopt(i)
         end do
         elst = fopt0

         !write out long summary
         if (root_write) then
             !calculate path-coordinate for maximum structure
             call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xopt,fvalue)
             call tssearch_utils_write_summary(1,ts_stdout,f0,f1,elst,fvalue)
             call tssearch_utils_write_forces(ndim,ts_stdout,xopt,gwrite)
         end if

         !check if only lst maximum is requested
         if (mtsmethod .eq. 1) then

           !calculate path-coordinate for maximum structure
           call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xopt,fvalue)

           !store maximum structure as penultimate frame
           call tssearch_utils_write(xopt,gwrite,elst,neg,output_file,&
                       tsstat,maxloop,fvalue)

           !calculate path-coordinate for product structure
           call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x1,fvalue)

           !store product as the last frame
           tsstat = 'PRO'
           call tssearch_utils_write(x1,gwrite_prod,eproduct,neg,output_file,&
                       tsstat,maxloop,fvalue)

           !transfer maximum structure to intermediate array and return
           !save current as transition state estimate
           do i = 1,ndim
            xm(i) = xopt(i)
            gm(i) = gopt(i)
           end do

           return

         end if 

        else

          !possibility of several maxima present between endpoints
          !write(ts_stdout,*) 
          !write(ts_stdout,'(a,f18.4)') ' Brent Tolerance :', tol_brent

          !call brent search
          call tssearch_algor_lstqstbrent(nat,ndim,idfix,xmin1,xmin2,1,pm,xm, &
                b0,a0,c0,fb0,fa0,fc0,numextr0,tol_brent,flopt0, &
                fopt0,period,xmax1,gmax1,neg,xclean1,xclean2,&
                emax_old,output_file)

          !call brent search
          call tssearch_algor_lstqstbrent(nat,ndim,idfix,xmin1,xmin2,1,pm,xm, &
                b1,a1,c1,fb1,fa1,fc1,numextr1,tol_brent,flopt1, &
                fopt1,period,xmax2,gmax2,neg,xclean1,xclean2,&
                emax_old,output_file)

          if (root_write) then
           write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,'(a)')         ' At least 2 maxima present between reactant and product'
           write(ts_stdout,'(a,f18.5)')   ' Maximum 1 : Position :        ',flopt0
           write(ts_stdout,'(a,f18.5,1x,a)') ' Maximum 1 : Energy   :        ',fopt0*energy_conv,energy_label_ts
           write(ts_stdout,'(a,f18.5)')   ' Maximum 2 : Position :        ',flopt1
           write(ts_stdout,'(a,f18.5,1x,a)') ' Maximum 2 : Energy   :        ',fopt1*energy_conv,energy_label_ts
           write(ts_stdout,'(a,i6)')      ' Cumulative number of energy/gradient calls: ',neg
           write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,*)
          end if

        end if

        !check distance between maxima in terms of f values
        if(dabs(flopt1-flopt0).gt.extremum_sep) then

           if (root_write) then
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,*)             ' Confirmation of at least 2 maxima along vector P'
            write(ts_stdout,*)             ' Path values of maxima separated by more than :',extremum_sep
            write(ts_stdout,'(a,f18.5,1x,a)') ' Maximum 1 : Energy:           ',fopt0*energy_conv,energy_label_ts
            write(ts_stdout,'(a,f18.5,1x,a)') ' Maximum 2 : Energy:           ',fopt1*energy_conv,energy_label_ts
            write(ts_stdout,'(a)')         ' Please check endpoints'
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           end if

           ierr = 1
           if (ierr/=0) call tssearch_utils_io_abort('Error detected in tssearch_algor_LST')
        end if

        !========================================!
        !start conjugate gradient refinement here!
        !========================================!
        if (root_write) then
         write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(ts_stdout,'(a)') ' Commencing conjugate gradient refinement'
         write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(ts_stdout,*)
        end if

        !transfer variables to new arrays
        do i = 1,ndim
          xx(i) = xopt(i)
          g(i) = gopt(i)
        end do
        flopt = flopt0
        fopt = elst

        !calculate tangent and special path hessian (lst flag set = 1)
        if (root_write) then
         write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(ts_stdout,'(a)') ' Estimating local change at LST maximum'
         write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(ts_stdout,*)
        end if

        call tssearch_algor_tgthess(nat,ndim,idfix,1,xx,xmin1,xmin2,flopt,g,xm,pm, &
                      period,g_rms,tang,g_tang,dgdf,dedf,neg)

        !setup tangent to path at lst maximum
        !dgdf is the rate of change of the gradient along the path
        sg  = 0.0_dprec 
        g2  = 0.0_dprec 
        s02 = 0.0_dprec
        do i = 1,ndim
          if (idfix(i).eq.1) then
            g(i) = 0.0_dprec
          end if
          s0(i) = tang(i)
          sg  = sg + s0(i)*dgdf(i)
          s02 = s02 + s0(i)*s0(i)
          g2  = g2 + g(i)*g(i)
        end do
        g_rms  = dsqrt(g2/nat)

        !setup special path hessian
        do i = 1,ndim
           h0(i) = dgdf(i)/sg
        end do

        !setup initial conjugate direction for refinement
        inloop = 0
        g2 = 0.0_dprec
        hg = 0.0_dprec
        f_move = 1000000.0_dprec
        do i = 1,ndim
          gwrite(i) = g(i)
          if (idfix(i).eq.1) then
            g(i) = 0.0_dprec
          end if
          g2 = g2 + g(i)**2
          hg = hg + h0(i)*g(i)
        end do
        g_rms = dsqrt(g2/nat)

        if (root_write) then
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,'(a,f18.5)') ' Local change at LST maximum: ',dedf
          write(ts_stdout,'(a,f18.5,1x,a)') ' RMS gradient at LST maximum: ',g_rms*force_conv,force_label_ts
          !write(ts_stdout,'(a,f18.5,1x,a)') ' RMS gradient before CG refinement: ',g_rms*force_conv,force_label_ts
          write(ts_stdout,'(a,i6)')      ' Cumulative number of energy/gradient calls: ',neg
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,*)
        end if

        !convergence checks
        !if ((g_rms <= grdmin) .or. (grad_ratio <= eps_projmin) ) then
        if ((g_rms <= grdmin)) then
          converged = .TRUE. 
          tsfound = .TRUE.

          do i = 1,ndim
           xm(i) = xx(i)  
           gm(i) = g(i)
          end do
          elst = fopt

          !calculate path-coordinate of intermediate
          call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xx,fvalue)

          !call energy-gradient one more time on the transition state
          call tssearch_utils_energy_gradient(ndim,xx,gwrite,elst)

          !store intermediate frame as penultimate
          tsstat = 'TS '
          call tssearch_utils_write(xm,gwrite,elst,neg,output_file,&
                       tsstat,maxloop,fvalue)
          fvm = fvalue

          !calculate path-coordinate of product
          call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x1,fvalue)

          !store product frame at the end
          tsstat = 'PRO'
          call tssearch_utils_write(x1,gwrite_prod,eproduct,neg,output_file,&
                       tsstat,maxloop,fvalue)

          !write summary of calculation
          if (root_write) then
           call tssearch_utils_write_summary(3,ts_stdout,f0,f1,fopt,fvm)
           call tssearch_utils_write_forces(ndim,ts_stdout,xm,gwrite)
          end if

          return
        end if

        !prepare initial direction
        !project out path direction from optimization space
        do i = 1,ndim
           s(i) = -g(i) + hg*s0(i)
        end do

        !===========================================================!
        !outer conjugate gradient loop starts here                  !
        !===========================================================!
        !perform conjugate direction search to find minimum
        do while(.not.converged)

         tsstat = 'CG '
         maxloop = 0

         inloop = inloop + 1

         f_old = fopt  
         g2_old = g2
         do i = 1,ndim

          gwrite(i) = g(i)
          if (idfix(i).eq.1) then
            g(i) = 0.0_dprec
          end if

          x_old(i) = xx(i)
          g_old(i) = g(i) 
         end do

        if (root_write) then
         write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(ts_stdout,'(a,i6)')  ' Conjugate gradient step:                    ',inloop
         write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(ts_stdout,*)
        end if

         !start line minimization
         istatus = 0 
         angmax = 90.0_dprec

         if (root_write) then
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,'(a)') ' Commencing line minimization along search path'
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         end if

         !find path-coordinate
         call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xx,fvalue)

         !store frame
         call tssearch_utils_write(x_old,gwrite,f_old,neg,output_file,&
                       tsstat,maxloop,fvalue)

         !write path-coordinate and energy of structure to output file
         if (root_write) then
           call tssearch_utils_write_summary(8,ts_stdout,f_old,0.0_dprec,0.0_dprec,fvalue)
         end if

         !call line minimizer
         call tssearch_algor_linemin(ndim,idfix,fopt,g,xx,s,f_move,&
                        angmax,angle,istatus,neg)

         !safeguard: use -s if angle is greater than angmax
         if (istatus .eq. 3) then   
            do i = 1,ndim
              s(i) = -s(i)
            end do

            if (root_write) then
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,'(a)') ' Search vector pointing uphill'
             write(ts_stdout,'(a)') ' Will change sign of vector and perform line search'
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            end if

            !call line minimizer
            call tssearch_algor_linemin(ndim,idfix,fopt,g,xx,s,f_move,&
                        angmax,angle,istatus,neg)
         end if

         if (istatus.eq.1) then
          if (root_write) then
           write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,'(a)') ' Line search successful'
           write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,*)
          end if
         end if

         !gradient and displacement after line minimization
         g2 = 0.0_dprec
         f_move = f_old - fopt
         x_move = 0.0_dprec
         do i = 1,ndim

           gwrite(i) = g(i)

           if (idfix(i).eq.1) then
              g(i) = 0.0_dprec
           end if

           x_move = x_move + (xx(i)-x_old(i))**2
           g2 = g2 + g(i)**2
         end do
         g_rms = dsqrt(g2/nat)
         disp_rms = dsqrt(x_move/nat)

         if (root_write) then
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,'(a,f18.5,1x,a)') ' RMS gradient:                 ',g_rms*force_conv,force_label_ts
          write(ts_stdout,'(a,f18.5,1x,a)') ' Energy:                       ',fopt*energy_conv,energy_label_ts
          write(ts_stdout,'(a,f18.5,1x,a)') ' Energy change:                ',f_move*energy_conv,energy_label_ts
          write(ts_stdout,'(a,i6)')      ' Cumulative number of energy/gradient calls: ',neg
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,*)
         end if

         do i = 1,ndim
           xopt(i) = xx(i)
           gopt(i) = g(i)
           xm(i)   = xx(i)
           gm(i)   = g(i)
         end do
         fm = fopt
         elst = fm
     
         !check convergence condition

         !if ((g_rms <= grdmin) .or. (grad_ratio <= eps_projmin) ) then
         if ((g_rms <= grdmin)) then
            converged = .TRUE.
            tsfound = .TRUE.

            !calculate path-coordinate of transition state
            call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xm,fvalue)  

            !call energy-gradient one more time on the transition state
            call tssearch_utils_energy_gradient(ndim,xm,gwrite,fm)

            !store transition state as penultimate
            tsstat = 'TS '
            call tssearch_utils_write(xm,gwrite,fm,neg,output_file,&
                       tsstat,maxloop,fvalue)

            fvm = fvalue

            !calculate path-coordinate of product
            call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x1,fvalue)   

            !store product as final
            tsstat = 'PRO'
            call tssearch_utils_write(x1,gwrite_prod,eproduct,neg,output_file,&
                       tsstat,maxloop,fvalue)

            !write summary of calculation
            if (root_write) then
             call tssearch_utils_write_summary(3,ts_stdout,f0,f1,fm,fvm)
             call tssearch_utils_write_forces(ndim,ts_stdout,xm,gwrite)
            end if 

            return
         end if

         !do we need to prepare for a qst cycle
         if (.not.converged) then
           s0dotg = 0.0_dprec
           do i = 1,ndim
             s0dotg = s0dotg + s0(i)*g(i)
           end do

           !projection along the path
           s0dotg2 = s0dotg*s0dotg/g2

           if (root_write) then
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,'(a,f12.5)')    ' Projection along local tangent at maximum: ',s0dotg2
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,*)
           end if

           qstcycle = .FALSE.  

           if (inloop.ge.maxinloop)  then
             qstcycle = .TRUE.
             if (root_write) then
              write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              write(ts_stdout,'(a)')        ' Maximum number of CG steps reached'
              write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              write(ts_stdout,*)
             end if
           end if

           if (s0dotg2 .gt. eta) then
             qstcycle = .TRUE.
             if (root_write) then
              write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              write(ts_stdout,'(a)')        ' System may be falling off the ridge'
              write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              write(ts_stdout,*)
             end if
           end if

           if (istatus .ne. 1) then
             qstcycle = .TRUE.
             if (root_write) then
              write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              write(ts_stdout,'(a)')        ' Line minimization unsuccessful'
              write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              write(ts_stdout,*)
             end if
           end if

           if (root_write) then
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,'(a,i6)')       ' Conjugate gradient step:                    ',inloop
            write(ts_stdout,'(a,i6)')       ' Maximum number of conjugate gradient steps: ',maxinloop
            write(ts_stdout,'(a,i6)')       ' Cumulative number of energy/gradient calls: ',neg
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,*)
           end if

           !if we need a qstcycle ?
           if (qstcycle) then    

             !some spade work before transiting
             !save current as transition state estimate
             do i = 1,ndim            
               xopt(i) = xx(i)
               gopt(i) = g(i)
               xm(i)   = xopt(i)
               gm(i)   = g(i)
             end do
             elst = fopt

             !calculate path-coordinate of structure
             call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xm,fvalue)

             !store frame before transiting to qst
             call tssearch_utils_write(xm,gwrite,elst,neg,output_file,&
                       tsstat,maxloop,fvalue)

             !check whether run is just Halgren-Lipscomb or LST/Optimization
             if (mtsmethod .eq. 12) then
   
                !calculate path-coordinate and store transition state as penultimate
                call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xm,fvalue)
                tsstat = 'TS '
                call tssearch_utils_write(xm,gwrite,elst,neg,output_file,&
                       tsstat,maxloop,fvalue)

                fvm = fvalue

                !calculate path-coordinate and store product frame as last
                call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x1,fvalue) 
                tsstat = 'PRO'
                call tssearch_utils_write(x1,gwrite_prod,eproduct,neg,output_file,&
                       tsstat,maxloop,fvalue)

                !write summary of calculation
                if (root_write) then
                 call tssearch_utils_write_summary(4,ts_stdout,f0,f1,elst,fvm)
                 call tssearch_utils_write_forces(ndim,ts_stdout,xm,gwrite)
                end if

                return

             else

                tsstat = 'QST'
                maxloop = 1

                if (root_write) then
                 write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                 write(ts_stdout,'(a)') ' Quadratic synchronous transit cycle required'
                 write(ts_stdout,'(a)') ' Transferring to QST routine'
                 write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                 write(ts_stdout,*)
                end if

                return

             end if            ! MTSMETHOD
           end if              ! QSTCYCLE

           !setup next conjugate search direction: fletcher-reeves or polak-ribiere
           if (root_write) then
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,'(a)') ' Setting up next conjugate search direction'
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,*)
           end if

           hg = 0.0_dprec
           do i = 1,ndim
                hg = hg + h0(i)*g(i)
           end do

           !polak-ribiere or flectcher-reeves
           if (polak) then
               beta_cg = 0.0_dprec   
               do i = 1,ndim
                  beta_cg = beta_cg + g(i) * (g(i)-g_old(i))    ! polak-ribiere
               end do
           else 
               beta_cg = 0.0_dprec   
               do i = 1,ndim
                  beta_cg = beta_cg + g(i) * g(i)               ! fletcher-reeves
               end do
           end if

           beta_cg = beta_cg / g2_old

           !set new direction
           do i = 1,ndim
                s(i) = -g(i) + hg*s0(i) + beta_cg*s(i)
           end do  

         end if                                       ! conjugate gradient convergence check
        end do                                        ! inner conjugate gradient loop

 666    format(3f15.9)

        return
  end subroutine tssearch_algor_LST

  subroutine tssearch_algor_QST(mtsmethod,nat,ndim,period,impose,&
                            x0,g0,ereactant,x1,g1,eproduct,&
                            xm,gm,eqst,xclean1,xclean2,nfix,idfix,&
                            tsfound,neg,output_file,tssearch_qst_max_iter,&
                            tssearch_cg_max_iter,tssearch_force_tol,&
                            tssearch_disp_tol)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, Dec 5 2001                                !
    !=========================================================================!
    use tssearch_utils, only: ts_stdout, energy_conv, energy_label_ts, &
         force_conv, root_write, force_label_ts, beta, tssearch_utils_write, &
         tssearch_utils_cart_to_frac, tssearch_utils_io_abort, &
         tssearch_utils_energy_gradient, tssearch_utils_write_summary, &
         tssearch_utils_write_forces, tssearch_utils_get_images

    implicit none

    ! Arguments
    integer       :: mtsmethod,nat,ndim,neg,nfix
    integer       :: idfix(ndim)
    real(kind=dprec) :: x0(ndim),g0(ndim),ereactant
    real(kind=dprec) :: x1(ndim),g1(ndim),eproduct
    real(kind=dprec) :: xclean1(ndim),xclean2(ndim)
    real(kind=dprec) :: xm(ndim),gm(ndim),eqst
    logical       :: tsfound,impose,period
    character(len=*),  intent(in)    :: output_file
    integer :: tssearch_qst_max_iter,tssearch_cg_max_iter
    real(kind=dprec) :: tssearch_force_tol,tssearch_disp_tol

    
    ! Local variables
    integer       :: ii
    integer       :: inloop,maxinloop,maxoutloop
    integer       :: numextr0
    integer       :: ierr,i,istatus
    integer       :: outloop

    real(kind=dprec) :: tolsrh
    real(kind=dprec) :: s0dotg2,disp_rms,s02
    real(kind=dprec) :: p(ndim),xopt(ndim),gopt(ndim)
    real(kind=dprec) :: s0dotg,p_frac(ndim)
    real(kind=dprec) :: xmin1(ndim),xmin2(ndim),gmin1(ndim),gmin2(ndim)
    real(kind=dprec) :: fmin1,fmin2,fm,f_ini,f_fin
    real(kind=dprec) :: lower_bound
    real(kind=dprec) :: beta_cg,s0(ndim),tang(ndim),h0(ndim),sg,dgdf(ndim),hg
    real(kind=dprec) :: g_tang,angle,f_old,eta,dedf
    real(kind=dprec) :: pm,g_rms,g2,f_move,g2_old,x_move
    real(kind=dprec) :: x_old(ndim),g_old(ndim),angmax
    real(kind=dprec) :: grdmin
    real(kind=dprec) :: rmsval
    real(kind=dprec) :: tol_brent,flopt0,fopt0
    real(kind=dprec) :: xx(ndim),g(ndim)
    real(kind=dprec) :: dreac,dprod,s(ndim),flopt
    real(kind=dprec) :: f0, f1, fopt
    real(kind=dprec) :: tol, pnorm
    real(kind=dprec) :: a, b, c, fa, fb, fc
    real(kind=dprec) :: gdum(ndim),xdum(ndim),fdum
    real(kind=dprec) :: Emax_old,out_of_bounds
    real(kind=dprec) :: fvalue,fvm,rescale,uu

    logical       :: converged,polak
    logical       :: newmaxcycle

    !====================================================!
    !set necessary variables                             !
    !====================================================!
    maxoutloop          = tssearch_qst_max_iter
    maxinloop           = tssearch_cg_max_iter
    grdmin              = tssearch_force_tol
    tolsrh              = tssearch_disp_tol
    out_of_bounds       = 1.0e+20_dprec
    tol                 = tolsrh
    Emax_old            = out_of_bounds
    eta                 = fallmax
    converged           = .FALSE.
    polak               = .TRUE.
    tsfound             = .FALSE.

    rescale             = 0.0_dprec
    tsstat              = 'QST'
    maxloop             = 1
    lower_bound         = 0.01_dprec
    !======================================================!

    !Write method header
    if (root_write) then
     write(ts_stdout,'(a)') '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(ts_stdout,'(a)') '+++++++++++++++++++++++ Quadratic synchronous transit method +++++++++++++++++++'
     write(ts_stdout,'(a)') '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(ts_stdout,*)
    end if

    !bootstrap superposition and periodicity
    if (period) then
      impose = .FALSE.
    end if

    !allocate gwrite,gwrite_reac,gwrite_prod
    if (mtsmethod == 3) then
    allocate(gwrite(ndim),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating gwrite in tssearch_algor_QST')
    allocate(gwrite_reac(ndim),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating gwrite_reac in tssearch_algor_QST')
    allocate(gwrite_prod(ndim),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating gwrite_prod in tssearch_algor_QST')
    end if

    !fix atom and superpose incompatibility
    if (.not.period .and. impose .and. nfix.gt.0) then
       ierr = 1
       if (root_write) then 
         write(ts_stdout,'(a)') ' Molecular superposition and fix atoms are incompatible'
       end if
       if (ierr/=0) call tssearch_utils_io_abort('Error detected in tssearch_algor_QST')
    end if

    !calculate p vector, check structures and setup tolerance for maximum search
    !write(ts_stdout,*)
    pnorm = 0.0_dprec
    do ii = 1,ndim
      p(ii) = x1(ii)-x0(ii)
      pnorm = pnorm + p(ii)*p(ii)
    end do
    pnorm=dsqrt(pnorm)

    !print out difference vector
    !first convert to fractional
    call tssearch_utils_cart_to_frac(p,p_frac)
    if (root_write) then
      write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(ts_stdout,'(a,a)') ' Difference vector: Fractional components '
      write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(ts_stdout,666)(p_frac(i),i=1,ndim)
      write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(ts_stdout,*)
    end if

    !check search tolerance
    call tssearch_algor_dflmin(ndim,p,tolsrh,tol)

    !set image information
    !if this a COMPLETELSTQST we don't need to do this again
    !this needs to be done in this routine only for straight-away QST start
    if (mtsmethod == 3) then
      if (period) then
        call tssearch_utils_get_images(rcut)
        if (root_write) then
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,'(a)') ' Periodic cell setup complete'
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,*)
        end if
      end if
    end if

    !scan test
    !calculate the pm value
    !call tssearch_algor_rmscalc(nat,ndim,x0,xm,dreac)
    !call tssearch_algor_rmscalc(nat,ndim,x1,xm,dprod)
    !pm = dreac/(dreac+dprod)
    !ngeom = 5
    !call tssearch_algor_scan_test(nat,ndim,period,ngeom,x0,x1,xm,pm, &
    !                       iopt,idfix,xclean1,xclean2,neg,output_file)
    !return

    !do we need to superpose molecular structures ?
    if (.not.period) then
      if (impose) then
         call tssearch_algor_rmscalc(nat,ndim,x0,x1,rmsval)
         if (root_write) then
           write(ts_stdout,'(a)')       ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,'(a,f18.4)') ' RMS deviation of structures before superposition :',rmsval
         end if
         call tssearch_algor_supmol(nat,ndim,x0,x1,rmsval)
         if (root_write) then
           write(ts_stdout,'(a,f18.4)') ' RMS deviation of structures after superposition  :',rmsval
           write(ts_stdout,'(a)')       ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,*)
         end if
      else
         !write(ts_stdout,'(a)') ' No coordinate superposition performed'
         !write(ts_stdout,*)
         call tssearch_algor_rmscalc(nat,ndim,x0,x1,rmsval)
         if (root_write) then
           write(ts_stdout,'(a)')       ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,'(a,f18.4)') ' RMS deviation of structures :',rmsval
           write(ts_stdout,'(a)')       ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,*)
         end if
      end if
    end if

    !==============================!
    !reactant structure information!
    !==============================!

    !calculate the path-coordinate for the reactant
    call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x0,fvalue)

    !check if we already have the reactant energy
    if (dabs(ereactant) .gt. zero) then
         f0 = ereactant
    else

     !calculate the energy-gradient of the reactant
     call tssearch_utils_energy_gradient(ndim,x0,g0,f0)
     ereactant = f0

     !update the energy counter
     neg = neg + 1

     !zero out gradients on fixed atoms
     do i = 1,ndim
       gwrite_reac(i) = g0(i)
       if (idfix(i).eq.1) then
            g0(i) = 0.0_dprec
       end if
     end do

    endif

    !store reactant information
    call tssearch_utils_write(x0,gwrite_reac,ereactant,neg,output_file,&
                       tsstat,maxloop,fvalue)

    !write path-coordinate and energy of structure to output file
    if (root_write) then
     call tssearch_utils_write_summary(5,ts_stdout,f0,0.0_dprec,0.0_dprec,fvalue)
     call tssearch_utils_write_forces(ndim,ts_stdout,x0,gwrite_reac)
    end if

    !==============================!
    !product structure information !
    !==============================!

    !calculate the path-coordinate for the product
    call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x1,fvalue)

    !check if we already have the product energy
    if (dabs(eproduct) .gt. zero) then
        f1 = eproduct
    else

     !calculate the energy-gradient of the product
     call tssearch_utils_energy_gradient(ndim,x1,g1,f1)
     eproduct  = f1

     !update the energy counter
     neg = neg + 1

     !zero out gradients on fixed atoms
     do i = 1,ndim
       gwrite_prod(i) = g1(i)
       if (idfix(i).eq.1) then
          g1(i) = 0.0_dprec
       end if
     end do

    endif

    !store product information
    call tssearch_utils_write(x1,gwrite_prod,eproduct,neg,output_file,&
                       tsstat,maxloop,fvalue)

    !write path-coordinate and energy of structure to output file
    if (root_write) then
     call tssearch_utils_write_summary(6,ts_stdout,f1,0.0_dprec,0.0_dprec,fvalue)
     call tssearch_utils_write_forces(ndim,ts_stdout,x1,gwrite_prod)
    end if

    !get intermediate structure information

    !superpose the structures if needed
    if (.not.period) then
         if (impose) then
            call tssearch_algor_supmol(nat,ndim,x0,xm,dreac) 
            call tssearch_algor_supmol(nat,ndim,x1,xm,dprod) 
            pm = dreac/(dreac+dprod)
         else
            call tssearch_algor_rmscalc(nat,ndim,x0,xm,dreac)
            call tssearch_algor_rmscalc(nat,ndim,x1,xm,dprod)
            pm = dreac/(dreac+dprod)
         end if
    else
         call tssearch_algor_rmscalc(nat,ndim,x0,xm,dreac)
         call tssearch_algor_rmscalc(nat,ndim,x1,xm,dprod)
         pm = dreac/(dreac+dprod)
    end if

    !calculate the path-coordinate for the intermediate
    call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xm,fvalue)  ! STORE INTERMEDIATE

    !check if we already have the intermediate energy
    if ((dabs(eqst) .gt. zero)) then
         fopt = eqst
         fm = eqst
    else

         call tssearch_utils_energy_gradient(ndim,xm,gm,fm)
         fopt   = fm
         eqst = fm

         !update counter
         neg = neg + 1

         !zero out gradients on fixed atoms
         do i = 1,ndim
          gwrite(i) = gm(i)
          if (idfix(i).eq.1) then
            gm(i) = 0.0_dprec
          end if
         end do

     end if

     !store intermediate information
     call tssearch_utils_write(xm,gwrite,fm,neg,output_file,&
                       tsstat,maxloop,fvalue)

     !write path-coordinate and energy of structure to output file
     if (root_write) then
       call tssearch_utils_write_summary(7,ts_stdout,fm,0.0_dprec,0.0_dprec,fvalue)
       call tssearch_utils_write_forces(ndim,ts_stdout,xm,gwrite)
     end if

     !calculate g_rms and zero out forces on fixed atoms
     g2 = 0.0_dprec
     do i = 1,ndim
       gwrite(i) = gm(i)
       if (idfix(i).eq.1) then
         gm(i) = 0.0_dprec
       end if
       g2 = g2 + gm(i)**2
     end do
     g_rms = dsqrt(g2/nat)

     !check to see if our starting intermediate is a transition state structure by chance
     if (g_rms .le. grdmin) then
           converged = .TRUE.
           tsfound = .TRUE.

           !calculate path-coordinate of trsnsition state
           call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xm,fvalue)

           !call energy-gradient one more time on the transition state
           call tssearch_utils_energy_gradient(ndim,xm,gwrite,fm)

           !store transition state as penultimate frame
           tsstat = 'TS '
           call tssearch_utils_write(xm,gwrite,fm,neg,output_file,&
                       tsstat,maxloop,fvalue)
           fvm = fvalue

           !calculate path-coordinate of product
           call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x1,fvalue)

           !last frame is product
           tsstat = 'PRO'
           call tssearch_utils_write(x1,gwrite_prod,eproduct,neg,output_file,&
                       tsstat,maxloop,fvalue)

           !write summary of calculation
           if (root_write) then
             call tssearch_utils_write_summary(3,ts_stdout,f0,f1,fm,fvm)
             call tssearch_utils_write_forces(ndim,ts_stdout,xm,gwrite)
           end if 

           return

     end if           ! grdmin

     !zero out gradients on fixed atoms
     do i=1,ndim
       if (idfix(i) .eq. 1) then
           g0(i)   = 0.0_dprec
           g1(i)   = 0.0_dprec
           gm(i)   = 0.0_dprec
           gopt(i) = 0.0_dprec
           g(i)    = 0.0_dprec
       endif
     enddo

     !some preparation before starting our maximization 
     fmin1 = f0
     fmin2 = f1
     fm = fopt
     do i=1,ndim
         xmin1(i) = x0(i)
         gmin1(i) = g0(i)
         xmin2(i) = x1(i)
         gmin2(i) = g1(i)
     end do

     !set initial and final f values
     f_ini = 0.0_dprec
     f_fin = 1.0_dprec

     !===========================================================!
     !outer maximization loop starts here                        !
     !===========================================================!
     !perform qst maximization loop

       outloop = 0
 100   continue
       outloop = outloop + 1
       
       tsstat = 'QST'
       maxloop = outloop

       !intermediate structure energy
       fopt = eqst
       fm = eqst

       !setup an initial point for line maximization
       !perturb the intermediate point
       uu = pm + 0.1_dprec

       !optimize the point
       call tssearch_algor_optipnt(nat,ndim,period,uu,xdum,xmin1,xmin2,pm,xm,2,idfix)

       !get the energy and gradient for this point
       call tssearch_utils_energy_gradient(ndim,xdum,gdum,fdum)

       !update energy-gradient counter
       neg = neg + 1

       !zero out gradients on fixed atoms
       do i=1,ndim
        if (idfix(i) .eq. 1) then
          gmin1(i)  = 0.0_dprec
          gmin2(i)  = 0.0_dprec
          gopt(i)= 0.0_dprec
          gm(i)  = 0.0_dprec
          gdum(i) = 0.0_dprec
        endif
       enddo

       !maximize along qst path
       if (root_write) then
        write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(ts_stdout,'(a)') ' Commencing maximum search'
        write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(ts_stdout,*)
       end if

       !set maximization flags
       numextr0= -1                ! maximization flag
       a = f_ini                   ! initial f
       b = uu                      ! intermediate f -> perturbed f
       c = f_fin                   ! final f
       fa = fmin1                  ! energy at initial f
       fb = fdum                   ! energy at intermediate f
       fc = fmin2                  ! energy at final f

       !calculate p vector
       do i=1,ndim
         p(i) = xmin2(i) - xmin1(i)
       end do

       !check search tolerance
       call tssearch_algor_dflmin(ndim,p,tolsrh,tol)
       tol_brent=dmin1(tol,lower_bound)

       !write(ts_stdout,*)
       !write(ts_stdout,'(a,f18.4)') ' Brent tolerance : ', tol_brent

       !call maximization routine
       !maximize along qst path
       call tssearch_algor_lstqstbrent(nat,ndim,idfix,xmin1,xmin2,2,pm,xm, &
                b,a,c,fb,fa,fc,numextr0,tol_brent,flopt0, &
                fopt0,period,xopt,gopt,neg,xclean1,xclean2, &
                emax_old, output_file)

       emax_old = fopt0
       eqst = fopt0

       !transfer variables to new arrays
       do i = 1,ndim
          xx(i) = xopt(i)
          g(i) = gopt(i)
          gwrite(i) = g(i)
       end do
       flopt = flopt0
       fopt = eqst

       !write summary of calculation
       if (root_write) then
          !calculate path-coordinate for maximum structure
          call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xopt,fvalue)
          call tssearch_utils_write_summary(2,ts_stdout,f0,f1,eqst,fvalue)
          call tssearch_utils_write_forces(ndim,ts_stdout,xopt,gwrite)
       end if

       !find the rms gradient at the maximum
       g2  = 0.0_dprec
       do i = 1,ndim
        if (idfix(i).eq.1) then
          g(i) = 0.0_dprec
        end if
        g2  = g2 + g(i)*g(i)
       end do

! 2121  continue

       !===========================================================!
       !conjugate refinement starts here                           !
       !===========================================================!
       if (root_write) then
        write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(ts_stdout,'(a)') ' Commencing conjugate gradient refinement'
        write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(ts_stdout,*)
       end if

       !calculate tangent and special path hessian: lst flag = 2
       if (root_write) then
         write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(ts_stdout,'(a)') ' Estimating local change at QST maximum'
         write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(ts_stdout,*)
       end if

       call tssearch_algor_tgthess(nat,ndim,idfix,2,xx,xmin1,xmin2,flopt,g,xm,pm, &
                period,g_rms,tang,g_tang,dgdf,dedf,neg)

       !setup tangent to path at the qst maximum
       !dgdf is the rate of change of the gradient along the path
       sg  = 0.0_dprec
       g2  = 0.0_dprec
       s02 = 0.0_dprec
       do i = 1,ndim

          if (idfix(i).eq.1) then
            g(i) = 0.0_dprec
          end if

          s0(i) = tang(i)
          sg  = sg + s0(i)*dgdf(i)
          s02 = s02 + s0(i)*s0(i)
          g2  = g2 + g(i)*g(i)
        end do
        g_rms  = dsqrt(g2/nat)

        !setup special path hessian
        do i = 1,ndim
           h0(i) = dgdf(i)/sg
           !write(ts_stdout,*) dgdf(i),sg,h0(i)
        end do

        !setup initial conjugate direction for refinement
        inloop = 0
        g2 = 0.0_dprec
        hg = 0.0_dprec
        f_move = 1000000.0_dprec
        do i = 1,ndim
          if (idfix(i).eq.1) then
            g(i) = 0.0_dprec
          end if
          g2 = g2 + g(i)**2
          hg = hg + h0(i)*g(i)
        end do
        g_rms = dsqrt(g2/nat)

        !write(ts_stdout,*)
        !write(ts_stdout,*) 'hg ',hg

        if (root_write) then
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,'(a,f18.5)') ' Local change at QST maximum: ',dedf
          write(ts_stdout,'(a,f18.5,1x,a)') ' RMS gradient at QST maximum: ',g_rms*force_conv,force_label_ts
          !write(ts_stdout,'(a,f18.5,1x,a)') ' RMS gradient before CG refinement: ',g_rms*force_conv,force_label_ts
          write(ts_stdout,'(a,i6)')    ' Cumulative number of energy/gradient calls: ',neg
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,*)
        end if

        !convergence checks
        !if ((g_rms <= grdmin) .or. (grad_ratio <= eps_projmin) ) then
        if ((g_rms <= grdmin)) then
          converged = .TRUE.
          tsfound = .TRUE.

          do i = 1,ndim
           xm(i) = xx(i)
           gm(i) = g(i)
          end do
          eqst = fopt

          !calculate path-coordinate of intermediate
          call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xx,fvalue)

          !call energy-gradient one more time on the transition state
          call tssearch_utils_energy_gradient(ndim,xx,gwrite,eqst)

          !store intermediate frame as penultimate
          tsstat = 'TS '
          call tssearch_utils_write(xm,gwrite,eqst,neg,output_file,&
                       tsstat,maxloop,fvalue)
          fvm = fvalue

          !calculate path-coordinate of product
          call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x1,fvalue)

          !store product frame at the end
          tsstat = 'PRO'
          call tssearch_utils_write(x1,gwrite_prod,eproduct,neg,output_file,&
                       tsstat,maxloop,fvalue)

          !write summary of calculation
          if (root_write) then
             call tssearch_utils_write_summary(3,ts_stdout,f0,f1,fopt,fvm)
             call tssearch_utils_write_forces(ndim,ts_stdout,xm,gwrite)
          end if 

          return
        end if

        !project out path direction from optimization space
        do i = 1,ndim
           s(i) = -g(i) + hg*s0(i)
        end do

        !===========================================================!
        !inner conjugate gradient loop starts here                  !
        !===========================================================!
        !perform conjugate direction search to find minimum
        do while(.not.converged)

         tsstat = 'CG '

         inloop = inloop + 1

         f_old = fopt
         g2_old = g2
         do i = 1,ndim

           gwrite(i) = g(i)
           if (idfix(i).eq.1) then
             g(i) = 0.0_dprec
           end if

           x_old(i) = xx(i)
           g_old(i) = g(i)
         end do

         if (root_write) then
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,'(a,i6)')  ' Conjugate gradient step:                    ',inloop
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,*)
         end if

         !start line minimization
         istatus = 0
         angmax = 90.0_dprec

         if (root_write) then
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(ts_stdout,'(a)') ' Commencing line minimization along search path'
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         end if

         !calculate the path-coordinate
         call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xx,fvalue)

         !store the frame
         call tssearch_utils_write(x_old,gwrite,f_old,neg,output_file,&
                       tsstat,maxloop,fvalue)

         !write path-coordinate and energy of structure to output file
         if (root_write) then
          call tssearch_utils_write_summary(8,ts_stdout,f_old,0.0_dprec,0.0_dprec,fvalue)
         end if

         !call line minimizer
         call tssearch_algor_linemin(ndim,idfix,fopt,g,xx,s,f_move,&
                        angmax,angle,istatus,neg)

         !safeguard: use -s if angle is greater than angmax
         if (istatus .eq. 3) then
            do i = 1,ndim
              s(i) = -s(i)
            end do

            if (root_write) then
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,'(a)') ' Search vector pointing uphill'
             write(ts_stdout,'(a)') ' Will change sign of vector and perform line search'
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            end if

            !call line minimizer
            call tssearch_algor_linemin(ndim,idfix,fopt,g,xx,s,f_move,&
                           angmax,angle,istatus,neg)
         end if

         if (istatus.eq.1) then
          if (root_write) then
           write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,'(a)') ' Line search successful'
           write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,*)
          end if
         end if

         !gradient and displacement after line minimization
         g2 = 0.0_dprec
         f_move = f_old - fopt
         x_move = 0.0_dprec
         do i = 1,ndim

           gwrite(i) = g(i)
           if (idfix(i).eq.1) then
             g(i) = 0.0_dprec
           end if

           x_move = x_move + (xx(i)-x_old(i))**2
           g2 = g2 + g(i)**2
         end do
         g_rms = dsqrt(g2/nat)
         disp_rms = dsqrt(x_move/nat)

         if (root_write) then
         write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(ts_stdout,'(a,f18.5,1x,a)') ' RMS gradient:                 ',g_rms*force_conv,force_label_ts
         write(ts_stdout,'(a,f18.5,1x,a)') ' Energy:                       ',fopt*energy_conv,energy_label_ts
         write(ts_stdout,'(a,f18.5,1x,a)') ' Energy change:                ',f_move*energy_conv,energy_label_ts
         write(ts_stdout,'(a,i6)')      ' Cumulative number of energy/gradient calls: ',neg
         write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(ts_stdout,*)
         end if

         do i = 1,ndim
           xopt(i) = xx(i)
           gopt(i) = g(i)
           xm(i)   = xx(i)
           gm(i)   = g(i)
         end do
         fm = fopt
         eqst = fm

         !convergence checks
         !if ((g_rms <= grdmin) .or. (grad_ratio <= eps_projmin) ) then
         if ((g_rms <= grdmin)) then
            converged = .TRUE.
            tsfound = .TRUE.

            !path coordinate of transition state
            call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xm,fvalue)

            !call energy-gradient one more time on the transition state
            call tssearch_utils_energy_gradient(ndim,xm,gwrite,eqst)

            !store transition state is penultimate frame
            tsstat = 'TS '
            call tssearch_utils_write(xm,gwrite,eqst,neg,output_file,&
                       tsstat,maxloop,fvalue)

            fvm = fvalue

            !path coordinate of product
            call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x1,fvalue)

            !product is last frame
            tsstat = 'PRO'
            call tssearch_utils_write(x1,gwrite_prod,eproduct,neg,output_file,&
                       tsstat,maxloop,fvalue)

            !write summary of calculation
            if (root_write) then
              call tssearch_utils_write_summary(3,ts_stdout,f0,f1,fm,fvm)
              call tssearch_utils_write_forces(ndim,ts_stdout,xm,gwrite)
            end if 

            return

         else if (outloop.ge.maxoutloop) then

            converged = .TRUE.
            if (root_write) then
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,'(a)')        ' QST cycle limit reached'
            write(ts_stdout,'(a,i6)')     ' QST cycle:                                  ',outloop
            write(ts_stdout,'(a,i6)')     ' Maximum number of QST cycles:               ',maxoutloop
            write(ts_stdout,'(a)')        ' Transition state search failed'
            write(ts_stdout,'(a)')        ' Please check the endpoints'
            write(ts_stdout,'(a,i6)')     ' Cumulative number of energy/gradient calls: ',neg
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,*)
            end if

            !path coordinate of transition state
            call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xm,fvalue)

            !call energy-gradient one more time on structure
            call tssearch_utils_energy_gradient(ndim,xm,gwrite,fm)

            !store transition state is penultimate frame
            tsstat = 'TS '
            call tssearch_utils_write(xm,gwrite,fm,neg,output_file,&
                       tsstat,maxloop,fvalue)

            fvm = fvalue

            !path coordinate of product
            call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x1,fvalue)

            !product is last frame
            tsstat = 'PRO'
            call tssearch_utils_write(x1,gwrite_prod,eproduct,neg,output_file,&
                       tsstat,maxloop,fvalue)

            !write summary of calculation
            if (root_write) then
              call tssearch_utils_write_summary(4,ts_stdout,f0,f1,fm,fvm)
            end if

            return
         end if

         !set up another conjugate direction ? or prepare for a maximization ?
         if (.not.converged) then
           s0dotg = 0.0_dprec
           do i = 1,ndim
              s0dotg = s0dotg + s0(i)*g(i)
           end do
           s0dotg2 = s0dotg*s0dotg/g2

           if (root_write) then
           write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,'(a,f12.5)')    ' Projection along local tangent at maximum: ',s0dotg2
           write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,*)
           end if

           newmaxcycle = .FALSE.
           if (inloop.ge.maxinloop)  then
             newmaxcycle = .TRUE.
             if (root_write) then
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,'(a)')        ' Maximum number of CG steps reached'
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,*)
             end if
           end if

           if (s0dotg2 .gt. eta) then
             newmaxcycle = .TRUE.
             if (root_write) then
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,'(a)')        ' System may be falling off the ridge'
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,*)
             end if
           end if
     
           if (istatus .ne.1) then
             newmaxcycle = .TRUE.
             if (root_write) then
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,'(a)')        ' Line minimization unsuccessful'
             write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
             write(ts_stdout,*)
             end if
           end if

           if (root_write) then
           write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,'(a,i6)')     ' QST cycle                                 : ',outloop
           write(ts_stdout,'(a,i6)')     ' CG step                                   : ',inloop
           write(ts_stdout,'(a,i6)')     ' Maximum number of QST cycles              : ',maxoutloop
           write(ts_stdout,'(a,i6)')     ' Maximum number of conjugate gradient steps: ',maxinloop
           write(ts_stdout,'(a,i6)')     ' Cumulative number of energy/gradient calls: ',neg
           write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(ts_stdout,*)
           end if

           !calculate path-coordinate and write frame
           call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xm,fvalue)
           call tssearch_utils_write(xm,gwrite,fm,neg,output_file,&
                       tsstat,maxloop,fvalue)

           !we need another maximization cycle
           if (newmaxcycle) then

            if (root_write) then
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,'(a)') ' New maximization cycle required'
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,*)
            end if

            tsstat = 'QST'
            maxloop = outloop + 1

            !reset xmin1 and xmin2: move the reference points
            if (rescale.ne.0.0_dprec) then
              if (.not.period) then
                 if (impose) then
                   call tssearch_algor_supmol(nat,ndim,xmin1,xm,dreac) 
                   call tssearch_algor_supmol(nat,ndim,xmin2,xm,dprod) 
                   pm = dreac/(dreac+dprod)
                 else
                   call tssearch_algor_rmscalc(nat,ndim,xmin1,xm,dreac)
                   call tssearch_algor_rmscalc(nat,ndim,xmin2,xm,dprod)
                   pm = dreac/(dreac+dprod)
                 end if
              else
                  call tssearch_algor_rmscalc(nat,ndim,xmin1,xm,dreac)
                  call tssearch_algor_rmscalc(nat,ndim,xmin2,xm,dprod)
                  pm = dreac/(dreac+dprod)
              end if

              f_ini = rescale * pm
              fvalue = f_ini

              !optimize point and calculate the path-coordinate
              call tssearch_algor_optipnt(nat,ndim,period,f_ini,x_old,xmin1,xmin2,pm,xm,2,idfix)
              call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x_old,fvalue)

              do i = 1,ndim
                  xmin1(i) = x_old(i)
              end do

              !calculate the energy and gradient and write record
              call tssearch_utils_energy_gradient(ndim,xmin1,gmin1,fmin1)
              neg = neg + 1
              call tssearch_utils_write(xmin1,gmin1,fmin1,neg,output_file,&
                       tsstat,maxloop,fvalue)

              f_fin = 1.0_dprec - rescale*(1.0_dprec-pm)
              fvalue = f_fin
              call tssearch_algor_optipnt(nat,ndim,period,f_fin,x_old,xmin1,xmin2,pm,xm,2,idfix) 
              call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,x_old,fvalue)

              do i = 1,ndim
                xmin2(i) = x_old(i)
              end do 

              !calculate the energy and gradient and write record
              call tssearch_utils_energy_gradient(ndim,xmin2,gmin2,fmin2)
              neg = neg + 1
              call tssearch_utils_write(xmin2,gmin2,fmin2,neg,output_file,&
                       tsstat,maxloop,fvalue)

              !loop back to maximization loop
              goto 100

            else

              if (.not.period) then
                  if (impose) then
                    call tssearch_algor_supmol(nat,ndim,xmin1,xm,dreac) 
                    call tssearch_algor_supmol(nat,ndim,xmin2,xm,dprod) 
                    pm = dreac/(dreac+dprod)
                  else
                    call tssearch_algor_rmscalc(nat,ndim,xmin1,xm,dreac)
                    call tssearch_algor_rmscalc(nat,ndim,xmin2,xm,dprod)
                    pm = dreac/(dreac+dprod)
                  end if
                  call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xm,fvalue)
                  call tssearch_utils_write(xm,gwrite,fm,neg,output_file,&
                       tsstat,maxloop,fvalue)
              else
                  call tssearch_algor_rmscalc(nat,ndim,xmin1,xm,dreac)
                  call tssearch_algor_rmscalc(nat,ndim,xmin2,xm,dprod)
                  pm = dreac/(dreac+dprod)
                  call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xm,fvalue)
                  call tssearch_utils_write(xm,gwrite,fm,neg,output_file,&
                       tsstat,maxloop,fvalue)
              end if

              !loop back to maximization loop
              goto 100

            end if   !rescale check
           end if    !maximization check

           !setup next conjugate search direction: fletcher-reeves or polak-ribiere
           if (root_write) then
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,'(a)') ' Setting up next conjugate search direction'
            write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(ts_stdout,*)
           end if

           hg = 0.0_dprec
           do i = 1,ndim
                hg = hg + h0(i)*g(i)
           end do

           !polak-ribiere or flectcher-reeves
           if (polak) then
             beta_cg = 0.0_dprec
             do i = 1,ndim
                beta_cg = beta_cg + g(i) * (g(i)-g_old(i))    ! polak-ribiere
             end do
           else
             beta_cg = 0.0_dprec
             do i = 1,ndim
                beta_cg = beta_cg + g(i) * g(i)               ! fletcher-reeves
             end do
           end if

           beta_cg = beta_cg / g2_old

           !set new direction
           do i = 1,ndim
                s(i) = -g(i) + hg*s0(i) + beta*s(i)
           end do

         end if                                    ! conjugate gradient convergence check
        end do                                     ! inner conjugate gradient loop

 666    format(3f15.9)

        return
  end subroutine tssearch_algor_QST

  subroutine tssearch_algor_LST_Test!(mtsmethod,nat,ndim,period,impose,&
                            !x_reac,g_reac,e_reac,x_prod,g_prod,e_prod,&
                            !x_intm,g_intm,e_intm,xclean1,xclean2,nfix,idfix,&
                            !tsfound,neg,output_file,tssearch_qst_max_iter,&
                            !tssearch_cg_max_iter,tssearch_force_tol,&
                            !tssearch_disp_tol)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, Nov 27 2001                               !
    !=========================================================================!
    use tssearch_utils, only: ts_stdout, clength, nzcell, nycell, nxcell, &
         ncell, gamma, blength, alpha, beta, alength, other_images, &
         tssearch_utils_get_images

    implicit none

    ! Arguments
    !integer :: nat,mtsmethod
    !integer :: ndim
    !logical :: period,impose

    !real(kind=dprec) :: x_reac(:)
    !real(kind=dprec) :: g_reac(:)
    !real(kind=dprec) :: e_reac
    !real(kind=dprec) :: x_prod(:)
    !real(kind=dprec) :: g_prod(:)
    !real(kind=dprec) :: e_prod
    !real(kind=dprec) :: x_intm(:)
    !real(kind=dprec) :: g_intm(:)
    !real(kind=dprec) :: e_intm
    !real(kind=dprec) :: xclean1(:),xclean2(:)
    !integer       :: idfix(:)
    !logical       :: tsfound
    !integer       :: neg,nfix
    !character(len=*),  intent(in)    :: output_file
    !integer :: tssearch_qst_max_iter,tssearch_cg_max_iter
    !real(kind=dprec) :: tssearch_force_tol,tssearch_disp_tol

    ! Local variables

    write(ts_stdout,*) 'alpha ',alpha
    write(ts_stdout,*) 'beta ',beta
    write(ts_stdout,*) 'gamma ',gamma
    write(ts_stdout,*) 'alength ',alength
    write(ts_stdout,*) 'blength ',blength
    write(ts_stdout,*) 'clength ',clength

    call tssearch_utils_get_images(rcut)

    write(ts_stdout,*) 'nxcell ',nxcell
    write(ts_stdout,*) 'nycell ',nycell
    write(ts_stdout,*) 'nzcell ',nzcell
    write(ts_stdout,*) 'ncell ',ncell
    write(ts_stdout,*) 'other_images ',other_images

    call tssearch_algor_periodic_test
    call tssearch_algor_message_test

    return
  end subroutine tssearch_algor_LST_Test

  subroutine tssearch_algor_QST_Test!(mtsmethod,nat,ndim,period,impose,&
                            !x_reac,g_reac,e_reac,x_prod,g_prod,e_prod,&
                            !x_intm,g_intm,e_intm,xclean1,xclean2,nfix,idfix,&
                            !tsfound,neg,output_file,tssearch_qst_max_iter,&
                            !tssearch_cg_max_iter,tssearch_force_tol,&
                            !tssearch_disp_tol)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, Nov 27 2001                               !
    !=========================================================================!

    implicit none

    ! Arguments
    !integer :: nat,mtsmethod
    !integer :: ndim
    !logical :: period,impose
    !real(kind=dprec) ::x_reac(:)
    !real(kind=dprec) ::g_reac(:)
    !real(kind=dprec) ::e_reac
    !real(kind=dprec) ::x_prod(:)
    !real(kind=dprec) ::g_prod(:)
    !real(kind=dprec) ::e_prod
    !real(kind=dprec) ::x_intm(:)
    !real(kind=dprec) ::g_intm(:)
    !real(kind=dprec) ::e_intm
    !real(kind=dprec) ::xclean1(:),xclean2(:)
    !integer       :: idfix(:)
    !logical       :: tsfound
    !integer       :: neg,nfix
    !character(len=*),  intent(in)    :: output_file
    !integer :: tssearch_qst_max_iter,tssearch_cg_max_iter
    !real(kind=dprec) :: tssearch_force_tol,tssearch_disp_tol

    ! Local variables

    call tssearch_algor_message_test

    return
  end subroutine tssearch_algor_QST_Test


  !-------------------------------------------------------------------------!
  !              P R I V A T E   R O U T I N E S                            !
  !-------------------------------------------------------------------------!

  subroutine tssearch_algor_jacobi(n,np,a,d,v,b,z)
  !=========================================================================!
  ! Matrix diagonalization of a real symmetric matrix by Jacobi rotations   !
  ! Based upon "Numerical Recipes, Second Edition"                          !
  !-------------------------------------------------------------------------!
  ! Arguments:                                                              !
  !    n, input  ! logical dimension of the matrix to be diagonalized       !
  !   np, input  ! physical dimension of the matrix storage area            !
  !    a, input  ! matrix to be diagonalized                                !
  !                    (only upper triangle and diagonal are required)      !
  !    d, output ! eigenvalues in ascending order                           !
  !    v, output ! eigenvectors of the matrix                               !
  !    b         ! temporary work array                                     !
  !    z         ! temporary work array                                     !
  !-------------------------------------------------------------------------!
  ! Written by Niri Govind,  Sept 01 2001                                   !
  !=========================================================================!
    use tssearch_utils, only : ts_stdout

    implicit none

    ! Arguments

    integer          :: n,np
    real(kind=dprec)    :: a(np,np),d(np),v(np,np),b(np),z(np)

    ! Local variables

    integer          ::  i,j,k,ip,iq,nrot,maxrot
    real(kind=dprec)    ::  sm,tresh,s,c,t,theta,tau,h,g,p

    !initialization
    maxrot = 100
    nrot = 0
    do ip = 1, n
       do iq = 1, n
          v(ip,iq) = 0.0_dprec
       end do
       v(ip,ip) = 1.0_dprec
    end do
    do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0_dprec
    end do

    !perform jacobi rotations
      do i = 1, maxrot
        sm = 0.0_dprec
        do ip = 1, n-1
           do iq = ip+1, n
              sm = sm + abs(a(ip,iq))
           end do
        end do
        if (sm .eq. 0.0_dprec)  goto 10
        if (i .lt. 4) then
           tresh = 0.2_dprec*sm / n**2
        else
           tresh = 0.0_dprec
        end if
        do ip = 1, n-1
           do iq = ip+1, n
              g = 100.0_dprec * abs(a(ip,iq))
              if (i.gt.4 .and. abs(d(ip))+g.eq.abs(d(ip)) &
                     .and. abs(d(iq))+g.eq.abs(d(iq))) then
                  a(ip,iq) = 0.0_dprec
              else if (abs(a(ip,iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h)+g .eq. abs(h)) then
                     t = a(ip,iq) / h
                  else
                     theta = 0.5_dprec*h / a(ip,iq)
                     t = 1.0_dprec / (abs(theta)+sqrt(1.0_dprec+theta**2))
                     if (theta .lt. 0.0_dprec)  t = -t
                  end if
                  c = 1.0_dprec / sqrt(1.0_dprec+t**2)
                  s = t * c
                  tau = s / (1.0_dprec+c)
                  h = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0_dprec
                  do j = 1, ip-1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = ip+1, iq-1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = iq+1, n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g - s*(h+g*tau)
                     a(iq,j) = h + s*(g-h*tau)
                  end do
                  do j = 1, n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g - s*(h+g*tau)
                     v(j,iq) = h + s*(g-h*tau)
                  end do
                  nrot = nrot + 1
              end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0_dprec
         end do
      end do

 10   continue

      !print warning if not converged
      if (nrot .eq. maxrot) then
         write(ts_stdout,*) ' Warning: tssearch_algor_jacobi -Matrix diagonalization not converged'
      end if

      !eigenvalue and eigenvector sorting
      do i = 1, n-1
         k = i
         p = d(i)
         do j = i+1, n
            if (d(j) .lt. p) then
               k = j
               p = d(j)
            end if
         end do
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
               p = v(j,i)
               v(j,i) = v(j,k)
               v(j,k) = p
            end do
         end if
      end do

      return
    end subroutine tssearch_algor_jacobi

    subroutine tssearch_algor_pathcoord(nat,ndim,period,x1,x2,xi,fvalue)
    !=========================================================================!
    ! Calculates the path coordinate of a structure                           !
    ! fvalue = rmsdevfromreactant/(rmsdevfromreactant + rmsdevfromproduct)    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  nat   , input  ! number of atoms                                       !
    !  ndim  , input  ! number of degrees of freedom  ! 3*nat                 !
    !  period, input  ! if periodic model or not                              !
    !  x1    , input  ! endpoint 1                                            !
    !  x2    , input  ! endpoint 2                                            !
    !  xi    , input  ! some intermediate point                               !
    !  fvalue, output ! path coordinate                                       !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Aug 21 2000                                    !
    !=========================================================================!

    implicit none

    ! Arguments

    integer          :: nat,ndim
    logical          :: period
    real(kind=dprec)    :: x1(ndim),x2(ndim),xi(ndim),fvalue
    real(kind=dprec)    :: xdum(ndim)

    ! Local variables

    integer          :: i
    real(kind=dprec)    :: dreac,dprod

    !read in intermediate coordinate
    do i = 1,ndim
        xdum(i)  = xi(i)
    end do

    !superpose molecular coordinates and do not if periodic
    if (.not.period) then
        call tssearch_algor_supmol(nat,ndim,x1,xdum,dreac)
        call tssearch_algor_supmol(nat,ndim,x2,xdum,dprod)
    else
        call tssearch_algor_rmscalc(nat,ndim,x1,xdum,dreac)
        call tssearch_algor_rmscalc(nat,ndim,x2,xdum,dprod)
    end if

    fvalue = dreac/(dreac + dprod)

    return
    end subroutine tssearch_algor_pathcoord

    subroutine tssearch_algor_conjdir(ndim,cgtype,s,g,g_old)
    !=========================================================================!
    ! Calculates the path coordinate of a structure                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  ndim     , input  ! number of degrees of freedom  ! 3*nat              !
    !  cgtype   , input  ! cgtype ! polak-ribiere(1) or fletcher-reeves(2)    !
    !  s        , input  ! direction                                          !
    !  g        , input  ! gradient                                           !
    !  g_old    , input  ! old gradient                                       !
    !  s        , output ! updated direction                                  !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Sept 05 2000                                   !
    !=========================================================================!

    use tssearch_utils, only: tssearch_utils_io_abort
    implicit none

    ! Arguments

    integer          :: ndim,cgtype
    real(kind=dprec)    :: s(ndim),g(ndim),g_old(ndim)

    ! Local variables

    integer          :: i
    real(kind=dprec)    :: beta,g2_old,g2

    ! old norm
    g2_old = 0.0_dprec
    do i = 1,ndim
       g2_old  = g2_old + g_old(i)*g_old(i)
    end do

    beta = 0.0_dprec
    ! Polak-Ribiere factor
    if (cgtype.eq.1) then
       do i = 1,ndim
          beta = beta + g(i) * (g(i)-g_old(i))
       end do
    else if (cgtype.eq.2) then
       ! Fletcher-Reeves factor
       do i = 1,ndim
          beta = beta + g(i) * g(i)
       end do
    else
       call tssearch_utils_io_abort('Unknown conjugate gradient type found in tssearch_algor_conjdir')
    end if

    beta = beta / g2_old

    ! Update direction
    g2 = 0.0_dprec
    do i = 1,ndim
        g2  = g2 + g(i)*g(i)
        s(i) = -g(i) + beta*s(i)
    end do

    return
    end subroutine tssearch_algor_conjdir

    subroutine tssearch_algor_projdir(ndim,hg,s0,s)
    !=========================================================================!
    ! Projects out a specific direction from a vector                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  ndim     , input  ! number of degrees of freedom  ! 3*nat              !
    !  s        , input  ! direction                                          !
    !  hg       , input  ! projecton pre-factor                               !
    !  s0       , input  ! direction to be projected out                      !
    !  s        , output ! updated direction                                  !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  09/2000                                        !
    !=========================================================================!

    implicit none

    ! Arguments

    integer          :: ndim
    real(kind=dprec)    :: hg,s0(ndim),s(ndim)

    ! Local variables

    integer          :: i

    ! Remove direction
    do i = 1,ndim
        s(i) = s(i) + hg*s0(i)
    end do

    return
    end subroutine tssearch_algor_projdir

    function tssearch_algor_wgtfunc(rcut,rr)
    !=========================================================================!
    ! Modified transit function weight                                        !
    ! Based on T. Halgren and H. Lipscomb, Chem. Phys. Lett (1977)            !
    ! It cuts off at rcut                                                     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  rcut     , input  ! cutoff distance                                    !
    !  rr       , input  ! distance between two atoms                         !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  09/2000                                        !
    !=========================================================================!

    implicit none
    
    ! Arguments

    real(kind=dprec) :: rcut,rr
    real(kind=dprec) :: tssearch_algor_wgtfunc

    ! Local variables

    real(kind=dprec) :: r0,rr2,rr4,rr6,rcut4,wgtfunc

    ! Setup some useful quantities
    r0    = 0.0_dprec
    rr2   = rr*rr
    rr4   = rr2*rr2
    rr6   = rr2*rr4
    rcut4 = rcut*rcut*rcut*rcut

    ! Modified Halgren-Lipscomb weight function
    if (rr.le.r0) then
       wgtfunc = (1.0_dprec / rr4   -  1.0_dprec / rcut4)
    else if (rr.gt.r0.and.rr.le.rcut) then
       wgtfunc = 1.0_dprec / rr4   -  1.0_dprec / rcut4
    else !if (rr.gt.rcut) then
       wgtfunc = 0.0_dprec
    end if
    tssearch_algor_wgtfunc = wgtfunc

    return
    end function tssearch_algor_wgtfunc

    subroutine tssearch_algor_supmol(nat,ndim,r1,r2,rmsval)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  09/2000                                        !
    !=========================================================================!

    implicit none

    integer           ::  i,ndim,nat

    real(kind=dprec)     ::  r1(ndim),r2(ndim)
    real(kind=dprec)     ::  x1(nat),y1(nat),z1(nat)
    real(kind=dprec)     ::  x2(nat),y2(nat),z2(nat)
    real(kind=dprec)     ::  xcent,ycent,zcent
    real(kind=dprec)     ::  rmsval
    real(kind=dprec)     ::  xnew,ynew,znew
    real(kind=dprec)     ::  xxyx,xxyy,xxyz,xyyx,xyyy,xyyz,xzyx,xzyy,xzyz
    real(kind=dprec)     ::  umat(3,3),work1(4),work2(4)
    real(kind=dprec)     ::  q(4),d(4),c(4,4),v(4,4)

    !fill up some temporary arrays
    do i = 1,nat
      x1(i) = r1(3*i-2)
      y1(i) = r1(3*i-1)
      z1(i) = r1(3*i)
      x2(i) = r2(3*i-2)
      y2(i) = r2(3*i-1)
      z2(i) = r2(3*i)
    end do

    !centroid of second structure
    xcent = 0.0_dprec
    ycent = 0.0_dprec
    zcent = 0.0_dprec
    do i = 1,nat
        xcent = xcent + x2(i)
        ycent = ycent + y2(i)
        zcent = zcent + z2(i)
    end do
    xcent = xcent / nat
    ycent = ycent / nat
    zcent = zcent / nat

    !shift structure 2 to origin
    do i = 1,nat
        x2(i) = x2(i) - xcent
        y2(i) = y2(i) - ycent
        z2(i) = z2(i) - zcent
    end do

    !centroid of first structure
    xcent = 0.0_dprec
    ycent = 0.0_dprec
    zcent = 0.0_dprec
    do i = 1,nat
        xcent = xcent + x1(i)
        ycent = ycent + y1(i)
        zcent = zcent + z1(i)
    end do

    !preserve centroid of first structure as reference
    xcent = xcent / nat
    ycent = ycent / nat
    zcent = zcent / nat

    !shift structure 1 to origin
    do i = 1,nat
        x1(i) = x1(i) - xcent
        y1(i) = y1(i) - ycent
        z1(i) = z1(i) - zcent
    end do

    !calculate optimum quaternion q
    !upper triangle of quadratic matrix
    xxyx = 0.0_dprec
    xxyy = 0.0_dprec
    xxyz = 0.0_dprec
    xyyx = 0.0_dprec
    xyyy = 0.0_dprec
    xyyz = 0.0_dprec
    xzyx = 0.0_dprec
    xzyy = 0.0_dprec
    xzyz = 0.0_dprec
    do i = 1,nat
       xxyx = xxyx + x1(i)*x2(i)
       xxyy = xxyy + y1(i)*x2(i)
       xxyz = xxyz + z1(i)*x2(i)
       xyyx = xyyx + x1(i)*y2(i)
       xyyy = xyyy + y1(i)*y2(i)
       xyyz = xyyz + z1(i)*y2(i)
       xzyx = xzyx + x1(i)*z2(i)
       xzyy = xzyy + y1(i)*z2(i)
       xzyz = xzyz + z1(i)*z2(i)
    end do

    !setup c
    c(1,1) = xxyx + xyyy + xzyz
    c(1,2) = xzyy - xyyz
    c(2,2) = xxyx - xyyy - xzyz
    c(1,3) = xxyz - xzyx
    c(2,3) = xxyy + xyyx
    c(3,3) = xyyy - xzyz - xxyx
    c(1,4) = xyyx - xxyy
    c(2,4) = xzyx + xxyz
    c(3,4) = xyyz + xzyy
    c(4,4) = xzyz - xxyx - xyyy

    !diagonalize c
    call tssearch_algor_jacobi(4,4,c,d,v,work1,work2)

    !extract the desired quaternion
    q(1) = v(1,4)
    q(2) = v(2,4)
    q(3) = v(3,4)
    q(4) = v(4,4)

    !create the rotation matrix
    umat(1,1) = q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2
    umat(2,1) = 2.0_dprec * (q(2) * q(3) - q(1) * q(4))
    umat(3,1) = 2.0_dprec * (q(2) * q(4) + q(1) * q(3))
    umat(1,2) = 2.0_dprec * (q(3) * q(2) + q(1) * q(4))
    umat(2,2) = q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2
    umat(3,2) = 2.0_dprec * (q(3) * q(4) - q(1) * q(2))
    umat(1,3) = 2.0_dprec * (q(4) * q(2) - q(1) * q(3))
    umat(2,3) = 2.0_dprec * (q(4) * q(3) + q(1) * q(2))
    umat(3,3) = q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2

    !rotate structure 2 to fit structure 1
    do i = 1,nat
       xnew = x2(i)*umat(1,1) + y2(i)*umat(1,2) + z2(i)*umat(1,3)
       ynew = x2(i)*umat(2,1) + y2(i)*umat(2,2) + z2(i)*umat(2,3)
       znew = x2(i)*umat(3,1) + y2(i)*umat(3,2) + z2(i)*umat(3,3)
       x2(i) = xnew
       y2(i) = ynew
       z2(i) = znew
    end do

    !translate first structure back to original position
    !refer structure 2 back to structure 1
    do i = 1,nat
       x1(i) = x1(i) + xcent
       y1(i) = y1(i) + ycent
       z1(i) = z1(i) + zcent
    end do
    do i = 1,nat
       x2(i) = x2(i) + xcent
       y2(i) = y2(i) + ycent
       z2(i) = z2(i) + zcent
    end do

    !return updated coordinates
    do i = 1,nat
      r1(3*i-2) = x1(i)
      r1(3*i-1) = y1(i)
      r1(3*i)   = z1(i)
      r2(3*i-2) = x2(i)
      r2(3*i-1) = y2(i)
      r2(3*i)   = z2(i)
    end do

    !calculate rmsvalue between structures after transformation
    call tssearch_algor_rmscalc(nat,ndim,r1,r2,rmsval)

    return
    end subroutine tssearch_algor_supmol

    subroutine tssearch_algor_rmscalc(nat,ndim,x1,x2,rmsval)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Sept 10 2001                                   !
    !=========================================================================!

    implicit none

    integer       :: i,ndim,nat

    real(kind=dprec) :: rmsfit,rmsterm,rmsval
    real(kind=dprec) :: xr,yr,zr,dist2
    real(kind=dprec) :: weight,norm
    real(kind=dprec) :: x1(ndim),x2(ndim)

    !compute the rms deviation of structures
    rmsfit = 0.0_dprec
    norm = 0.0_dprec
    do i = 1,nat
        weight = 1.0_dprec
        xr = x1(3*i-2) - x2(3*i-2)
        yr = x1(3*i-1) - x2(3*i-1)
        zr = x1(3*i)   - x2(3*i)
        dist2 = xr**2 + yr**2 + zr**2
        norm = norm + weight
        rmsterm = dist2 * weight
        rmsfit = rmsfit + rmsterm
    end do
    rmsval = sqrt(rmsfit/norm)

    return
    end subroutine tssearch_algor_rmscalc

    subroutine tssearch_algor_optipnt(nat,ndim,period,f,xx,x0,x1,pm,xm,iopt,idfix)
    !=========================================================================!
    ! Calculates optimized linear and quadratic synchronized coordinates      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  nat     , input  ! number of atoms                                     !
    !  ndim    , input  ! number of raw degrees of freedom 3*nat              !
    !  period  , input  ! periodic or non-periodic model                      !
    !  f       , input  ! interpolation factor                                !
    !  xx      , output ! optimized interpolated point                        !
    !  x0      , input  ! endpoint 1                                          !
    !  x1      , input  ! endpoint 2                                          !
    !  pm      , input  ! intermediate path value based on rms deviations     !
    !  xm      , input  ! intermediate point for quadratic interpolation      !
    !  iopt    , input  ! lst (1) mode or qst (2) mode                        !
    !  idfix   , input  ! fixed atom information                              !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Aug 26 2000                                    !
    !=========================================================================!

    implicit none

    ! Arguments

    integer          :: nat,ndim,iopt,idfix(ndim)
    logical          :: period
    real(kind=dprec)    :: f,xx(ndim),x0(ndim),x1(ndim),xm(ndim),pm

    ! Local variables

    integer          :: i
    real(kind=dprec)    :: value,grdmin

    !Initialize optimization parameters
    grdmin = 1.0e-10_dprec

    !Interpolate the positions
    do i = 1,ndim
         xx(i) = (1.0_dprec-f)*x0(i) + f*x1(i)
    end do

    !Optimize the synchronous transit function
    call tssearch_algor_varmet(nat,ndim,period,xx,x0,x1,idfix,f,pm,xm, &
             iopt,value,grdmin)

    !Need to design a more generic routine to replace this
!   call tssearch_utils_safe(ndim,nat,period,x0,x1,xx,s,idfix,0.0_dprec,0,2)

    return
    end subroutine tssearch_algor_optipnt

    subroutine tssearch_algor_tgthess(nat,ndim,idfix,iopt,xx,x0,x1,f,g,xm,pm, &
                         period,g_rms,tang,g_tang,dgdf,dedf,neg)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  nat     , input  ! number of atoms                                     !
    !  ndim    , input  ! number of raw degrees of freedom 3*nat              !
    !  period  , input  ! periodic or non-periodic model                      !
    !  f       , input  ! interpolation factor                                !
    !  xx      , output ! optimized interpolated point                        !
    !  x0      , input  ! endpoint 1                                          !
    !  x1      , input  ! endpoint 2                                          !
    !  pm      , input  ! intermediate path value based on rms deviations     !
    !  xm      , input  ! intermediate point for quadratic interpolation      !
    !  iopt    , input  ! lst (1) mode or qst (2) mode                        !
    !  idfix   , input  ! fixed atom information                              !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Sept 09 2000                                   !
    !=========================================================================!
    use tssearch_utils, only: tssearch_utils_energy_gradient

    implicit none

    integer       ::  i,ndim,nat,idfix(ndim),iopt,neg

    real(kind=dprec) ::  g_rms,g_tang,f,tang_norm
    real(kind=dprec) ::  g2,f0,ff,fb,dedf
    real(kind=dprec) ::  delta,energyf,energyb
    real(kind=dprec) ::  x0(ndim),xf(ndim),xb(ndim),x1(ndim)
    real(kind=dprec) ::  g(ndim),gf(ndim),gb(ndim)
    real(kind=dprec) ::  tang(ndim),dgdf(ndim),pm,xm(ndim)
    real(kind=dprec) ::  xx(ndim)

    logical       ::  period

    !parameters
    delta = 0.01_dprec
    f0 = f

    !gradient norm
    g2 = 0.0_dprec
    do i = 1,ndim
       if (idfix(i).eq.1) then
          g(i) = 0.0_dprec
       end if
       g2 = g2 + g(i)**2
    end do
    g_rms = sqrt(g2/dble(nat))

    !perturb forward
    do i = 1,ndim
       xf(i) = xx(i)
    end do
    ff = f0 + delta

    !write(ts_stdout,*)
    !write(ts_stdout,'(a,f18.4)')     ' Forward Step : ',ff

    !optimize point
    call tssearch_algor_optipnt(nat,ndim,period,ff,xf,x0,x1,pm,xm,iopt,idfix)

    !calculate the energy and gradient
    call tssearch_utils_energy_gradient(ndim,xf,gf,energyf)

    !increment energy-gradient counter
    neg = neg + 1

    !perturb backward
    do i = 1,ndim
       xb(i) = xx(i)
    end do
    fb = f0 - delta

    !write(ts_stdout,*)
    !write(ts_stdout,'(a,f18.4)')     ' Backward Step : ',fb

    !optimize point
    call tssearch_algor_optipnt(nat,ndim,period,fb,xb,x0,x1,pm,xm,iopt,idfix)

    !calculate the energy and gradient
    call tssearch_utils_energy_gradient(ndim,xb,gb,energyb)

    !increment energy-gradient counter
    neg = neg + 1

    !reinstate f
    f = f0

    !path tangent
    tang_norm = 0.0_dprec
    do i = 1,ndim
       if (idfix(i).eq.1) then
           gf(i) = 0.0_dprec
           gb(i) = 0.0_dprec
       end if
       tang(i) = xf(i) - xb(i)
       tang_norm = tang_norm + tang(i)*tang(i)
       dgdf(i) = gf(i) - gb(i)
    end do
    tang_norm = dsqrt(tang_norm)

    !projected gradient
    g_tang = 0.0_dprec
    do i = 1,ndim
       tang(i) = tang(i) / tang_norm     !normalize tangent
       g_tang = g_tang + g(i)*tang(i)    !tangent projection
    end do
    g_tang = g_tang / sqrt(dble(nat))

    dedf = (energyf - energyb)/2._dprec/delta

    return
    end subroutine tssearch_algor_tgthess

    subroutine tssearch_algor_periodic_adjust(xr,yr,zr,i)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Sept 09 2000                                   !
    !=========================================================================!
    use tssearch_utils, only: icell, gamma_cos, beta_term, gamma_term, &
         gamma_sin, beta_cos, alength, blength, clength

    implicit none

    ! Arguments
    integer, intent(in)           :: i
    real(kind=dprec), intent(inout)  :: xr,yr,zr

    ! Local variables
    real(kind=dprec)  :: xmove,ymove,zmove
    real(kind=dprec)  :: xfrac,yfrac,zfrac
    real(kind=dprec)  :: xbox,ybox,zbox
    
    !setup some local variables
    xbox  = alength
    ybox  = blength
    zbox  = clength

    !shift reference to images if needed
    if (i.eq.0) then
       xmove = 0.0_dprec
       ymove = 0.0_dprec
       zmove = 0.0_dprec
    else
       xmove = icell(1,i) * xbox
       ymove = icell(2,i) * ybox
       zmove = icell(3,i) * zbox
    end if

    !transform to fractional
    zfrac = zr / gamma_term
    yfrac = (yr - zfrac*beta_term) / gamma_sin
    xfrac = xr - yfrac*gamma_cos - zfrac*beta_cos

    !shift cell
    xfrac = xfrac + xmove
    yfrac = yfrac + ymove
    zfrac = zfrac + zmove

    !revert back to cartesian
    xr = xfrac + yfrac*gamma_cos + zfrac*beta_cos
    yr = yfrac*gamma_sin + zfrac*beta_term
    zr = zfrac * gamma_term

    return
    end subroutine tssearch_algor_periodic_adjust

    subroutine tssearch_algor_varmet(nat,ndim,period,x0,xmin1,xmin2,&
                  idfix,fl,pm,xm,iopt,f0,grdmin)
    !=========================================================================!
    ! Davidon's variable metric algorithm                                     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, Nov 30 2001                               !
    !=========================================================================!
    use tssearch_utils, only: ts_stdout, tssearch_utils_io_abort

    implicit none

    ! Arguments
    logical        :: period
    integer        :: ndim,nat,idfix(ndim),iopt
    real(kind=dprec)  :: x0(ndim),xmin1(ndim),xmin2(ndim),pm,xm(ndim),fl,f0
    real(kind=dprec)  :: grdmin

    ! Local variables
    integer        :: i,j,mdim,ncalls
    integer        :: niter,nbigang,maxbigang
    integer        :: maxiter,nextiter,ierr

    real(kind=dprec)  :: epsilon
    real(kind=dprec)  :: f,fprime,f0old,f0prime,srchnorm
    real(kind=dprec)  :: sgangle,sg,snorm,lambda,cosang
    real(kind=dprec)  :: fmove,xmove,gnorm,grms,rms
    real(kind=dprec)  :: m2,n2,u2,v,mu,mw,us,qk0
    real(kind=dprec)  :: a,b,b0,c,alpha_varmet,gamma_varmet,delta_varmet
    real(kind=dprec)  :: x(ndim),x0old(ndim)
    real(kind=dprec)  :: g(ndim),search(ndim)
    real(kind=dprec)  :: s(ndim),w(ndim)
    real(kind=dprec)  :: k(ndim),k0(ndim)
    real(kind=dprec)  :: m(ndim),n(ndim),u(ndim)
    real(kind=dprec)  :: p(ndim),q(ndim),jq(ndim)
    real(kind=dprec)  :: scale(ndim),fctmin
    real(kind=dprec)  :: stpmax,jmatrix_guess,angmax

    !save jmatrix for next entry
!   real(kind=dprec), dimension(:,:), allocatable, save  :: jmatrix
    real(kind=dprec)  :: jmatrix(ndim,ndim)

    logical        :: restart,converged,set_scale

    !initialization and set-up for the optimization
    if (ndim .gt. 3*nat) then
         write(ts_stdout,*) 'Error: Number of dimensions greater than 3 times number of atoms'
         ierr = 1
         if (ierr/=0) call tssearch_utils_io_abort('Error detected in tssearch_algor_varmet')
    end if

    !initialize parameters
    mdim = ndim
    rms = sqrt(dble(ndim))
    rms = rms / sqrt(3.0_dprec)
    maxbigang = 2
    epsilon = 1.0e-08_dprec
    restart = .true.
    converged = .false.
    fctmin = -1000000.0_dprec
    maxiter = 1000
    nextiter = 1
    stpmax = stpmax_varmet
    jmatrix_guess = 0.4_dprec
    angmax = 180.0_dprec
    nbigang = 0

    !initialize scale
    set_scale = .false.
    do i = 1, ndim
      scale(i) = 1.0_dprec
    end do
    if (.not. set_scale) then
     do i = 1, ndim
       if (scale(i) .eq. 0.0_dprec)  scale(i) = 1.0_dprec
     end do
    end if

!    !allocate jmatrix
!    allocate(jmatrix(ndim,ndim),stat=ierr)
!    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating jmatrix in tssearch_algor_varmet')
 
    !get initial function and gradient values
    niter = nextiter - 1
    maxiter = niter + maxiter
    do i = 1,ndim
      x0old(i) = x0(i)
    end do

    !call the energy-gradient calculator
    f0 = tssearch_algor_transit(ndim,period,idfix,x0,xmin1,xmin2,fl,pm,xm,iopt,g)

    !initialize function-gradient call counter
    ncalls = 1

    !zero out gradients on fixed atoms
    do i = 1, ndim
      if (idfix(i).eq.1) then
          g(i) = 0.0_dprec
      end if
    end do

    !store initial function value
    f0old = f0

    !start: setup the initial jmatrix
    !outer loop
    do while (.not. converged)

     if (restart) then

          !jmatrix is initialized to a diagonal matrix
          do j = 1, mdim
             do i = 1, ndim
                jmatrix(i,j) = 0.0_dprec
             end do
          end do

          !fill diagonal elements
          do j = 1, mdim
             jmatrix(j,j) = jmatrix_guess
          end do

          !calculate k0 and w
          do j = 1, mdim
             k0(j) = 0.0_dprec
             do i = 1, ndim
                k0(j) = k0(j) + jmatrix(i,j)*g(i)
             end do
             w(j) = k0(j)
          end do

          restart = .false.
     end if

     !calculate rms gradient, rms change, change in energy
     gnorm = 0.0_dprec
     grms = 0.0_dprec
     do i = 1, ndim
        gnorm = gnorm + g(i)**2
        grms = grms + (g(i)*scale(i))**2
     end do
     gnorm = sqrt(gnorm)
     grms = sqrt(grms) / rms
     xmove = 0.0_dprec

     if (niter .ne. 0) then
        do i = 1, ndim
           xmove = xmove + ((x0(i)-x0old(i))/scale(i))**2
           x0old(i) = x0(i)
        end do
        xmove = sqrt(xmove) / rms
        fmove = f0old - f0
        f0old = f0
     end if

     !check norm of gradient, change in function, iteration number
     !stopping condition
     if (grms.lt.grdmin .or. f0.lt.fctmin .or. niter.ge.maxiter) then
          converged = .true.
          goto 700
     end if

     !step 1: commence next iteration
     niter = niter + 1
     sg = 0.0_dprec
     snorm = 0.0_dprec
     do j = 1, mdim
        s(j) = -k0(j)
        snorm = snorm + s(j)**2
        sg = sg - s(j)*g(j)
     end do
     f0prime = -snorm
     snorm = sqrt(snorm)
     cosang = sg / (snorm*gnorm)
     cosang = min(1.0_dprec,max(-1.0_dprec,cosang))
     sgangle = acos(cosang) * 180.0_dprec/pi
     if (sgangle .gt. angmax) then
        nbigang = nbigang + 1
     else
        nbigang = 0
     end if
     lambda = 2.0_dprec

     !comments after step 1 in davidon's paper
     if (4.0_dprec*(f0-fctmin) .lt. -f0prime) then
         do j = 1, mdim
             s(j) = -s(j) * (4.0_dprec*(f0-fctmin)/f0prime)
         end do
         f0prime = -4.0_dprec * (f0-fctmin)
     end if

     !step 2: find the next starting point

 200  continue

     !initialize starting point
     do i = 1, ndim
          search(i) = 0.0_dprec
     end do

     !calculate x = x0 + Js  -> calculate Js
     do j = 1, mdim
       do i = 1, ndim
           search(i) = search(i) + jmatrix(i,j)*s(j)
       end do
     end do

     srchnorm = 0.0_dprec
     do i = 1, ndim
        srchnorm = srchnorm + search(i)**2
     end do
     srchnorm = sqrt(srchnorm)

     if (srchnorm .gt. stpmax) then
        do j = 1, mdim
           s(j) = (stpmax/srchnorm) * s(j)
        end do
        do i = 1, ndim
           search(i) = (stpmax/srchnorm) * search(i)
        end do
        f0prime = (stpmax/srchnorm) * f0prime
        lambda = 0.5_dprec
     end if

     !check if -f0prime < epsilon
     !if true, stop
     if (-f0prime .lt. epsilon) then
        converged = .true.
        goto 700
     end if

     !check fixed atom information
     !update position: x = x0 + Js
     do i = 1, ndim
        if (idfix(i).eq.1) then
           x(i) = x0(i)
        else
           x(i) = x0(i) + search(i)
        end if
     end do

     !call function-gradient calculator
     f = tssearch_algor_transit(ndim,period,idfix,x,xmin1,xmin2,fl,pm,xm,iopt,g)

     !update counter
     ncalls = ncalls + 1

     !zero out gradients on fixed atoms
     do i = 1, ndim
        if (idfix(i).eq.1) then
           g(i) = 0.0_dprec
        end if
     end do

     !is f > f0 ?
     !if true, repeat step 2. otherwise jump to step 3
     if (f .ge. f0) then
         do j = 1, mdim
            s(j) = 0.5_dprec * s(j)
         end do
         f0prime = 0.5_dprec * f0prime
         lambda = 0.5_dprec
          goto 200        !step 2
     end if

     !step 3: check if we need to update the step or take another one
     !calculate k
     do j = 1, mdim
        k(j) = 0.0_dprec
        do i = 1, ndim
           k(j) = k(j) + jmatrix(i,j)*g(i)
        end do
     end do

     !set fprime and b0 
     fprime = 0.0_dprec
     do j = 1, mdim
        fprime = fprime + k(j)*s(j)
     end do
     b0 = fprime - f0prime

     !calculate m, update quantities as per step 3 in paper
     do j = 1, mdim
        m(j) = s(j) + k0(j) - k(j)
        k0(j) = k(j)
     end do
     do i = 1, ndim
        x0(i) = x(i)
     end do
     f0 = f
     f0prime = fprime

     !is b0 < epsilon ? if not, jump to step 4
     !if true, update s and return to step 2
     if (b0 .lt. epsilon) then
          do j = 1, mdim
             s(j) = s(j) * lambda
          end do
          f0prime = f0prime * lambda
          goto 200      !step 2
     end if

     !step 4:  do we need to update
     if (nbigang .ge. maxbigang) then
        restart = .true.
        goto 700
     end if

     m2 = 0.0_dprec
     do j = 1, mdim
        m2 = m2 + m(j)**2
     end do

     !is m2 < epsilon. is yes, return to step 1 (start)
     if (m2 .lt. epsilon) then
        goto 700
     end if

     v = 0.0_dprec
     do j = 1, mdim
        v = v + m(j)*s(j)
     end do

     !calculate mu
     mu = v - m2

     !calculate u
     mw = 0.0_dprec
     do j = 1, mdim
        mw = mw + m(j)*w(j)
     end do
     do j = 1, mdim
        u(j) = w(j) - m(j)*(mw/m2)
     end do
     u2 = 0.0_dprec
     do j = 1, mdim
        u2 = u2 + u(j)**2
     end do

     !check condition 
     if (m2*u2 .ge. epsilon) then
 
     !step 4a
          us = 0.0_dprec
          do j = 1, mdim
             us = us + u(j)*s(j)
          end do
          do j = 1, mdim
             n(j) = u(j)*(us/u2)
          end do
          n2 = us * us/u2
     else
          do j = 1, mdim
             n(j) = 0.0_dprec
          end do
          n2 = 0.0_dprec
     end if

     !step 5:  test inner product of projected s and del-g
     b = n2 + mu * v/m2
     if (b .lt. epsilon) then
         do j = 1, mdim
            n(j) = s(j) - m(j)*(v/m2)
         end do
         n2 = b0 - mu * v/m2
         b = b0
     end if

     !step 6: calculate gamma_varmet and delta_varmet
     if (mu*v .ge. m2*n2) then
         gamma_varmet = 0.0_dprec
         delta_varmet = sqrt(v/mu)
     else
     !step 6a
         a = b - mu
         c = b + v
         gamma_varmet = sqrt((1.0_dprec-mu*v/(m2*n2))/(a*b))
         delta_varmet = sqrt(c/a)
         if (c .lt. a) then
            gamma_varmet = -gamma_varmet
         end if
     end if

     !step 7: update jmatrix
     alpha_varmet = v + mu*delta_varmet + m2*n2*gamma_varmet
     do j = 1, mdim
      p(j) = m(j)*(delta_varmet-n2*gamma_varmet) + n(j)*(gamma_varmet*v)
      q(j) = m(j)*((1.0_dprec+n2*gamma_varmet)/alpha_varmet) - n(j)*(gamma_varmet * mu/alpha_varmet)
      w(j) = m(j)*(n2*(1.0_dprec+gamma_varmet*mu*v)/alpha_varmet) - &
                     n(j)*((1.0_dprec+delta_varmet)*mu*v/alpha_varmet)
     end do

     ! k0 = k0 + p(qT)k0         !qT = q transpose
     qk0 = 0.0_dprec
     do j = 1, mdim
        qk0 = qk0 + q(j)*k0(j)
     end do
     do j = 1, mdim
        k0(j) = k0(j) + p(j)*qk0
     end do

     ! J = J + Jq(pT)            !pT = p transpose
     do i = 1, ndim
        jq(i) = 0.0_dprec
     end do
     do j = 1, mdim
        do i = 1, ndim
           jq(i) = jq(i) + jmatrix(i,j)*q(j)
        end do
     end do
     do j = 1, mdim
          do i = 1, ndim
             jmatrix(i,j) = jmatrix(i,j) + jq(i)*p(j)
          end do
     end do
     if (n2 .le. 0.0_dprec) then
        do j = 1, mdim
           w(j) = k0(j)
        end do
     end if

 700  continue

    end do  !outer do while

    return
    end subroutine tssearch_algor_varmet

    function tssearch_algor_transit(ndim,period,idfix,xx,xmin1,xmin2,t,pm,xm,iopt,g)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Sept 15 2000                                   !
    !=========================================================================!
    use tssearch_utils, only: ncell, other_images

    implicit none

    logical           :: period
    integer           :: i,j,ix,iy,iz,jx,jy,jz,ndim,n,k,iopt

    integer           :: idfix(ndim)
    real(kind=dprec)     :: value,xx(ndim),g(ndim),t
    real(kind=dprec)     :: xci,yci,zci,xcd,ycd,zcd,xmin1(ndim)
    real(kind=dprec)     :: x1i,y1i,z1i,x1d,y1d,z1d,xmin2(ndim)
    real(kind=dprec)     :: x2i,y2i,z2i,x2d,y2d,z2d,xm(ndim),pm
    real(kind=dprec)     :: xmi,ymi,zmi,xmd,ymd,zmd
    real(kind=dprec)     :: term,termx,termy,termz,gamma_hg
    real(kind=dprec)     :: r1,r2,rc,rm,ri,ri4,rd,wi,wc,wd
    real(kind=dprec)     :: tq,pq,cutoff,cutoff2
    
    real(kind=dprec)     :: tssearch_algor_transit
      
    integer           :: mode

    !zero out the synchronous transit function and gradients
    value = 0.0_dprec
    n=ndim/3
    do i = 1,ndim
      g(i) = 0.0_dprec
    end do
    tq = 1.0_dprec - t

    xmi = xm(1) ! just to satisfy FTNCHECK
    ymi = xm(2)
    zmi = xm(3)

    !set the cutoff distance for interatomic distances
    cutoff = 1000.0_dprec
    cutoff2 = cutoff**2

    !set the type of synchronous transit path to be used
    if (iopt .eq. 1) then
        mode = 1                      !linear (LST)
        pq=0 ! this variable has to be initialised for iopt=1
    else
        mode = 2                      !quadratic (QST)
        pq = 1.0_dprec - pm
    end if

    !check if periodic case or not
    if (period) then

    !periodic part

    !part involving primary cell
      do i = 1, n-1
         iz = 3 * i
         iy = iz - 1
         ix = iz - 2
         xci = xx(ix)
         yci = xx(iy)
         zci = xx(iz)
         x1i = xmin1(ix)
         y1i = xmin1(iy)
         z1i = xmin1(iz)
         x2i = xmin2(ix)
         y2i = xmin2(iy)
         z2i = xmin2(iz)
         if (mode .eq. 2) then               ! QST
            xmi = xm(ix)
            ymi = xm(iy)
            zmi = xm(iz)
         end if
         do j = i+1, n
            jz = 3 * j
            jy = jz - 1
            jx = jz - 2
            xcd = xci - xx(jx)
            ycd = yci - xx(jy)
            zcd = zci - xx(jz)

            x1d = x1i - xmin1(jx)
            y1d = y1i - xmin1(jy)
            z1d = z1i - xmin1(jz)

            x2d = x2i - xmin2(jx)
            y2d = y2i - xmin2(jy)
            z2d = z2i - xmin2(jz)

            rc = xcd**2 + ycd**2 + zcd**2
            r1 = x1d**2 + y1d**2 + z1d**2
            r2 = x2d**2 + y2d**2 + z2d**2

            if (min(rc,r1,r2) .lt. cutoff2) then
               rc = sqrt(rc)
               r1 = sqrt(r1)
               r2 = sqrt(r2)
               ri = tq*r1 + t*r2
                if (mode .eq. 2) then          ! QST
                  xmd = xmi - xm(jx)
                  ymd = ymi - xm(jy)
                  zmd = zmi - xm(jz)

                  rm = sqrt(xmd**2+ymd**2+zmd**2)
                  gamma_hg = (rm-pq*r1-pm*r2) / (pm*pq)
                  ri = ri + gamma_hg*t*tq
               end if
               ri4 = ri**4
               rd = rc - ri
               value = value + (rd**2)*tssearch_algor_wgtfunc(rcut,ri)
               term = 2.0_dprec * (rd/rc)*tssearch_algor_wgtfunc(rcut,ri)
               termx = term * xcd
               termy = term * ycd
               termz = term * zcd
               g(ix) = g(ix) + termx
               g(iy) = g(iy) + termy
               g(iz) = g(iz) + termz
               g(jx) = g(jx) - termx
               g(jy) = g(jy) - termy
               g(jz) = g(jz) - termz
            end if
         end do
      end do

      !part involving image cells
      if (other_images) then

      do i = 1, n
         iz = 3 * i
         iy = iz - 1
         ix = iz - 2
         xci = xx(ix)
         yci = xx(iy)
         zci = xx(iz)
         x1i = xmin1(ix)
         y1i = xmin1(iy)
         z1i = xmin1(iz)
         x2i = xmin2(ix)
         y2i = xmin2(iy)
         z2i = xmin2(iz)
         if (mode .eq. 2) then             ! QST
            xmi = xm(ix)
            ymi = xm(iy)
            zmi = xm(iz)
         end if
         do j = i, n
            jz = 3 * j
            jy = jz - 1
            jx = jz - 2

            do k=1,ncell

               xcd = xci - xx(jx)
               ycd = yci - xx(jy)
               zcd = zci - xx(jz)

               call tssearch_algor_periodic_adjust(xcd,ycd,zcd,k)

               x1d = x1i - xmin1(jx)
               y1d = y1i - xmin1(jy)
               z1d = z1i - xmin1(jz)

               call tssearch_algor_periodic_adjust(x1d,y1d,z1d,k)

               x2d = x2i - xmin2(jx)
               y2d = y2i - xmin2(jy)
               z2d = z2i - xmin2(jz)

               call tssearch_algor_periodic_adjust(x2d,y2d,z2d,k)

               rc = xcd**2 + ycd**2 + zcd**2
               r1 = x1d**2 + y1d**2 + z1d**2
               r2 = x2d**2 + y2d**2 + z2d**2

               if (min(rc,r1,r2) .lt. cutoff2) then
                 rc = sqrt(rc)
                 r1 = sqrt(r1)
                 r2 = sqrt(r2)
                 ri = tq*r1 + t*r2

                 if (mode .eq. 2) then           ! QST
                   xmd = xmi - xm(jx)
                   ymd = ymi - xm(jy)
                   zmd = zmi - xm(jz)

                   call tssearch_algor_periodic_adjust(xmd,ymd,zmd,k)

                   rm = sqrt(xmd**2+ymd**2+zmd**2)
                   gamma_hg = (rm-pq*r1-pm*r2) / (pm*pq)
                   ri = ri + gamma_hg*t*tq
                 end if
                 ri4 = ri**4
                 rd = rc - ri
                 value = value + (rd**2)*tssearch_algor_wgtfunc(rcut,ri)
                 term = 2.0_dprec * (rd/rc)*tssearch_algor_wgtfunc(rcut,ri)
                 termx = term * xcd
                 termy = term * ycd
                 termz = term * zcd
                 g(ix) = g(ix) + termx
                 g(iy) = g(iy) + termy
                 g(iz) = g(iz) + termz
                 g(jx) = g(jx) - termx
                 g(jy) = g(jy) - termy
                 g(jz) = g(jz) - termz
              end if
            end do
         end do
      end do
      end if

      !break degeneracies
      do i = 1, ndim

         if (idfix(i).eq.1) then
           g(i) = 0.0_dprec
         end if

         wc = xx(i)
         wi = tq*xmin1(i) + t*xmin2(i)
         wd = wc - wi
         value = value + 0.000001_dprec*wd**2
         g(i) = g(i) + 0.000002_dprec*wd
      end do

      tssearch_algor_transit = value


      else                          !periodic cell check


      !molecular part

      !portion based on interpolated interatomic distances

      do i = 1, n-1
         iz = 3 * i
         iy = iz - 1
         ix = iz - 2
         xci = xx(ix)
         yci = xx(iy)
         zci = xx(iz)
         x1i = xmin1(ix)
         y1i = xmin1(iy)
         z1i = xmin1(iz)
         x2i = xmin2(ix)
         y2i = xmin2(iy)
         z2i = xmin2(iz)
         if (mode .eq. 2) then         ! QST
            xmi = xm(ix)
            ymi = xm(iy)
            zmi = xm(iz)
         end if
         do j = i+1, n
            jz = 3 * j
            jy = jz - 1
            jx = jz - 2
            xcd = xci - xx(jx)
            ycd = yci - xx(jy)
            zcd = zci - xx(jz)
            x1d = x1i - xmin1(jx)
            y1d = y1i - xmin1(jy)
            z1d = z1i - xmin1(jz)
            x2d = x2i - xmin2(jx)
            y2d = y2i - xmin2(jy)
            z2d = z2i - xmin2(jz)
            rc = xcd**2 + ycd**2 + zcd**2
            r1 = x1d**2 + y1d**2 + z1d**2
            r2 = x2d**2 + y2d**2 + z2d**2
            if (min(rc,r1,r2) .lt. cutoff2) then
              rc = sqrt(rc)
              r1 = sqrt(r1)
              r2 = sqrt(r2)
              ri = tq*r1 + t*r2
               if (mode .eq. 2) then          ! QST
                  xmd = xmi - xm(jx)
                  ymd = ymi - xm(jy)
                  zmd = zmi - xm(jz)
                  rm = sqrt(xmd**2+ymd**2+zmd**2)
                  gamma_hg = (rm-pq*r1-pm*r2) / (pm*pq)
                  ri = ri + gamma_hg*t*tq
               end if
               ri4 = ri**4
               rd = rc - ri
               value = value + rd**2/ri4
               term = 2.0_dprec * rd/(ri4*rc)
               termx = term * xcd
               termy = term * ycd
               termz = term * zcd
               g(ix) = g(ix) + termx
               g(iy) = g(iy) + termy
               g(iz) = g(iz) + termz
               g(jx) = g(jx) - termx
               g(jy) = g(jy) - termy
               g(jz) = g(jz) - termz
            end if
         end do
      end do

      !break degeneracies
      do i = 1, ndim

         if (idfix(i).eq.1) then
           g(i) = 0.0_dprec
         end if

         wc = xx(i)
         wi = tq*xmin1(i) + t*xmin2(i)
         wd = wc - wi
         value = value + 0.000001_dprec*wd**2
         g(i) = g(i) + 0.000002_dprec*wd
      end do

      tssearch_algor_transit = value

      end if                          !period

      return
    end function tssearch_algor_transit

    subroutine tssearch_algor_lstqstbracket(nat,ndim,idfix,xmin1,xmin2,f0,aa,bb,cc, &
                 iopt, xqst, pqst, flmax, tol, stepinc, step0, &
                 faa, fbb, fcc, numextr, period, gdum, ifirst, neg, &
                 xclean1, xclean2, output_file)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! tries to find a bracket for a nearest to x max or min
    ! along the line  x+fl*p, fl belongs to [0.0,flmax],
    ! 
    ! INPUT:
    ! ndim - dimension of x,p
    ! x(ndim) - initial point for line search
    ! p(ndim) - direction of line search
    ! f0 - function value at x
    ! flmax - upper bound for a,b,c
    ! tol - tolerance for a,b,c
    ! step0 - initial step
    ! stepinc - the factor of increasing step during bracket search
    ! numextr = 1  - bracketing a min
    !         =-1  - bracketing a max
    ! ifirst  =0/1 first or next call to the routine
    ! xqst,pqst
    !
    ! OUTPUT:
    ! aa,bb,cc - bracketing triplet
    ! faa,fbb,fcc - corresponding function values
    ! numextr = 0  - bracketing failed
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Sept 20 2000                                   !
    !=========================================================================!
    use tssearch_utils, only: root_write, scf_fail, ts_stdout, &
         tssearch_utils_energy_gradient, tssearch_utils_write, &
         tssearch_utils_write_summary, tssearch_utils_io_abort

    implicit none

    logical          :: period

    integer          :: nat,ndim,neg,ii
    integer          :: idfix(ndim)
    integer          :: ifirst,ifailed,numextr
    integer          :: iopt

    real(kind=dprec)    :: gdum(ndim)
    real(kind=dprec)    :: xx(ndim)
    real(kind=dprec)    :: xmin1(ndim)
    real(kind=dprec)    :: xmin2(ndim)
    real(kind=dprec)    :: xclean1(ndim)
    real(kind=dprec)    :: xclean2(ndim)
    real(kind=dprec)    :: xqst(ndim)
    real(kind=dprec)    :: pqst
    real(kind=dprec)    :: flmax,step0,tol
    real(kind=dprec)    :: step_sk
    real(kind=dprec)    :: a,b,fa,f0
    real(kind=dprec)    :: aa,bb,faa,fbb
    real(kind=dprec)    :: step,stepinc,c
    real(kind=dprec)    :: fb,fc,fcc,cc
    real(kind=dprec)    :: eps,fvalue

    character(len=*) :: output_file

    !set tolerance
    eps = eps_lstqstbracket

    !initial checks
    if(flmax.le.dmin1(step0,tol))then
       write(ts_stdout,*) 'Error: Bracket Too Small ','No Room For Movement'
       call tssearch_utils_io_abort('Error detected in tssearch_algor_lstqstbracket')
    endif

    !initialize
    step_sk = 0.0_dprec
 
    !are you calling this routine for the first time or not ?
    if(ifirst.eq.0) then

         a=0.0_dprec
         b=step0
         fa=f0

         !interpolate and optimize coordinates
         call tssearch_algor_optipnt(nat,ndim,period,b,xx,xmin1,xmin2,pqst,xqst,&
                                              iopt,idfix)

         !calculate the path-coordinate
         call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xx,fvalue)

         !calculate the energy and gradient
         call tssearch_utils_energy_gradient(ndim,xx,gdum,fb)

         !write structure to a cumulative file
         call tssearch_utils_write(xx,gdum,fb,neg,output_file,&
                       tsstat,maxloop,fvalue)

         !update energy-gradient counter
         neg = neg + 1

         !write path-coordinate and energy of structure to output file
         if (root_write) then
           call tssearch_utils_write_summary(8,ts_stdout,fb,0.0_dprec,0.0_dprec,fvalue)
         end if

         !check for fixed atome information
         do ii = 1,ndim
           gwrite(ii) = gdum(ii)     !preserve a record of the complete gradient list
           !nxg: 11 May 2002: comment out
           !if (idfix(ii).eq.1) then
           !  gdum(ii) = 0.0_dprec
           !end if
         end do

    else

         a=aa
         b=bb
         fa= faa
         fb= fbb

    endif

    !update step
    step=step0*stepinc
    ifailed=0

 1  continue

    c=b+step + step_sk
    if(c.ge.flmax) then
        if(ifailed.eq.0)then
          c=flmax
          ifailed=1
        else
          numextr=0
          return
        endif
    endif

    !interpolate and optimize coordinates
    call tssearch_algor_optipnt(nat,ndim,period,c,xx,xmin1,xmin2,pqst,xqst,&
                             iopt,idfix)

    !calculate path-coordinate
    call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xx,fvalue)

    !calculate the energy and gradient
    call tssearch_utils_energy_gradient(ndim,xx,gdum,fc)

    !write structure to a cumulative file
    call tssearch_utils_write(xx,gdum,fc,neg,output_file,&
                       tsstat,maxloop,fvalue)

    !update energy-gradient counter
    neg = neg + 1

    !write path-coordinate and energy of structure to output file
    if (root_write) then
       call tssearch_utils_write_summary(8,ts_stdout,fc,0.0_dprec,0.0_dprec,fvalue)
    end if

    !check for fixed atome information
    do ii = 1,ndim
         gwrite(ii) = gdum(ii)  !preserve a record of the complete gradient list
         if (idfix(ii).eq.1) then
           gdum(ii) = 0.0_dprec
         end if
    end do

    !check if energy calculation failed
    if(scf_fail) then
          step_sk = step_sk + step
          go to 111
    else
          step_sk = 0.0_dprec
    endif

    !check energy differences for triplet
    if(((fc-fb).gt.eps .and. (fb-fa).lt.-eps .and. numextr.eq.1).or. &
            ((fc-fb).lt.-eps.and.(fb-fa).gt.eps.and.numextr.eq.-1))then

    aa=a
    bb=b
    cc=c
    faa=fa
    fbb=fb
    fcc=fc
    return

    else 

    a=b
    b=c
    fa=fb
    fb=fc

    endif

111 continue
    step=step*stepinc
    goto 1

    end subroutine tssearch_algor_lstqstbracket

    subroutine tssearch_algor_lstqstbrent(nat,ndim,idfix,xmin1,xmin2,iopt, &
               pqst,xqst,fl0,a,b,f0,fa,fb,numextr,tol,flopt,fopt, &
               period,xdum,gdum,neg,xclean1,xclean2,emax_old,output_file)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! this is a modification of function brent from 'Numerical Recipes':
    !
    ! INPUT:
    !  ndim - dimension of vectors xvec,pvec
    !  xvec(ndim) - current point
    !  pvec(ndim) - current direction of line search (line: xvec+fl*pvec)
    !  a < fl0  < b - bracketing triplet
    !  f0 - function value at fl0 
    !  fa,fb - function values at a and b
    !  numextr - (-1: bracketing of max),(1:bracketing of min)
    !  tol - accuracy for fl (NOTE difference with Num. Rec. where tol is for
    !                         function, not argument, as here)
    ! OUTPUT:
    !  flopt - optimal fl
    !  fopt - corresponding function value
    !  a<flopt<b - new bracketing triplet
    !  fa,fb - updated
    !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Sept 20 2000                                   !
    !=========================================================================!
    use tssearch_utils, only: root_write, ts_stdout, tssearch_utils_write, &
         tssearch_utils_write_summary, tssearch_utils_energy_gradient
 
    implicit none

    integer        :: ndim,iter,itmax,numextr

    real(kind=dprec)  :: xm,x,r,q,p,etemp,d,u,fl0,fx,f0,cgold,emax_old
    real(kind=dprec)  :: fa,fb,w,a,fw,v,b,fv,e
    real(kind=dprec)  :: xmin1(ndim),xmin2(ndim)
    real(kind=dprec)  :: xclean1(ndim),xclean2(ndim)
    real(kind=dprec)  :: flopt,fopt,fu,eps

    parameter(itmax=50,cgold=0.381966_dprec)

    logical          :: period

    integer          :: iopt,ii
    integer          :: neg
    integer          :: nat

    real(kind=dprec)    :: tol
    real(kind=dprec)    :: gdum(ndim)
    real(kind=dprec)    :: xqst(ndim),pqst
    real(kind=dprec)    :: xdum(ndim)

    integer          :: idfix(ndim)

    real(kind=dprec)    :: tol1, tol2
    real(kind=dprec)    :: fvalue
    character(len=*) :: output_file

    !set tolerance
    eps = eps_lstqstbrent
    d = 0

    !setup for maximization
    x=fl0
    fx=f0*numextr
    if((fa*numextr - fb*numextr) .lt. -eps) then
       w=a
       fw=fa*numextr
       v=b
       fv=fb*numextr
    else
       w=b
       fw=fb*numextr
       v=a
       fv=fa*numextr
    endif

    tol1=tol
    tol2=2.*tol1
    e=dmax1(tol1*5.0_dprec,0.25_dprec*(b-a))

    do iter=1,itmax
       xm=0.5_dprec*(a+b)

       if(dabs(x-xm).le.(tol2-0.5_dprec*(b-a))) goto 3
       if(dabs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-W)*r
          q=2.0*(q-r)
          if(q.gt.0.0_dprec) p = -p
          q=dabs(q)
          etemp=e
          e=d
          if(dabs(p).ge.dabs(0.5_dprec*q*etemp).or.p.le.q*(a-x).or. p.ge.q*(b-x)) goto 1
          d=p/q
          u=x+d
          
          if(u-a.lt.tol2.or.b-u.lt.tol2) d= dsign(tol1,xm-x)
          if(u-a.lt.tol2) then
             d=a-x+tol2
          elseif (b-u.lt.tol2) then
             d=b-x-tol2
          endif
          goto 2
       endif
1      if(x.ge.xm) then
          e=a-x
       else
          e=b-x
       endif
       d=cgold*e
       
2      if(dabs(d).ge.tol1) then
          u=x+d
       else
          u=x+dsign(tol1,d)
       endif
       
       
       !interpolate and optimize coordinates
       call tssearch_algor_optipnt(nat,ndim,period,u,xdum,xmin1,xmin2,pqst,xqst,iopt,idfix)
       
       !calculate path-coordinate
       call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xdum,fvalue)
       
       !calculate the energy and gradient
       call tssearch_utils_energy_gradient(ndim,xdum,gdum,fu) 
       
        !write structure to record
       call tssearch_utils_write(xdum,gdum,fu,neg,output_file,&
            tsstat,maxloop,fvalue)
       
       !update energy counter
       neg = neg + 1
       
       !write path-coordinate and energy of structure to output file
       if (root_write) then
          call tssearch_utils_write_summary(8,ts_stdout,fu,0.0_dprec,0.0_dprec,fvalue)
       end if
       
       !check fixed atom information
       do ii = 1,ndim
          gwrite(ii) = gdum(ii) !preserve a record of the complete gradient list
          !nxg: 11 May 2002: comment out
          !if (idfix(ii).eq.1) then
          !  gdum(ii) = 0.0_dprec
          !end if
       end do
       
       !this if block is needed only in transition state searches
       !it interrupts the search for the maximum if the maximum on 
       !the previous iteration is already exceeded - the step will 
       !be discarded anyway
       
       if((fu-emax_old).ge.eps)then
          flopt=u
          fopt=fu
          return
       endif
       
       fu=fu*numextr
       
       if((fu - fx).le.-eps) then
          if(u.ge.x) then
             a=x
             fa=fx*numextr
          else
             b=x
             fb=fx*numextr
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
          
       else
          
          if(u.lt.x) then
             a=u
             fa=fu*numextr
          else
             b=u
             fb=fu*numextr
          endif
          
          if((fu - fw).le.-eps .or. w.eq.x) then
             v=w
             fv=fw
             w=u
             fw=fu
          else if((fu - fv).le.-eps .or. v.eq.x.or.v.eq.w) then
             v=u
             fv=fu
          endif
          
       endif
       
    enddo

       if (iter.ge.itmax) then
          if (root_write) then
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write (ts_stdout,*)    ' Exceeded maximum iterations in tssearch_algor_lstqstbrent'
          write (ts_stdout,*)    ' Returning best estimate'
          write(ts_stdout,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          end if
          flopt=u
          fopt=fu
          return
       end if

   3   flopt=x

       !interpolate and optimize coordinates
       call tssearch_algor_optipnt(nat,ndim,period,flopt,xdum,xmin1,xmin2,pqst,xqst,iopt,idfix)

       !calculate path-coordinate
       call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xdum,fvalue)
  
       !calculate the energy and gradient
       call tssearch_utils_energy_gradient(ndim,xdum,gdum,fx)

       !write record
       call tssearch_utils_write(xdum,gdum,fx,neg,output_file,&
                       tsstat,maxloop,fvalue)

       fopt=fx

       !check fixed atom information
       do ii = 1,ndim
           gwrite(ii) = gdum(ii) !preserve a record of the complete gradient list
           !nxg: 11 May 2002: comment out
           !if (idfix(ii).eq.1) then
           !  gdum(ii) = 0.0_dprec
           !end if
       end do

       !update energy-gradient counter
       neg = neg + 1

       !write path-coordinate and energy of structure to output file
       if (root_write) then
          call tssearch_utils_write_summary(8,ts_stdout,fx,0.0_dprec,0.0_dprec,fvalue)
       end if

       return
    end subroutine tssearch_algor_lstqstbrent


    subroutine tssearch_algor_linemin(ndim,idfix,f,g,x,p,f_move, &
                     angmax,angle,istatus,neg) 
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Sept 25 2000                                   !
    !=========================================================================!
    use tssearch_utils, only: tssearch_utils_energy_gradient, &
         tssearch_utils_io_abort

    implicit none
  
    integer          :: i,intmax,istatus,intpln,ndim,idfix(ndim)
    integer          :: neg
    integer          :: niter,nitmax

    real(kind=dprec)    :: epsilon,stpmax,stpmin,gamma_ratio
    real(kind=dprec)    :: angmax,angle,x(ndim),g(ndim),p(ndim),s(ndim)
    real(kind=dprec)    :: f,f_move,s_norm,g_norm,cosang
    real(kind=dprec)    :: step,parabola,cube,cubstp,sss,ttt,step_new
    real(kind=dprec)    :: f_0,f_1,f_a,f_b,f_c,sg_0,sg_1,sg_a,sg_b,sg_c
    real(kind=dprec)    :: eps

    !initialize necessary parameters
    istatus = 0
    epsilon = 0.1_dprec
    eps =  eps_linemin
    stpmax = stpmax_linemin
    stpmin = stpmin_linemin

    if (angmax .eq. 0.0_dprec)  angmax = 180.0_dprec
    intmax = 5         !maximum number of interpolation steps allowed
    nitmax = 20        !maximum number of iterations allowed
    niter = 0

    !write some information to the output file
    !write(ts_stdout,*)
    !write(ts_stdout,'(a,f18.4)')   ' Line Search Maximum Step                      : ',stpmax
    !write(ts_stdout,*) 'rescale :',rescale
    !write(ts_stdout,*)

    !setup search vector
    do i = 1,ndim
      s(i) = p(i)
    end do

    !calculate the length of the gradient and the search vector
    g_norm = 0.0_dprec
    s_norm = 0.0_dprec
    do i = 1,ndim
        gwrite(i) = g(i) !preserve a record of the complete gradient list
        if (idfix(i).eq.1) then
              g(i) = 0.0_dprec
        end if
        g_norm = g_norm + g(i)*g(i)
        s_norm = s_norm + s(i)*s(i)
    end do
    g_norm = dsqrt(g_norm)
    s_norm = dsqrt(s_norm)

    !normalize the search vector and calculate projected gradient
    f_0 = f
    sg_0 = 0.0_dprec
    do i = 1,ndim
        s(i) = s(i) / s_norm
        sg_0 = sg_0 + s(i)*g(i)
    end do

    !check the angle between the search direction and negative gradient vectors
    cosang = -sg_0 / g_norm
    cosang = min(1.0_dprec,max(-1.0_dprec,cosang))
    angle = acos(cosang) * 180.0_dprec/pi
   
    !check if the angle is greater that 180 degs
    if (angle .gt. angmax) then
         istatus = 3
         return
    end if

    !set the initial stepsize
    step = 2.0_dprec * abs(f_move/sg_0)
    step = min(step,s_norm)
    if (step .gt. stpmax)  step = stpmax
    if (step .lt. stpmin)  step = stpmin

    !parabolic extrapolation
 10 continue         !outer iteration loop

    niter = niter + 1

    intpln = 0
    f_b = f_0
    sg_b = sg_0

    !update last point with latest and take another step
 20 continue

    !check number of iterations
    !if (niter.gt.nitmax) then
    !     istatus = 1          ! return with the best available point
    !     return
    !end if

    !update iteration counter
    !niter = niter + 1

    f_a = f_b
    sg_a = sg_b
    step_new = step

!   call tssearch_utils_safe(ndim,nat,period,xclean1,xclean2,x,s,idfix,step_new,1,1)

    !check fixed atom information 
    do i = 1,ndim
       if (idfix(i).eq.1) then
             x(i) = x(i)
       else
             x(i) = x(i) + step_new*s(i)
       end if
    end do

    !calculate function and gradient and test condition
    call tssearch_utils_energy_gradient(ndim,x,g,f_b)
    !f = f_b

    !update energy-gradient counter
    neg = neg + 1

    sg_b = 0.0_dprec
    do i = 1,ndim
       gwrite(i) = g(i)          !preserve a record of the complete gradient list
       if (idfix(i).eq.1) then
             g(i) = 0.0_dprec
       end if
       sg_b = sg_b + s(i)*g(i)
    end do
    gamma_ratio = dabs(sg_b/sg_0)

    !test the convergence condition
    if (gamma_ratio.le.epsilon .and. (f_b - f_a).lt.-eps) then
             f = f_b
             if (istatus .eq. 0)  istatus = 1            ! minimization successful
             return
    end if

    !start interpolation if slope changes or function increases
    if (sg_a*sg_b.lt.0.0_dprec .or. (f_a - f_b).lt.-eps)  goto 30
          step = 2.0_dprec * step              
          if (sg_b .gt. sg_a) then
               parabola = (f_a-f_b) / (sg_b-sg_a)
               if (parabola .gt. 2.0_dprec*step)  parabola = 2.0_dprec * step
               if (parabola .lt. 0.5_dprec*step)  parabola = 0.5_dprec * step 
               step = parabola
          end if
          if (step .gt. stpmax)  step = stpmax
          goto 20

    !cubic interpolation procedure
 30 continue
         intpln = intpln + 1

         !check number of iterations
         !if (niter.gt.nitmax .or. intpln.gt.intmax) then
         !  istatus = 1          ! return with the best available point
         !  return
         !end if

         sss = 3.0_dprec*(f_b-f_a)/step - sg_a - sg_b
         ttt = sss*sss - sg_a*sg_b
         if (ttt .lt. 0.0_dprec) then
               f = f_b
               istatus = 5          ! interpolation problem: safeguard against negative square root
               return
         end if
         ttt = sqrt(ttt)

         cube = step * (sg_b+ttt+sss)/(sg_b-sg_a+2.0_dprec*ttt)
         if (cube.lt.0.0_dprec .or. cube.gt.step) then
               f = f_b
               istatus = 5          ! interpolation problem
               return
         end if

         !update step
         step_new = cube

!        call tssearch_utils_safe(ndim,nat,period,xclean1,xclean2,x,s,idfix,step_new,-1,1)

         !check fixed atom information and update coordinate
         do i = 1,ndim
           if (idfix(i).eq.1) then
             x(i) = x(i)
           else
             x(i) = x(i) - step_new*s(i)
           end if
         end do

         !get new function value and gradient and test condition ---
         call tssearch_utils_energy_gradient(ndim,x,g,f_c)

         !update energy-gradient counter
         neg = neg + 1

         !check convergence condition again
         sg_c = 0.0_dprec
         do i = 1,ndim
              gwrite(i) = g(i) !preserve a record of the complete gradient list
              if (idfix(i).eq.1) then
                g(i) = 0.0_dprec
              end if
              sg_c = sg_c + s(i)*g(i)
         end do
         gamma_ratio = dabs(sg_c/sg_0)
         if (gamma_ratio .le. epsilon) then
              f = f_c
              if (istatus .eq. 0)  istatus = 1         ! successfully minimized
              return
         end if

         !next pair of bracketing points
         if ((f_c - f_a) .le. -eps .or. (f_c - f_b) .le. -eps) then
               cubstp = min(abs(cube),abs(step-cube))
               if (cubstp.ge.stpmin .and. intpln.lt.intmax) then

                 if (sg_a*sg_b .lt. 0.0_dprec) then

                    if (sg_a*sg_c .lt. 0.0_dprec) then
                        f_b = f_c
                        sg_b = sg_c
                        step = step - cube
                    else
                        f_a = f_c
                        sg_a = sg_c
                        step = cube

                        step_new = cube
!                       call tssearch_utils_safe(ndim,nat,period,xclean1,xclean2,x,s,idfix,step_new,1,1)

                        do i = 1,ndim
                         if (idfix(i).eq.1) then
                           x(i) = x(i)
                         else
                           x(i) = x(i) + step_new*s(i)
                         end if
                        end do

                    end if
            
                 else
           
                    if (sg_a*sg_c.lt.0.0_dprec .or. (f_a - f_c) .le. -eps) then
                       f_b = f_c
                       sg_b = sg_c
                       step = step - cube
                    else
                       f_a = f_c
                       sg_a = sg_c
                       step = cube

                       step_new = cube
!                      call tssearch_utils_safe(ndim,nat,period,xclean1,xclean2,x,s,idfix,step_new,1,1)

                       do i = 1,ndim
                         if (idfix(i).eq.1) then
                           x(i) = x(i)
                         else
                           x(i) = x(i) + step_new*s(i)
                         end if
                       end do

                    end if

                 end if

                 goto 30
               end if 
         end if  

         !interpolation failed: reset to best current point
            
         f_1 = min(f_a,f_b,f_c)

!        if (f_1 .eq. f_a) then
         if (abs(f_1 - f_a) .le. eps) then

              sg_1 = sg_a

              step_new = cube-step
!             call tssearch_utils_safe(ndim,nat,period,xclean1,xclean2,x,s,idfix,step_new,1,1)

              do i = 1,ndim
                if (idfix(i).eq.1) then
                  x(i) = x(i)
                else
                  x(i) = x(i) + step_new*s(i)
                end if
              end do

!        else if (f_1 .eq. f_b) then
         else if (abs(f_1 - f_b) .le. eps) then

              sg_1 = sg_b

              step_new = cube
!             call tssearch_utils_safe(ndim,nat,period,xclean1,xclean2,x,s,idfix,step_new,1,1)

              do i = 1,ndim
                if (idfix(i).eq.1) then
                  x(i) = x(i)
                else
                  x(i) = x(i) + step_new*s(i)
                end if
              end do

!        else if (f_1 .eq. f_c) then
         else if (abs(f_1 - f_c) .le. eps) then

              sg_1 = sg_c

         else
             sg_1=0 ! Quit
             call tssearch_utils_io_abort('Failure in TSSEARCH_ALGOR_LINEMIN')
         end if 

         !restart from best point with smaller stepsize

         if ((f_1 - f_0) .gt. eps) then

              call tssearch_utils_energy_gradient(ndim,x,g,f)
              neg = neg + 1

              do i=1,ndim
               gwrite(i) = g(i) !preserve a record of the complete gradient list
               if (idfix(i).eq.1) then
                 g(i) = 0.0_dprec
               end if
              end do

              istatus = 5          ! interpolation error
              return
         end if

         f_0 = f_1
         sg_0 = sg_1
         if (sg_1 .gt. 0.0_dprec) then
            do i = 1,ndim
              s(i) = -s(i)
            end do
         end if
         step = max(cube,step-cube) / 10.0_dprec

         if (step .lt. stpmin)  step = stpmin

         !if restarted once, return with best point
         if (istatus .eq. 2) then
             call tssearch_utils_energy_gradient(ndim,x,g,f)
             neg = neg + 1
             do i=1,ndim
               gwrite(i) = g(i) !preserve a record of the complete gradient list
               if (idfix(i).eq.1) then
                 g(i) = 0.0_dprec
               end if
             end do
             istatus = 4           ! bad interpolation
             return
         else
             istatus = 2           ! perform the search again
             if (niter.lt.nitmax) then    ! to make sure we don't iterate for ever
               goto 10
             else 
              !istatus = 1          ! return with best available point
              !call tssearch_utils_energy_gradient(ndim,x,g,f) !get its energy
              !neg = neg + 1
              istatus = 2
              return
             end if
         end if

    end subroutine tssearch_algor_linemin

    subroutine tssearch_algor_dflmin(ndim,p,tol_in,tol_out)
    !=========================================================================!
    !calculate accuracy up to which flopt should be calculated for given p    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Dec 4 2001                                     !
    !=========================================================================!
    use tssearch_utils, only : ts_stdout, tssearch_utils_io_abort

    implicit none

    integer,intent(in) :: ndim
    real(kind=dprec),intent(in) :: p(ndim),tol_in
    real(kind=dprec),intent(out) :: tol_out
 
    real(kind=dprec) :: pmax
    integer :: i

    !find maximum value in array
!   pmax=VecMax(ndim,p,idum)
    pmax = 0.0_dprec
    do i = 1 , ndim
      pmax = max(pmax,abs(p(i)))
    end do

    if(pmax.gt.1e-9_dprec) then
        tol_out=tol_in/pmax
    else
        write(ts_stdout,'(a)') 'Error: Geometries are too close'
        !write(ts_stdout,'(a,f11.6)') ' Maximum component of p is ',pmax
        write(ts_stdout,'(a)') 'Please check the structures'
        call tssearch_utils_io_abort('Error detected in tssearch_algor_dflmin')
    endif

    return
  end subroutine tssearch_algor_dflmin

  subroutine tssearch_algor_scan_test(nat,ndim,period,ngeom,x0,x1,xm,pm, &
                                iopt,idfix,xclean1,xclean2,neg,output_file) 
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Dec 10 2001                                    !
    !=========================================================================!
    use tssearch_utils, only: ts_stdout, tssearch_utils_write, &
         tssearch_utils_energy_gradient

    implicit none

    ! Arguments
    integer :: nat,ndim,ngeom,idfix(:),iopt,neg
    logical :: period
    real(kind=dprec) ::x0(:),x1(:),xm(:),pm,xclean1(:),xclean2(:)
    character(len=*),  intent(in)    :: output_file

    ! Local variables
    integer :: ii
    real(kind=dprec) :: fscan,xx(ndim),gg(ndim),ff,fvalue
 
    ! Scan test
    if (iopt == 1 ) then
     write(ts_stdout,*)
     write(ts_stdout,*) ' Start a LST scan'
     write(ts_stdout,*)
     write(ts_stdout,*) 'tsstat ',tsstat
     write(ts_stdout,*) 'maxloop ',maxloop
     write(ts_stdout,*) 'output_file ',output_file
    else
     write(ts_stdout,*)
     write(ts_stdout,*) ' Start a QST scan'
     write(ts_stdout,*)
     write(ts_stdout,*) 'tsstat ',tsstat
     write(ts_stdout,*) 'maxloop ',maxloop
     write(ts_stdout,*) 'output_file ',output_file
    end if

    do ii = 0, ngeom
       fscan = 1.0_dprec/ngeom * ii * 1.0_dprec
       call tssearch_algor_optipnt(nat,ndim,period,fscan,xx,x0,x1,pm,xm,iopt,idfix)
       call tssearch_algor_pathcoord(nat,ndim,period,xclean1,xclean2,xx,fvalue)
       call tssearch_utils_energy_gradient(ndim,xx,gg,ff)
       call tssearch_utils_write(xx,gg,ff,neg,output_file,&
                       tsstat,maxloop,fvalue)
    end do

  end subroutine tssearch_algor_scan_test

  subroutine tssearch_algor_periodic_test
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind,  Dec 05 2001                                    !
    !=========================================================================!
    use tssearch_utils, only: ts_stdout, gamma_cos, gamma_sin, alength, &
         blength, clength, beta_term, beta_sin, beta_cos, alpha_cos, gamma_term

    implicit none

    ! Local variables
    real(kind=dprec)  :: xbox,ybox,zbox,xbox2,ybox2,zbox2

    !setup some local variables
    xbox  = alength
    ybox  = blength
    zbox  = clength
    xbox2 = 0.5_dprec*xbox
    ybox2 = 0.5_dprec*ybox
    zbox2 = 0.5_dprec*zbox

    write(ts_stdout,*) 'xxxxxxxxxxxxxxxxxxxxx'
    write(ts_stdout,*) alpha_cos
    write(ts_stdout,*) beta_sin
    write(ts_stdout,*) beta_cos
    write(ts_stdout,*) gamma_sin
    write(ts_stdout,*) gamma_cos
    write(ts_stdout,*) beta_term
    write(ts_stdout,*) gamma_term
    write(ts_stdout,*) xbox
    write(ts_stdout,*) ybox
    write(ts_stdout,*) zbox
    write(ts_stdout,*) xbox2
    write(ts_stdout,*) ybox2
    write(ts_stdout,*) zbox2
    write(ts_stdout,*) 'xxxxxxxxxxxxxxxxxxxxx'

    return
    end subroutine tssearch_algor_periodic_test
    
    subroutine tssearch_algor_message_test
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, Nov 27 2001                               !
    !=========================================================================!
    use tssearch_utils, only : ts_stdout

    implicit none

    write(ts_stdout,*) 'Have a good day!'

    return
    end subroutine tssearch_algor_message_test

end module tssearch_algor
