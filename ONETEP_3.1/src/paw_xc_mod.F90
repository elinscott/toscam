! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!       Projector-Augmented Wave Exchange-Correlation module     !
!                                                                !
! This module
!----------------------------------------------------------------!
! This module was written by Nicholas Hine in May-June 2010.     !
!================================================================!

module paw_xc

  use constants, only: DP, stdout

  implicit none

  private

  public :: paw_xc_pot_rad_LM
  public :: paw_xc_exc_dc

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_xc_pot_rad_LM(vxc,exc,den,r,rab,rmax,npts,npts_max,ns, &
       nLM_max,lmax,den_work,pot_work,inter)

    !=================================================================!
    ! This subroutine calculates the sphere exchange-correlation      !
    ! potential v_xc for a radial charge distribution given as a sum  !
    ! over its moments LM.                                            !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 02/06/10.                           !
    !=================================================================!

    use constants, only: SQRT_PI
    use gaunt_coeff, only: realgaunt
    use services, only: services_radial_integral_rmax
    use xc, only: xc_radial

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    integer, intent(in) :: ns
    integer, intent(in) :: nLM_max
    integer, intent(in) :: npts_max
    integer, intent(in) :: lmax
    real(kind=DP), intent(out) :: vxc(npts_max,ns,nLM_max)
    real(kind=DP), intent(inout) :: exc
    real(kind=DP), intent(in) :: den(npts_max,ns,nLM_max)
    real(kind=DP), intent(in) :: r(npts_max)
    real(kind=DP), intent(in) :: rab(npts_max)
    real(kind=DP), intent(in) :: rmax
    real(kind=DP), intent(out) :: den_work(npts_max,ns,ns,4)
    real(kind=DP), intent(out) :: pot_work(npts_max,ns,ns,6)
    real(kind=DP), intent(out) :: inter(npts_max)

    ! Local Variables
    integer :: i,is1,is2
    integer :: lup,mup,iLM
    real(kind=DP), parameter :: sqrt_4pi = 2.0_DP*SQRT_PI
    real(kind=DP), parameter :: inv_sqrt_4pi = 1.0_DP/sqrt_4pi
    real(kind=DP), parameter :: delta = 1.0e-4_DP
    real(kind=DP), parameter :: told = 1.0e-14_DP
    logical :: full_expansion

    ! Initialisation
    vxc(:,:,:) = 0.0_DP
    den_work(:,:,:,:) = 0.0_DP
    pot_work(:,:,:,:) = 0.0_DP
    inter(:) = 0.0_DP
    full_expansion = .false.

    ! TODO
    ! pot_work(:,s1, 1,1) = v^s1_xc(n_00)
    ! pot_work(:,s1,s2,2) = (d v^s1_xc(n)/d n^s2)(n_s)
    ! pot_work(:,s1,s2,3) = (d^2 v^s1_xc / d n^s2^2)(n_s)
    ! pot_work(:,:,:,4) = exc(n_s)
    ! den_work(:,s1, 1,1) = \sum_L>0 (n^s1_LM)(n^s2_LM)
    ! den_work(:,s1,s2,2) = \sum_L'>0 \sum_L''>0 n^s1_L'M' n^s2_L''M'' G^LM_L'M'L''M''
    ! den_work(:,s1,s2,3) = \sum_L>0 \sum_L'>0 \sum_L''>0 n^s1_LM n^s2_L'M' n^s3_L''M'' G^LM_L'M'L''M''

    ! =============================================!
    ! Find v_xc(n_s(r)*(1+delta))                  !
    ! =============================================!

    ! Place (1+delta)*n^s1_s(r) in den_work(:,s1,1,1)
    do is1=1,ns
       ! Copy n_s(r) to den_work
       do is2=1,ns
          den_work(1:npts,is2,1,1) = den(1:npts,is2,1)*inv_sqrt_4pi
       end do
       ! Scale relevant spin-channel
       den_work(1:npts,is1,1,1) = den_work(1:npts,is1,1,1)*(1.0_DP + delta)
       ! Find v_xc(n_s(r)*(1+delta)), put into pot_work(:,:,:,5)
       call xc_radial(npts,npts_max,ns,den_work(1:npts_max,1:ns,1,1), &
            den_work(1:npts_max,1:ns,1,2),pot_work(1:npts_max,1:ns,is1,5))
    end do

    ! =============================================!
    ! Find v_xc(n_s(r)*(1-delta))                  !
    ! =============================================!

    ! Place (1-delta)*n^s1_s(r) in den_work(:,s1,1,1)
    do is1=1,ns
       ! Copy n_s(r) to den_work
       do is2=1,ns
          den_work(1:npts,is2,1,1) = den(1:npts,is2,1)*inv_sqrt_4pi
       end do
       ! Scale relevant spin-channel
       den_work(1:npts,is1,1,1) = den_work(1:npts,is1,1,1)*(1.0_DP - delta)
       ! Find v_xc(n_s(r)*(1-delta)), put into pot_work(:,:,:,6)
       call xc_radial(npts,npts_max,ns,den_work(1:npts_max,1:ns,1,1), &
            den_work(1:npts_max,1:ns,1,2),pot_work(1:npts_max,1:ns,is1,6))
    end do

    ! =============================================!
    ! Find v_xc(n_s(r)) and e_xc(n_s(r)            !
    ! =============================================!

    ! Place n_s(r) in den_work(:,:,:,1)
    do is1=1,ns
       den_work(1:npts,is1,1,1) = den(1:npts,is1,1)*inv_sqrt_4pi
    end do

    ! Find v_xc(n_s(r)) and e_xc(n_s(r)), put into pot_work(:,:,:,1 and 4)
    call xc_radial(npts,npts_max,ns,den_work(1:npts_max,1:ns,1,1), &
         den_work(1:npts_max,1:ns,1,2),pot_work(1:npts_max,1:ns,1,1), &
         exc_rad=pot_work(1:npts_max,1,1,4))

    ! Put 1/(2*(\sum_s1 n^s1_s(r))*delta) into den_work(:,:,:,2)
    do i=1,npts
       if (sum(den_work(i,1:ns,1,1)) > told) then
          den_work(i,1,1,2) = 1.0_DP / (2.0_DP*sum(den_work(i,1:ns,1,1))*delta)
       else
          den_work(i,1,1,2) = 0.0_DP
       end if
    end do

    ! Turn energy e_xc[n_s(r)] back into into energy density eps_xc[n_s(r)]
    pot_work(1:npts,1,1,4) = pot_work(1:npts,1,1,4) * &
         den_work(1:npts,1,1,2)*delta*2.0_DP

    ! =============================================!
    ! Calculate d v_xc / dn and d^2 v_xc / dn^2    !
    ! =============================================!

    ! Put 1/(2*n^s1_s(r)*delta) into den_work(:,s1,1,2)
    do is1=1,ns
       do i=1,npts
          if (den_work(i,is1,1,1) > told) then
             den_work(i,is1,1,2) = 1.0_DP / (2.0_DP*den_work(i,is1,1,1)*delta)
          else
             den_work(i,is1,1,2) = 0.0_DP
          end if
       end do
    end do

    do is1=1,ns
       do is2=1,ns
          ! Find (dv^s1_xc / dn^s2)_{n=n_s}
          pot_work(1:npts,is1,is2,2) = &
               (pot_work(1:npts,is1,is2,5) - pot_work(1:npts,is1,is2,6)) * &
               den_work(1:npts,is2,1,2)

          ! Find (d2v^s1_xc / dn^s2^2)_{n=n_s}
          pot_work(1:npts,is1,is2,3) = &
               (pot_work(1:npts,is1,is2,5) + pot_work(1:npts,is1,is2,6) &
               - 2.0_DP*pot_work(1:npts,is1,1,1))*den_work(1:npts,is2,1,2)**2*4.0_DP

       end do
    end do

    ! =============================================!
    ! Calculate other density sums in workspace    !
    ! =============================================!

    ! Calculate den_work(:,s1,s2,1) = \sum_L>0 (n^s1_LM)(n^s2_LM)
    call internal_sum_nLM2

#if 0
    write(76,*)
    do is1=1,ns
       do is2=1,ns
           write(76,'(a)') &
                '#                 r,                den,&
                &                vxc,              epsxc,&
                &           dv_xc/dn,        d2v_xc/dn^2,&
                &           sum_nLM2,    1/(2*den*delta),&
                &      v_xc(1+delta),       v_xc(1-delta)'
          write(76,*) '# is1,is2: ',is1,is2
          do i=1,npts
             write(76,'(2f20.10,f20.15,2f20.10,e20.12,2f20.10,2f20.15)') &
                   r(i),den(i,is1,1), &
                   pot_work(i,is1,1,1),pot_work(i,1,1,4), &
                   pot_work(i,is1,is2,2),pot_work(i,is1,is2,3), &
                   den_work(i,is1,1,1),den_work(i,is2,1,2), &
                   pot_work(i,is1,is2,5),pot_work(i,is1,is2,6)
          end do
          write(76,*)
       end do
    end do
#endif

    ! =============================================!
    ! Calculate XC potential for L=0 component     !
    ! =============================================!

    do is1=1,ns
       vxc(1:npts,is1,1) = sqrt_4pi*pot_work(1:npts,is1,1,1)
       do is2=1,ns
          vxc(1:npts,is1,1) = vxc(1:npts,is1,1) + &
               0.5_DP/sqrt_4pi * pot_work(1:npts,is1,is2,3) * &
               den_work(1:npts,is2,is2,1)
       end do
       if (ns==2) then
          is2 = 1 + modulo(is1,2)
          vxc(1:npts,is1,1) = vxc(1:npts,is1,1) + &
               2.0_DP * 0.5_DP/sqrt_4pi * pot_work(1:npts,is1,is2,3) * &
               den_work(1:npts,is2,is2,1)
       end if
    end do

    ! =============================================!
    ! Calculate XC potential for L>0 components    !
    ! =============================================!

    do is1=1,ns
       do lup=1,lmax
          do mup=-lup,lup
             iLM = lup*lup + lup + 1 + mup
             do is2=1,ns
                vxc(1:npts,is1,iLM) = vxc(1:npts,is1,iLM) + &
                     den(1:npts,is1,iLM)*pot_work(1:npts,is1,is2,2)
             end do
             if (full_expansion) then
                call internal_sum_nLMnLM_G
                do is2=1,ns
                   vxc(1:npts,is1,1) = vxc(1:npts,is1,1) + &
                        0.5_DP * pot_work(1:npts,is1,is2,3) * &
                        den_work(1:npts,is2,is2,2)
                end do
                if (ns==2) then
                   is2 = 1 + modulo(is1,2)
                   vxc(1:npts,is1,1) = vxc(1:npts,is1,1) + &
                        2.0_DP * 0.5_DP * pot_work(1:npts,is1,is2,3) * &
                        den_work(1:npts,is2,is2,2)
                end if
             end if
          end do
       end do
    end do

#if 0
    if (ns==2) then
       do lup=0,lmax
          do mup=-lup,lup
             iLM = lup*lup + lup + 1 + mup
             do i=1,npts
                write(77,'(i3,5f20.12)') iLM,r(i),den(i,1,iLM),den(i,2,iLM),vxc(i,1,iLM),vxc(i,2,iLM)
             end do
          end do
       end do
    else
       do lup=0,lmax
          do mup=-lup,lup
             do i=1,npts
                iLM = lup*lup + lup + 1 + mup
                write(77,'(i3,3f20.12)') iLM,r(i),den(i,1,iLM),vxc(i,1,iLM)
             end do
          end do
       end do
    end if
#endif

    ! =============================================!
    ! Calculate XC energy                          !
    ! =============================================!

    ! Calculate \int_0^rc dr n_00(r) eps_xc[n_00(r)](r) r^2
    if (ns==1) then
       den_work(1:npts,1,1,4) = sqrt_4pi*den(1:npts,1,1) * &
            pot_work(1:npts,1,1,4)*r(1:npts)**2
    else
       den_work(1:npts,1,1,4) = sqrt_4pi*(den(1:npts,1,1) + &
            den(1:npts,2,1)) * pot_work(1:npts,1,1,4)*r(1:npts)**2
    end if
    exc = services_radial_integral_rmax(npts,rab,r,rmax, &
          den_work(1:npts,1,1,4),inter(1:npts))

    ! Calculate den_work(:,3) = \sum_L>0 \sum_L'>0 \sum_L''>0
    !                           n^s1_LM n^s2_L'M' n^s2_L''M'' G^LM_L'M'L''M''
    if (full_expansion) call internal_sum_nLMnLMnLM_G

    ! Calculate \int_0^rc \sum_L>0 n_LM(r)*n_LM(r) dv_xc / dn r^2 dr
    !         + \int_0^rc (expr. above in dw(3)) * d^2 v_xc / dn^2 r^2 dr
    den_work(:,:,1,4) = 0.0_DP
    do is1=1,ns
       do is2=1,ns
          den_work(1:npts,is1,1,4) = den_work(1:npts,is1,1,4) + &
               0.5_DP*den_work(1:npts,is1,is2,1) * &
               pot_work(1:npts,is1,is2,2)*r(1:npts)**2
       end do
       if (full_expansion) then
          den_work(1:npts,is1,1,4) = den_work(1:npts,is1,1,4) + &
               (1.0_DP/6.0_DP) * den_work(1:npts,is1,1,3) * &
               pot_work(1:npts,is1,is2,2) * r(1:npts)**2
       end if
       exc = exc + services_radial_integral_rmax(npts,rab,r,rmax, &
            den_work(1:npts,is1,1,4),inter(1:npts))
    end do

contains

    ! Calculate den_work(:,s1,s2,1) = \sum_L>0 (n^s1_LM)(n^s2_LM)
    subroutine internal_sum_nLM2

      den_work(1:npts,:,:,1) = 0.0_DP
      do is1=1,ns
         do is2=1,ns
            do lup=1,lmax
               do mup=-lup,lup
                  iLM = lup*lup + lup + 1 + mup
                  den_work(1:npts,is1,is2,1) = den_work(1:npts,is1,is2,1) + &
                       den(1:npts,is1,iLM)*den(1:npts,is2,iLM)
               end do
            end do
         end do
      end do

    end subroutine internal_sum_nLM2

    ! Calculate den_work(:,s1,s2,2) = \sum_L'>0 \sum_L''>0 n^s1_L'M' n^s2_L''M'' G^LM_L'M'L''M''
    subroutine internal_sum_nLMnLM_G

      ! Local Variables
      integer :: lupp,mupp,iLMp
      integer :: luppp,muppp,iLMpp
      real(kind=DP) :: rglm

      den_work(1:npts,:,:,2) = 0.0_DP
      do is1=1,ns
         do is2=1,ns
            do lupp=1,lmax
               do mupp=-lupp,lupp
                  iLMp = lupp*lupp + lupp + 1 + mupp
                  do luppp=1,lmax
                     do muppp=-luppp,luppp
                        iLMpp = luppp*luppp + luppp + 1 + muppp
                        rglm = realgaunt(lup,mup,lupp,mupp,luppp,muppp)
                        !if (abs(rglm)>gaunt_tol) &
                             den_work(1:npts,is1,is2,2) = den(1:npts,is1,iLMp) * &
                             den(1:npts,is2,iLMpp) * rglm
                     end do
                  end do
               end do
            end do
         end do
      end do

    end subroutine internal_sum_nLMnLM_G

    ! Calculate den_work(:,3) = \sum_L>0 \sum_L'>0 \sum_L''>0 n^s1_LM n^s2_L'M' n^s2_L''M'' G^LM_L'M'L''M''
    subroutine internal_sum_nLMnLMnLM_G

       ! Local Variables
       integer :: lupp,mupp,iLMp
       integer :: luppp,muppp,iLMpp
       real(kind=DP) :: rglm

       den_work(1:npts,:,:,3) = 0.0_DP
       do is1=1,ns
          do is2=1,ns
             do lup=1,lmax
                do mup=-lup,lup
                   iLM = lup*lup + lup + 1 + mup
                   do lupp=1,lmax
                      do mupp=-lupp,lupp
                         iLMp = lupp*lupp + lupp + 1 + mupp
                         do luppp=1,lmax
                            do muppp=-luppp,luppp
                               iLMpp = luppp*luppp + luppp + 1 + muppp
                               rglm = realgaunt(lup,mup,lupp,mupp,luppp,muppp)
                               !if (abs(rglm)>gaunt_tol) &
                                    den_work(1:npts,is1,is2,3) = &
                                    den_work(1:npts,is1,is2,3) + &
                                    den(1:npts,is1,iLM)*den(1:npts,is2,iLMp) * &
                                    den(1:npts,is2,iLMpp) * rglm
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do

    end subroutine internal_sum_nLMnLMnLM_G

  end subroutine paw_xc_pot_rad_LM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_xc_exc_dc(exc_dc,vxc,den,r,rab,rmax,npts,npts_max, &
       ns,nLM,nLM_max,work,inter)

    !=================================================================!
    ! This subroutine calculates the double counting correction to    !
    ! Exchange correlation energy for a spherical charge distribution !
    ! and a spherical xc potential.                                   !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 15/06/10.                           !
    !=================================================================!

    use services, only: services_radial_integral_rmax

    implicit none

    ! Arguments
    integer, intent(in) :: nLM
    integer, intent(in) :: npts
    integer, intent(in) :: ns
    integer, intent(in) :: nLM_max
    integer, intent(in) :: npts_max
    real(kind=DP), intent(out) :: exc_dc
    real(kind=DP), intent(in) :: vxc(npts_max,ns,nLM_max)
    real(kind=DP), intent(in) :: den(npts_max,ns,nLM_max)
    real(kind=DP), intent(in) :: r(npts_max)
    real(kind=DP), intent(in) :: rab(npts_max)
    real(kind=DP), intent(in) :: rmax
    real(kind=DP), intent(out) :: work(npts_max)
    real(kind=DP), intent(out) :: inter(npts_max)

    ! Local Variables
    integer :: is
    integer :: iLM

    ! Calculate XC energy double-counting term for each spin
    exc_dc = 0.0_DP
    do is=1,ns
       do iLM=1,nLM
          work(1:npts) = &
               den(1:npts,is,iLM) * vxc(1:npts,is,iLM) * r(1:npts)**2
          exc_dc = exc_dc + services_radial_integral_rmax(npts,rab,r,rmax, &
               work(1:npts),inter(1:npts))
       end do
    end do

  end subroutine paw_xc_exc_dc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module paw_xc
