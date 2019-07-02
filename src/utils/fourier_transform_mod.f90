module fourier_transform_mod

   use genvar

   ! INTERFACE realft
   !  MODULE PROCEDURE realft_,realft__
   ! END INTERFACE

   ! INTERFACE fourier_transform_1
   !  MODULE PROCEDURE fourier_transform_1__,fourier_transform_1_
   ! END INTERFACE

! contains
!
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
!
!
!  !-----------------------------------------------------------------!
!  ! this program is the kernel for an optical propagation program.
!  ! it does a series of 2d ffts followed by a multiplication.
!  ! it is modeled after the AFWL program HELP or High Energy
!  ! Laser Propagation.
!  !
!  ! it does the 2d fft by first doing a collection of 1d ffts
!  ! then a transpose followed by a second collection of 1d ffts.
!  !
!  ! the sections of the program that are commented out represent
!  ! different ways of doing the operations.  the most interesting
!  ! addition is using the using the subroutine shuff to generate
!  ! a nonuniform ordering for accessing the array.
!  !
!  ! the routines fourier_transform_1 was taken from the book
!  ! "numerical recipes in fortran, 1st addition."
!  ! however, the authors of that book derived their routine
!  !-----------------------------------------------------------------!
!
!
!    subroutine two_d_fft
!    implicit none
!
!         integer,parameter:: rsize=1024
!         integer kkkk
!         integer i,j,k,ijk,isign
!         real(kind=DP),allocatable:: x(:)
!         integer, allocatable:: index(:)
!         complex(kind=DP), allocatable:: a(:,:)
!         complex(kind=DP) tmp
!         complex(kind=DP) factor
!         real(kind=DP) gen,fft1,fft2,trans,totf,fact
!         real(kind=DP) t0,t1,t2,t3,t4,t5
!         integer OMP_GET_MAX_THREADS
!
!         factor=rsize
!         factor=1.d0/(factor)
!         isign=1
!         gen=0
!         fft1=0
!         fft2=0
!         trans=0
!         totf=0
!         fact=0
!         allocate(a(rsize,rsize))
!         allocate(x(rsize),index(rsize))
!         call shuff(index,rsize,8)
!         t0=ccm_time ()
!         do j=1,rsize
!            call random_number(x)
!            do i=1,rsize
!                 a(i,j)=cmplx(x(i),0.d0)
!            enddo
!         enddo
!         write(*,'(("(",g20.10,",",g20.10,")"))')a(rsize/2+1,rsize-2)
!
!         do 10 ijk=1,20
!             t1=ccm_time ()
! !$OMP PARALLEL DO SCHEDULE (RUNTIME)
!             do i=1,rsize
!                  call fourier_transform_1((/(real(a(kkkk,i)),aimag(a(kkkk,i)),kkkk=1,size(a(:,1)))/),rsize,isign)
!             enddo
! !$OMP END PARALLEL DO
!             t2=ccm_time ()
! !$OMP PARALLEL DO SCHEDULE (RUNTIME) PRIVATE(i,j,k,tmp)
!             do k=1,rsize
!                 i=k
!                 do j=i,rsize
!                     tmp=a(i,j)
!                     a(i,j)=a(j,i)
!                     a(j,i)=tmp
!                 enddo
!             enddo
! !$OMP END PARALLEL DO
!             t3=ccm_time ()
! !$OMP PARALLEL DO SCHEDULE (RUNTIME)
!             do i=1,rsize
!                  call fourier_transform_1((/(real(a(kkkk,i)),aimag(a(kkkk,i)),kkkk=1,size(a(:,1)))/),rsize,isign)
!             enddo
! !$OMP END PARALLEL DO
!             t4=ccm_time ()
! !$OMP PARALLEL DO SCHEDULE (RUNTIME)
!             do j=1,rsize
!             do i=1,rsize
!                   a(i,j)=factor*a(i,j)
!                enddo
!             enddo
! !$OMP END PARALLEL DO
!             t5=ccm_time ()
!             gen=gen+t1-t0
!             fft1=fft1+t2-t1
!             fft2=fft2+t4-t3
!             trans=trans+t3-t2
!             totf=totf+t5-t1
!             fact=fact+t5-t4
!             isign=isign*(-1)
!             write(*,'(i3)',advance="no")ijk
!     10  continue
!         write(*,*)
!         write(*,'(("(",g20.10,",",g20.10,")"))')a(rsize/2+1,rsize-2)
!         write(*,*)"number of  transforms",ijk-1
!         write(*,'("      fft1 time= ",f9.4)')fft1
!         write(*,'(" transpose time= ",f9.4)')trans
!         write(*,'("      fft2 time= ",f9.4)')fft2
!         write(*,'("   scaling time= ",f9.4)')fact
!         write(*,'("    total time = ",f9.4)',advance="no")totf
!         write(*,'(" for matrix of rsize",i6)')rsize
!         write(*,*)"THREADS   = ",OMP_GET_MAX_THREADS()
!         stop
! end subroutine
!
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
!
!                !--------------------------!
!
!       subroutine fourier_transform_1__(data,nn,isign)
!       implicit none
!       integer    :: j,nn,isign,i
!       complex(kind=DP) :: data(:)
!       real(kind=DP)    :: data_(2*nn)
!        if(size(data)/=nn) stop 'error fourier transform sizes do not match'
!        do i=1,nn
!         data_(2*i-1)=real(data(i))
!         data_(2*i  )=aimag(data(i))
!        enddo
!        call fourier_transform_1_(data_,nn,isign)
!        do i=1,nn
!         data(i)=data_(2*i-1)+imi*data_(2*i)
!        enddo
!       end subroutine
!
!                !--------------------------!
!
!       subroutine fourier_transform_1_(data,nn,isign)
!       implicit none
!       integer            :: i,j,isign,nn,n,m,mmax,istep
!       real(kind=DP), parameter :: two_pi = 2.d0*pi
!       real(kind=DP)            :: wr,wi,wpr,wpi,wtemp,theta,tempr,tempi
!       real(kind=DP)            :: data(2*nn)
!       n=2*nn
!       j=1
!       do 11 i=1,n,2
!         if(j.gt.i)then
!           tempr=data(j)
!           tempi=data(j+1)
!           data(j)=data(i)
!           data(j+1)=data(i+1)
!           data(i)=tempr
!           data(i+1)=tempi
!         endif
!         m=n/2
! 1       if ((m.ge.2).and.(j.gt.m)) then
!           j=j-m
!           m=m/2
!         go to 1
!         endif
!         j=j+m
! 11    continue
!       mmax=2
! 2     if (n.gt.mmax) then
!         istep=2*mmax
!         theta=two_pi/(isign*mmax)
!         wpr=-2.0d0*sin(0.5d0*theta)**2
!         wpi=sin(theta)
!         wr=1.0d0
!         wi=0.0d0
!         do 13 m=1,mmax,2
!           do 12 i=m,n,istep
!             j=i+mmax
!             tempr=(wr)*data(j)-(wi)*data(j+1)
!             tempi=(wr)*data(j+1)+(wi)*data(j)
!             data(j)=data(i)-tempr
!             data(j+1)=data(i+1)-tempi
!             data(i)=data(i)+tempr
!             data(i+1)=data(i+1)+tempi
! 12        continue
!           wtemp=wr
!           wr=wr*wpr-wi*wpi+wr
!           wi=wi*wpr+wtemp*wpi+wi
! 13      continue
!         mmax=istep
!       go to 2
!       endif
!       return
!       end subroutine
!
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
!
! subroutine shuff(index,m,n)
!     integer tens,ones
!     dimension index(:)
!     tens=n
!     ones=1
!     j=0
!     k=1
!     do while (j < m)
!        j=j+1
!        if(k.gt. m)then
!          write(*,*)k,m
!          stop
!        endif
!        index(k)=j
!        k=k+tens
!        if(k .gt. m)then
!          ones=ones+1
!          k=ones
!        endif
!     enddo
! end subroutine
!
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
!
!
! SUBROUTINE fourier_transform_rout_(data,nn,isign)
! INTEGER :: isign,nn
! REAL    :: data(2*nn)
! !--------------------------------------------------------------------------------------------!
! !Replaces data(1:2*nn) by its discrete Fourier transform, if isign is input as 1; or replaces
! !data(1:2*nn) by nn times its inverse discrete Fourier transform, if isign is input as -1.
! !data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn
! !MUST be an integer power of 2 (this is not checked for!).
! !--------------------------------------------------------------------------------------------!
! INTEGER i,istep,j,m,mmax,n
! REAL tempi,tempr
! REAL(8) theta,wi,wpi,wpr,wr,wtemp !Double precision for the trigonometric
!                                            !recurrences.
!   n=2*nn
!   j=1
!   do i=1,n,2  !This is the bit-reversal section of the routine.
!     if(j.gt.i)then
!       tempr=data(j)  !Exchange the two complex numbers.
!       tempi=data(j+1)
!       data(j)=data(i)
!       data(j+1)=data(i+1)
!       data(i)=tempr
!       data(i+1)=tempi
!     endif
!     m=nn
! 1   if ((m.ge.2).and.(j.gt.m)) then
!       j=j-m
!       m=m/2
!       goto 1
!     endif
!     j=j+m
!   end do
!   mmax=2      !Here begins the Danielson-Lanczos section of the routine.
! 2 if (n.gt.mmax) then  !Outer loop executed log2 nn times.
!     istep=2*mmax
!     theta=6.28318530717959d0/(isign*mmax) !Initialize for the trigonometric recurrence.
!     wpr=-2.d0*sin(0.5d0*theta)**2
!     wpi=sin(theta)
!     wr=1.d0
!     wi=0.d0
!     do m=1,mmax,2  !Here are the two nested inner loops.
!       do i=m,n,istep
!         j=i+mmax  !This is the Danielson-Lanczos formula:
!         tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
!         tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
!         data(j)=data(i)-tempr
!         data(j+1)=data(i+1)-tempi
!         data(i)=data(i)+tempr
!         data(i+1)=data(i+1)+tempi
!       end do
!       wtemp=wr  !Trigonometric recurrence.
!       wr=wr*wpr-wi*wpi+wr
!       wi=wi*wpr+wtemp*wpi+wi
!     end do
!     mmax=istep
!     goto 2  !Not yet done.
!   endif     !All done.
! return
! END subroutine
!
!
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
!
! !---------------------------------------------------------------------------------------------!
! !  Calculates the Fourier transform of a set of n real-valued data points. Replaces this data !
! !  (which is stored in array data(1:n)) by the positive frequency half of its complex Fourier !
! !  transform. The real-valued first and last components of the complex transform are returned !
! !  as elements data(1) and data(2), respectively. n must be a power of 2. This routine        !
! !  also calculates the inverse transform of a complex data array if it is the                 !
! !  transform of real data. (Result in this case must be multiplied by 2/n.)                   !
! !---------------------------------------------------------------------------------------------!
!
! SUBROUTINE realft__(data,n,isign)
! INTEGER :: isign,n
! REAL    :: data(n)
! INTEGER :: i,i1,i2,i3,i4,n2p3
! REAL    :: c1,c2,h1i,h1r,h2i,h2r,wis,wrs
! REAL(8) :: theta,wi,wpi,wpr,wr,wtemp !Double precision for the trigonometric recurrences.
!
!   theta=pi/dble(n/2)      !Initialize the recurrence.
!
!   c1=0.5
!   if (isign.eq.1) then
!     c2=-0.5
!     call fourier_transform_rout_(data,n/2,+1)                !The forward transform is here.
!   else
!     c2=0.5                                 !Otherwise set up for an inverse transform.
!     theta=-theta
!   endif
!   wpr=-2.0d0*sin(0.5d0*theta)**2
!   wpi=sin(theta)
!   wr=1.0d0+wpr
!   wi=wpi
!   n2p3=n+3
!   do i=2,n/4  !Case i=1 done separately below.
!     i1=2*i-1
!     i2=i1+1
!     i3=n2p3-i2
!     i4=i3+1
!     wrs=sngl(wr)
!     wis=sngl(wi)
!     h1r=c1*(data(i1)+data(i3))  !The two separate transforms are separated out of data.
!     h1i=c1*(data(i2)-data(i4))
!     h2r=-c2*(data(i2)+data(i4))
!     h2i=c2*(data(i1)-data(i3))
!     data(i1)=h1r+wrs*h2r-wis*h2i !Here they are recombined to form the true transform
!                                  !of the original real data.
!     data(i2)=h1i+wrs*h2i+wis*h2r
!     data(i3)=h1r-wrs*h2r+wis*h2i
!     data(i4)=-h1i+wrs*h2i+wis*h2r
!     wtemp=wr                     !The recurrence.
!     wr=wr*wpr-wi*wpi+wr
!     wi=wi*wpr+wtemp*wpi+wi
!   end do
!
!   if (isign.eq.1) then
!     h1r=data(1)
!     data(1)=h1r+data(2)
!     data(2)=h1r-data(2)      !Squeeze the first and last data together to get
!                              !them all within the original array.
!   else
!     h1r=data(1)
!     data(1)=c1*(h1r+data(2))
!     data(2)=c1*(h1r-data(2))
!     call fourier_transform_rout_(data,n/2,-1)  !This is the inverse transform for the case isign=-1.
!   end if
!   return
! END subroutine
!
!
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
!
!       SUBROUTINE realft_(DATA,N,ISIGN)
!       IMPLICIT REAL(8) (A-H,O-Z)
!       REAL(8) WR,WI,WPR,WPI,WTEMP,THETA
!       DIMENSION DATA(*)
!
!       THETA = 2.d0 * pi / 2.0D0/DBLE(N)
!       C1=0.5D0
!       IF (ISIGN.EQ.1) THEN
!         C2=-0.5D0
!         CALL fourier_transform_2(DATA,N,+1)
!       ELSE
!         C2=0.5D0
!         THETA=-THETA
!       ENDIF
!       WPR=-2.0D0*DSIN(0.5D0*THETA)**2
!       WPI=DSIN(THETA)
!       WR=1.0D0+WPR
!       WI=WPI
!       N2P3=2*N+3
!       DO 11 I=2,N/2+1
!         I1=2*I-1
!         I2=I1+1
!         I3=N2P3-I2
!         I4=I3+1
!         H1R=C1*(DATA(I1)+DATA(I3))
!         H1I=C1*(DATA(I2)-DATA(I4))
!         H2R=-C2*(DATA(I2)+DATA(I4))
!         H2I=C2*(DATA(I1)-DATA(I3))
!         DATA(I1)=H1R+WR*H2R-WI*H2I
!         DATA(I2)=H1I+WR*H2I+WI*H2R
!         DATA(I3)=H1R-WR*H2R+WI*H2I
!         DATA(I4)=-H1I+WR*H2I+WI*H2R
!         WTEMP=WR
!         WR=WR*WPR-WI*WPI+WR
!         WI=WI*WPR+WTEMP*WPI+WI
! 11    CONTINUE
!       IF (ISIGN.EQ.1) THEN
!         H1R=DATA(1)
!         DATA(1)=H1R+DATA(2)
!         DATA(2)=H1R-DATA(2)
!       ELSE
!         H1R=DATA(1)
!         DATA(1)=C1*(H1R+DATA(2))
!         DATA(2)=C1*(H1R-DATA(2))
!         CALL fourier_transform_2(DATA,N,-1)
!       ENDIF
!       RETURN
!       END subroutine
!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
!
!       SUBROUTINE fourier_transform_2(DATA,NN,ISIGN)
!       IMPLICIT REAL(8) (A-H,O-Z)
!       REAL(8)  :: WR,WI,WPR,WPI,WTEMP,THETA
!       DIMENSION DATA(*)
!       N=2*NN
!       J=1
!       DO 11 I=1,N,2
!         IF(J.GT.I)THEN
!           TEMPR=DATA(J)
!           TEMPI=DATA(J+1)
!           DATA(J)=DATA(I)
!           DATA(J+1)=DATA(I+1)
!           DATA(I)=TEMPR
!           DATA(I+1)=TEMPI
!         ENDIF
!         M=N/2
! 1       IF ((M.GE.2).AND.(J.GT.M)) THEN
!           J=J-M
!           M=M/2
!         GO TO 1
!         ENDIF
!         J=J+M
! 11    CONTINUE
!       MMAX=2
! 2     IF (N.GT.MMAX) THEN
!         ISTEP=2*MMAX
!         THETA=6.28318530717959D0/(ISIGN*MMAX)
!         WPR=-2.D0*DSIN(0.5D0*THETA)**2
!         WPI=DSIN(THETA)
!         WR=1.D0
!         WI=0.D0
!         DO 13 M=1,MMAX,2
!           DO 12 I=M,N,ISTEP
!             J=I+MMAX
!             TEMPR=WR*DATA(J)-WI*DATA(J+1)
!             TEMPI=WR*DATA(J+1)+WI*DATA(J)
!             DATA(J)=DATA(I)-TEMPR
!             DATA(J+1)=DATA(I+1)-TEMPI
!             DATA(I)=DATA(I)+TEMPR
!             DATA(I+1)=DATA(I+1)+TEMPI
! 12        CONTINUE
!           WTEMP=WR
!           WR=WR*WPR-WI*WPI+WR
!           WI=WI*WPR+WTEMP*WPI+WI
! 13      CONTINUE
!         MMAX=ISTEP
!       GO TO 2
!       ENDIF
!       RETURN
!       END subroutine
!
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************

end module
