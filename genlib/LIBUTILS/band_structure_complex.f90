module band_structure_complex_number

  use genvar
  use tetraedron
  USE com_metrop


contains



!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

SUBROUTINE tetra_weights(z, zekf, weight_irr, itt, max_metropolis_steps, p_min, p_max, ndim, nkp, ntet)
!*****************************************************************************!
!* Calculates k-point weights by Analytic tetrahedron method adapted for     *!
!* complex eigenvalues. Complex eigenvalues can not be simply sorted and     *!
!* therefore bands construction is not obvious. The method uses Metropolis   *!
!* algorithm to sort bands such that for each band distance between energies *!
!* acrros a tetrahedron is minimal.                                          *!
!* The result is stored in weight_irr and only irreducible                   *!
!* k-point weights are stored.                                               *!
!*****************************************************************************!
  IMPLICIT NONE
  !---------- Passed variables ----------
  complex(8), intent(in)  :: z
  complex(8), intent(in)  :: zekf(ndim, nkp)
  complex(8), intent(out) :: weight_irr(ndim, nkp)
  INTEGER, intent(in)     :: itt(5, ntet)
  INTEGER, intent(in)     :: ndim, nkp, ntet, p_min, p_max
  INTEGER, intent(in)     :: max_metropolis_steps
  !---------- Local variables ----------
  complex(8) :: zp(ndim,4), wgh(ndim,4)
  INTEGER    :: identity(ndim,4),indx(ndim,4)
  real(8)     :: dr, mtet
  INTEGER    :: ic, i, p, ikp, itet, rtet, j1, j2, j3
  
  DO ic=1,4
     DO p=1,ndim
        identity(p,ic)=p
     ENDDO
  ENDDO

  weight_irr = 0.0

  mtet = 0
  DO itet = 1,ntet
     mtet = mtet+itt(1,itet)
  ENDDO

  CALL Init_Tmetrop(ndim)

  DO itet = 1,ntet  ! Tetrahedron loop 
     ! Sets the vectors of corner energies to zp(p=1:ndim,ic=1:4)
     DO ic=1,4                        
        ikp = itt(ic+1,itet)
        zp(:,ic) = zekf(:,ikp)
     ENDDO
     indx = identity
     dr = TMetropolis(zp, max_metropolis_steps, ndim)
     IF (dr.GT.1e-6) THEN 
        CALL Permute_TEV_g(ndim, zp, indx)
     ENDIF

     DO p=1,ndim                                   ! Over all eigenvalues
        CALL tetrahedron(z, zp(p,:), wgh(p,:)) ! Actually calculates weights for all 4 corners of one tetrahedra 
        wgh(p,:) = wgh(p,:) / mtet
     ENDDO
     IF (dr.GT.1e-6) THEN
        CALL InversePermute_g(ndim, indx, wgh)
     ENDIF
     DO ic=1,4
        ikp = itt(ic+1,itet)
        ! Stores weights for irreducible k-points
        DO p=p_min,p_max
           weight_irr(p,ikp) = weight_irr(p,ikp) + wgh(p,ic)*itt(1,itet) ! multiplicity of this tetra
        ENDDO
     ENDDO
  ENDDO
  CALL Destroy_metrop
END SUBROUTINE tetra_weights


!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************


SUBROUTINE tetra_weights_int(z, zekf, weight_irr, itt, max_metropolis_steps, ndim, nkp, ntet)
!!! Same subroutine as above exept it computes integral over momentum as well as 
!!! integral over frequency. This is usefull to obtain density in case self-energy is constant
!!! This is the case at very high frequencies in DMFT or in LDA and LDA+U methods.
  IMPLICIT NONE
  !---------- Passed variables ----------
  complex(8), intent(in)  :: z
  complex(8), intent(in)  :: zekf(ndim, nkp)
  complex(8), intent(out) :: weight_irr(ndim, nkp)
  INTEGER, intent(in)     :: itt(5, ntet)
  INTEGER, intent(in)     :: ndim, nkp, ntet
  INTEGER, intent(in)     :: max_metropolis_steps
  !---------- Local variables ----------
  complex(8) :: zp(ndim,4), wgh(ndim,4)
  INTEGER    :: ic, i, p, ikp, itet, rtet, j1, j2, j3
  real(8)     :: mtet
  !---------- Initialize ----------
  weight_irr = 0.0

  mtet = 0
  DO itet = 1,ntet
     mtet = mtet+itt(1,itet)
  ENDDO

  DO itet = 1,ntet  ! Tetrahedron loop 
     ! Sets the vectors of corner energies to zp(p=1:ndim,ic=1:4)
     DO ic=1,4                        
        ikp = itt(ic+1,itet)
        zp(:,ic) = zekf(:,ikp)
     ENDDO
     DO p=1,ndim                                   ! Over all eigenvalues
        CALL tetrahedron_int(z, zp(p,:), wgh(p,:)) ! Actually calculates weights for all 4 corners of one tetrahedra 
        wgh(p,:) = wgh(p,:) / mtet
     ENDDO
     DO ic=1,4
        ikp = itt(ic+1,itet)
        ! Stores weights for irreducible k-points
        DO p=1,ndim
           weight_irr(p,ikp) = weight_irr(p,ikp) + wgh(p,ic)*itt(1,itet) ! multiplicity of this tetra
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE tetra_weights_int

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

SUBROUTINE tetra_zweights_int(om, om1, om2, zekf, weight, weight_int, itt, max_metropolis_steps, ndim, nkp, ntet)
! Before it was called ini_wtsm_dos_and_int
!*****************************************************************************!
!* Calculates k-point weights by Analytic tetrahedron method adapted for     *!
!* complex eigenvalues. Complex eigenvalues can not be simply sorted and     *!
!* therefore bands construction is not obvious. The method uses Metropolis   *!
!* algorithm to sort bands such that for each band distance between energies *!
!* acrros a tetrahedron is minimal.                                          *!
!* The result is stored in weight_irr. Only irreducible weights are stored   *!
!* Calculates both, the usual tetrahedron weights for spectral density       *!
!* as well as integral weights to calcule electron density                   *!
!*****************************************************************************!
  IMPLICIT NONE
  !---------- Passed variables ----------
  complex(8), intent(in)  :: om, om1, om2
  complex(8), intent(in)  :: zekf(ndim, nkp)
  complex(8), intent(out) :: weight(ndim, nkp), weight_int(ndim,nkp)
  INTEGER, intent(in)     :: itt(5, ntet)
  INTEGER, intent(in)     :: ndim, nkp, ntet
  INTEGER, intent(in)     :: max_metropolis_steps
  !---------- Local variables ----------
  complex(8) :: zp(ndim,4), wgh(ndim,4), wgh_int(ndim,4)
  INTEGER    :: identity(ndim,4),indx(ndim,4)
  real(8)     :: dr, mtet
  INTEGER    :: ic, i, p, ikp, itet, rtet, j1, j2, j3
  
  DO ic=1,4
     DO p=1,ndim
        identity(p,ic)=p
     ENDDO
  ENDDO

  weight = 0.0
  weight_int = 0.0

  mtet = 0
  DO itet = 1,ntet
     mtet = mtet+itt(1,itet)
  ENDDO

  CALL Init_Tmetrop(ndim)

  DO itet = 1,ntet  ! Tetrahedron loop 
     ! Sets the vectors of corner energies to zp(p=1:ndim,ic=1:4)
     DO ic=1,4                        
        ikp = itt(ic+1,itet)
        zp(:,ic) = zekf(:,ikp)
     ENDDO
     indx = identity
     ! Finds better reordering of energies
     dr = TMetropolis(zp, max_metropolis_steps, ndim)
     IF (dr.GT.1e-6) THEN 
        CALL Permute_TEV_g(ndim, zp, indx)
     ENDIF

     DO p=1,ndim                                   ! Over all eigenvalues
        CALL ztetrahedron_int(om1, om2, zp(p,:), wgh_int(p,:))
        CALL tetrahedron(om, zp(p,:), wgh(p,:)) ! Actually calculates weights for all 4 corners of one tetrahedra 
        wgh(p,:) = wgh(p,:) / mtet
        wgh_int(p,:) = wgh_int(p,:) / mtet
        ! Checks how bad the weights are
        IF (abs(wgh(p,1))+abs(wgh(p,2))+abs(wgh(p,3))+abs(wgh(p,4)) .GT. 500./mtet) THEN
           print *, p,itet
           write(6,'(A,2g17.8,1x,2g17.8,1x,2g17.8,1x,2g17.8)') &
           & 'Tezave ini_wtsm ', om-zp(p,1), om-zp(p,2), om-zp(p,3), om-zp(p,4)
           write(6,'(A,2f17.8,1x,2f17.8,1x,2f17.8,1x,2f17.8)') & 
           & 'Tezave ini_wtsm ', wgh(p,1)*mtet, wgh(p,2)*mtet, wgh(p,3)*mtet, wgh(p,4)*mtet
           DO i=1,4
              IF (abs(wgh(p,i)).GT.(100./mtet)) wgh(p,i)=0.0
           ENDDO
        ENDIF
     ENDDO
     IF (dr.GT.1e-6) THEN
        CALL InversePermute_g(ndim, indx, wgh)
        CALL InversePermute_g(ndim, indx, wgh_int)
     ENDIF
     DO ic=1,4
        ikp = itt(ic+1,itet)
        ! Stores weights for irreducible k-points
        DO p=1,ndim
           weight(p,ikp)     = weight(p,ikp)     + wgh(p,ic)*itt(1,itet)
           weight_int(p,ikp) = weight_int(p,ikp) + wgh_int(p,ic)*itt(1,itet)
        ENDDO
     ENDDO
  ENDDO
  CALL Destroy_metrop
END SUBROUTINE tetra_zweights_int

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

SUBROUTINE tetra_zweights_int_only(om1, om2, zekf, weight_int, itt, max_metropolis_steps, ndim, nkp, ntet)
! Before it was called ini_wtsm_dos_and_int
!*****************************************************************************!
!* Calculates k-point weights by Analytic tetrahedron method adapted for     *!
!* complex eigenvalues. Complex eigenvalues can not be simply sorted and     *!
!* therefore bands construction is not obvious. The method uses Metropolis   *!
!* algorithm to sort bands such that for each band distance between energies *!
!* acrros a tetrahedron is minimal.                                          *!
!* The result is stored in weight_irr. Only irreducible weights are stored   *!
!* Calculates both, the usual tetrahedron weights for spectral density       *!
!* as well as integral weights to calcule electron density                   *!
!*****************************************************************************!
  IMPLICIT NONE
  !---------- Passed variables ----------
  complex(8), intent(in)  :: om1, om2
  complex(8), intent(in)  :: zekf(ndim, nkp)
  complex(8), intent(out) :: weight_int(ndim,nkp)
  INTEGER, intent(in)     :: itt(5, ntet)
  INTEGER, intent(in)     :: ndim, nkp, ntet
  INTEGER, intent(in)     :: max_metropolis_steps
  !---------- Local variables ----------
  complex(8) :: zp(ndim,4), wgh_int(ndim,4)
  INTEGER    :: identity(ndim,4),indx(ndim,4)
  real(8)     :: dr, mtet
  INTEGER    :: ic, i, p, ikp, itet, rtet, j1, j2, j3
  
  DO ic=1,4
     DO p=1,ndim
        identity(p,ic)=p
     ENDDO
  ENDDO

  weight_int = 0.0

  mtet = 0
  DO itet = 1,ntet
     mtet = mtet+itt(1,itet)
  ENDDO

  CALL Init_Tmetrop(ndim)

  DO itet = 1,ntet  ! Tetrahedron loop 
     ! Sets the vectors of corner energies to zp(p=1:ndim,ic=1:4)
     DO ic=1,4                        
        ikp = itt(ic+1,itet)
        zp(:,ic) = zekf(:,ikp)
     ENDDO
     indx = identity
     ! Finds better reordering of energies
     dr = TMetropolis(zp, max_metropolis_steps, ndim)
     IF (dr.GT.1e-6) THEN 
        CALL Permute_TEV_g(ndim, zp, indx)
     ENDIF

     DO p=1,ndim                                   ! Over all eigenvalues
        CALL ztetrahedron_int(om1, om2, zp(p,:), wgh_int(p,:))
        wgh_int(p,:) = wgh_int(p,:) / mtet
        ! Checks how bad the weights are
     ENDDO
     IF (dr.GT.1e-6) THEN
        CALL InversePermute_g(ndim, indx, wgh_int)
     ENDIF
     DO ic=1,4
        ikp = itt(ic+1,itet)
        ! Stores weights for irreducible k-points
        DO p=1,ndim
           weight_int(p,ikp) = weight_int(p,ikp) + wgh_int(p,ic)*itt(1,itet)
        ENDDO
     ENDDO
  ENDDO
  CALL Destroy_metrop
END SUBROUTINE tetra_zweights_int_only

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

SUBROUTINE cmp_simple_integral_weights(weight, w1,w2,zekf,wk,ndim,nkp)
  IMPLICIT NONE
  complex(8), intent(out) :: weight(ndim,nkp)
  real(8), intent(in)     :: w1, w2
  complex(8), intent(in)  :: zekf(ndim, nkp)
  real(8), intent(in)     :: wk(nkp)
  INTEGER, intent(in)     :: ndim, nkp
  INTEGER                 :: ikp, p
  DO ikp=1,nkp
     DO p=1,ndim
        weight(p,ikp) = wk(ikp)*(log(cmplx(w2,1e-10,8)-zekf(p,ikp))-log(cmplx(w1,1e-10,8)-zekf(p,ikp)))
     ENDDO
  ENDDO
END SUBROUTINE cmp_simple_integral_weights

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************


end module
