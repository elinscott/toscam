MODULE overlap_module

   implicit none

   private

   !---------------------------------------------------------------------!
   ! COMPUTE THE WEIGHT OF THE GROUND STATE ON CERTAIN REFERENCE VECTORS !
   !---------------------------------------------------------------------!

   public :: compute_overlap

contains

   subroutine compute_overlap(AIM, beta, GS, vec_list)

      ! COMPUTE STATIC CORRELATIONS

      use eigen_sector_class, only: eigensectorlist_type, gsenergy, nlowest, &
         partition
      use linalg, only: conj, dexpc, MPLX
      use genvar, only: dp, log_unit, pi
      use aim_class, only: aim2impbathstate, aim_type
      use common_def, only: dump_message, find_rank
      use sector_class, only: dimen_func, state_func
      use eigen_class, only: eigen_type
      use readable_vec_class, only: readable_vec_type

      implicit none

      TYPE(AIM_type), INTENT(IN)             :: AIM
      REAL(DP), INTENT(IN)                  :: beta
      TYPE(eigensectorlist_type), INTENT(IN) :: GS
      TYPE(readable_vec_type), INTENT(INOUT) :: vec_list(:)
      TYPE(eigen_type), POINTER :: eigen => NULL()
      INTEGER                   :: isector, ieigen, iGS, nGS, ivec, iket, &
                                   AIMrank, IMPstate, BATHstate
      REAL(DP)                 :: boltz, Zpart, E0
      COMPLEX(DP)              :: coeff
      COMPLEX(DP), ALLOCATABLE :: overlap(:, :), roverlap(:, :)
      REAL(DP), ALLOCATABLE    :: phase_(:)
      INTEGER                   :: nvec, nkets, kets(100)

      CALL dump_message(TEXT="######################################")
      CALL dump_message(TEXT="### OVERLAP WITH REFERENCE VECTORS ###")
      CALL dump_message(TEXT="######################################")

      !--------------------------------------------!
      ! FIRST WE FIND ALL DISTINCT IMPURITY STATES !
      !--------------------------------------------!

      nvec = SIZE(vec_list)
      nkets = 0
      kets = -1

      DO ivec = 1, nvec
         DO iket = 1, vec_list(ivec)%nket
            IF (find_rank(vec_list(ivec)%state(iket), kets(1:nkets)) == 0) THEN
               nkets = nkets + 1
               kets(nkets) = vec_list(ivec)%state(iket)
            ENDIF
         ENDDO
      ENDDO

      !----------------------------------------------------!
      ! FINALLY COMPUTE THE OVERLAP WITH REFERENCE VECTORS !
      !----------------------------------------------------!

      Zpart = partition(beta, GS)
      E0 = GSenergy(GS)
      nGS = nlowest(GS)

      ALLOCATE (overlap(nvec, nGS))
      overlap = 0.0_DP
      iGS = 0

      DO isector = 1, GS%nsector
         DO ieigen = 1, GS%es(isector)%lowest%neigen
            iGS = iGS + 1
            eigen => GS%es(isector)%lowest%eigen(ieigen)
            boltz = DEXPc(-beta*(eigen%val - E0))/Zpart
            DO AIMrank = 1, dimen_func(GS%es(isector)%sector)
               CALL AIM2IMPBATHstate(IMPstate, BATHstate, state_func(AIMrank, &
                                                                     GS%es(isector)%sector), AIM%Nc, AIM%Nb)
               IF (find_rank(IMPstate, kets) /= 0) THEN
                  DO ivec = 1, nvec
                     DO iket = 1, vec_list(ivec)%nket
                        IF (vec_list(ivec)%state(iket) == IMPstate) THEN
                           coeff = vec_list(ivec)%coeff(iket)
                           overlap(ivec, iGS) = overlap(ivec, iGS) + &
                                                conj(coeff)*eigen%vec%rc(AIMrank)*boltz
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      !-----------------------------------------------------------------------!
      ! RENORMALIZE WITH OVERALL PHASE FACTOR SUCH THAT FIRST OVERLAP IS REAL !
      !-----------------------------------------------------------------------!

      IF (nvec >= 1) THEN
         ALLOCATE (phase_(nGS), roverlap(nvec, nGS))
         DO iGS = 1, nGS
            phase_(iGS) = ATAN2(AIMAG(overlap(1, iGS)), &
               real(overlap(1, iGS), kind=DP))
            roverlap(:, iGS) = overlap(:, iGS)*MPLX(-phase_(iGS))
         ENDDO
      ENDIF

      DO ivec = 1, nvec
         write (log_unit, *) &
            '----------------------------------------------------'
         write (log_unit, *) 'projection on state : ', ivec, &
            TRIM(ADJUSTL(vec_list(ivec)%title))
         write (log_unit, *) 'amplitude in each sectors   : '
         write (log_unit, '(1000f8.5)') (abs(overlap(ivec, iGS)), iGS=1, nGS)
         write (log_unit, *) 'phase in each sectors (pi)  : '
         write (log_unit, '(1000f8.5)') (ATAN2(AIMAG(roverlap(ivec, iGS)), &
               real(roverlap(ivec, iGS), kind=DP))/pi, iGS=1, nGS)
         write (log_unit, *) &
            '----------------------------------------------------'
      ENDDO

      DEALLOCATE (overlap)
      if (allocated(phase_)) DEALLOCATE (phase_, roverlap)

   end subroutine

end module
