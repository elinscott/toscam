MODULE green_class_computeAB

   use genvar, only: DP

   IMPLICIT NONE

   private

   public :: compute_greenab

contains

   SUBROUTINE compute_greenAB(greenAB, greenBA, greenAA, greenBB, AIM, beta, &
                              GS, Asector, Bsector, applyA, applyB, COMPUTE_DYN)

      use aim_class, only: aim_type
      use eigen_class, only: eigen_type
      use eigen_sector_class, only: eigensector_type, &
         eigensectorlist_type, gsenergy, partition
      use genvar, only: messages3, rank, size2
      use globalvar_ed_solver, only: FLAG_MPI_GREENS, &
         USE_TRANSPOSE_TRICK_MPI
      use green_class, only: green_type
      use mask_class, only: mask_type
      use mpi_mod, only: mpisum
      use sector_class, only: sector_type
      use green_class_compute_symmetric, only: reshuffle_vecAB

      implicit none

      TYPE(green_type), INTENT(INOUT) :: greenAB
      TYPE(green_type), INTENT(INOUT) :: greenBA
      TYPE(green_type), INTENT(IN)    :: greenAA
      TYPE(green_type), INTENT(IN)    :: greenBB
      TYPE(AIM_type), INTENT(IN)    :: AIM
      REAL(DP), INTENT(IN)    :: beta
      TYPE(eigensectorlist_type), INTENT(IN)    :: GS
      INTEGER                                   :: ii, isec_back

      INTERFACE
         ! RETURNS SECTOR OF A,A^+|0 >
         subroutine Asector(Asec, pm, sector)

            use sector_class, only: sector_type

            implicit none

            TYPE(sector_type), INTENT(INOUT) :: Asec
            TYPE(sector_type), INTENT(IN)    :: sector
            CHARACTER(LEN=1), INTENT(IN)   :: pm

         end subroutine

         ! RETURNS SECTOR OF B,B^+|0 >
         subroutine Bsector(Bsec, pm, sector)

            use sector_class, only: sector_type

            implicit none

            TYPE(sector_type), INTENT(INOUT) :: Bsec
            TYPE(sector_type), INTENT(IN)    :: sector
            CHARACTER(LEN=1), INTENT(IN)   :: pm

         end subroutine

         ! COMPUTES A,A^+|0 >
         subroutine applyA(Aes, pm, MASK, es, rankr)

            use eigen_sector_class, only: eigensector_type

            implicit none

            TYPE(eigensector_type), INTENT(INOUT) :: Aes
            CHARACTER(LEN=1), INTENT(IN)        :: pm
            TYPE(eigensector_type), INTENT(IN)    :: es
            INTEGER, INTENT(IN)                   :: rankr
            LOGICAL, INTENT(IN)                   :: MASK(:)

         end subroutine

         ! COMPUTES B,B^+|0 >
         subroutine applyB(Bes, pm, MASK, es, rankr)

            use eigen_sector_class, only: eigensector_type

            implicit none

            TYPE(eigensector_type), INTENT(INOUT) :: Bes
            CHARACTER(LEN=1), INTENT(IN)        :: pm
            TYPE(eigensector_type), INTENT(IN)    :: es
            INTEGER, INTENT(IN)                   :: rankr
            LOGICAL, INTENT(IN)                   :: MASK(:)
         END SUBROUTINE
      END INTERFACE

      LOGICAL, OPTIONAL, INTENT(IN)    :: COMPUTE_DYN
      TYPE(sector_type)                         :: Asec, Bsec
      TYPE(eigensector_type)                    :: Apm_es(2), Bpm_es(2)
      TYPE(eigensector_type), POINTER           :: es => NULL()
      TYPE(eigen_type), POINTER           :: eigen => NULL(), OPeigen => &
                                             NULL()
      REAL(DP)                                 :: E0
      INTEGER                                   :: apply_timer, &
                                                   compute_dyn_timer, &
                                                   compute_green_timer, i
      INTEGER                                   :: ipm, jpm, kpm, iph, iorb, &
                                                   jorb, iind, iw, i1, j1, j, &
                                                   i_, v(2), issz
      REAL(DP)                                 :: Zpart, boltz
      COMPLEX(DP)                              :: dyn(greenAB%Nw)
      LOGICAL                                   :: BA, compute_dyn_correl
      LOGICAL                                   :: ORBMASKvec(2, 2, greenAB%N)
      LOGICAL, SAVE                             :: first_call = .true.
      TYPE(mask_type)                           :: MASKAB(2, 2), MASKBA(2, 2)
      INTEGER                                   :: ieigen, isector, uup, ddn, &
                                                   itot, ssz, iiorb, jjorb, &
                                                   ktot, kk, jj, iieigen
      TYPE(eigen_type)                          :: ABpmsym
      TYPE(eigen_type), POINTER                 :: Ai => NULL()
      LOGICAL                                   :: NOT_COMMENSURATE
      INTEGER, allocatable                      :: indices_state_sector(:, :)
      INTEGER                                   :: iorb_f, jorb_f

      itot = AIM%bath%Nb + AIM%impurity%Nc

      compute_dyn_correl = .true.
      IF (PRESENT(COMPUTE_DYN)) compute_dyn_correl = COMPUTE_DYN

#ifdef _complex
      BA = .true.
#else
      BA = .false.
#endif

      call init_data()

      Zpart = partition(beta, GS)
      E0 = GSenergy(GS)

      write (*, *) ' --- PARSE LIST OF GS SECTORS FOR GREEN FUNCTION --- '
      write (*, *) ' ... compute AB green functions : ', greenAB%compute(1, 1:2), greenAB%compute(2, 1:2)

      if (FLAG_MPI_GREENS < 2) then
         DO isector = 1, GS%nsector
            write (*, *) 'initialize sector'
            call init_sector()
            write (*, *) '..done..'
            if (.not. (USE_TRANSPOSE_TRICK_MPI .and. NOT_COMMENSURATE)) then
               DO ipm = 1, 2
                  DO jpm = 1, 2
                     IF (greenAB%compute(ipm, jpm)) then
                        IF (associated(greenAB%correl(ipm, jpm)%MM%MASK%mat)) then
                           IF (ANY(greenAB%correl(ipm, jpm)%MM%MASK%mat)) then
                              write (*, *) 'COMPUTING ipm,jpm : ', ipm, jpm
                              DO iph = 1, 2
                                 write (*, *) 'particle-hole contribution', iph
                                 call print_label()
                                 call do_ab()
                                 call do_ba()
                              ENDDO
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      else
         call scan_indices
         write (*, *) ' === START MAIN LOOP === '
         do kk = rank + 1, ktot, size2
            call fix_indices()
            if (messages3) write (*, *) 'INIT SECTOR', rank
            call init_sector()
            call fix_indices()
            call print_label()
            call do_ab()
            call do_ba()
         enddo
      endif

      if (size2 > 1) then
         if (FLAG_MPI_GREENS > 0) then
            write (*, *) ' === SYNCHRONIZE GREEN FUNCTIONS === '
            do ii = 1, 2
               do jj = 1, 2
                  IF (GREENAB%compute(ii, jj)) then
                     IF (associated(greenAB%correl(ii, jj)%MM%MASK%mat)) then
                        IF (ANY(greenAB%correl(ii, jj)%MM%MASK%mat)) then
                           DO i_ = 1, size(greenAB%correlstat(ii, jj)%rc%vec(:))
                              call mpisum(greenAB%correlstat(ii, jj)%rc%vec(i_))
                              call mpisum(greenAB%correl(ii, jj)%vec(i_, :))
                              if (BA) then
                                 call mpisum(greenBA%correlstat(ii, jj)%rc%vec(i_))
                                 call mpisum(greenBA%correl(ii, jj)%vec(i_, :))
                              endif
                           enddo
                        endif
                     endif
                  endif
               enddo
            enddo
         endif
      endif

      CALL reshuffle_vecAB(greenAB, greenBA, greenAA, greenBB, &
                           compute_dyn_correl, BA)
      call clean_everything()

   contains

#include "green_class_computeAB.h"

   end subroutine

end module
