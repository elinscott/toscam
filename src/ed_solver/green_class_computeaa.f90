module green_class_computeAA

   use genvar, only: DBL

   implicit none

   private

   public :: compute_greenaa

contains

   SUBROUTINE compute_greenAA(green, AIM, beta, GS, Asector, applyA, &
        COMPUTE_DYN, keldysh_level, GS_out)

      use aim_class,                     only: aim_type
      use common_def,                    only: dump_message
      use eigen_class,                   only: eigen_type
      use eigen_sector_class,            only: eigensector_type, &
           eigensectorlist_type
      use genvar,                        only: log_unit, messages3, rank, size2
      use globalvar_ed_solver,           only: donot_compute_holepart, &
           FLAG_MPI_GREENS, USE_TRANSPOSE_TRICK_MPI
      use green_class,                   only: green_type
      use H_class,                       only: delete_H
      use mask_class,                    only: mask_type
      use mpirout,                       only: mpisum
      use sector_class,                  only: sector_type 
      use green_class_compute_symmetric, only: reshuffle_vecAA

      implicit none

      TYPE(green_type),           INTENT(INOUT) :: green
      TYPE(AIM_type),             INTENT(IN)    :: AIM
      REAL(DBL),                  INTENT(IN)    :: beta
      TYPE(eigensectorlist_type), INTENT(IN)    :: GS

      interface 
      ! RETURNS SECTOR OF A,A^+|0> 
         subroutine Asector(Asec,pm,sector)

            use sector_class, only: sector_type 

            implicit none

            TYPE(sector_type), INTENT(INOUT) :: Asec
            TYPE(sector_type), INTENT(IN)    :: sector
            CHARACTER(LEN = 1), INTENT(IN)   :: pm

         end subroutine

         subroutine applyA(Aes, pm, MASK, es, rankr)

            use eigen_sector_class, only: eigensector_type

            implicit none

            TYPE(eigensector_type), INTENT(INOUT) :: Aes
            CHARACTER(LEN = 1), INTENT(IN)        :: pm
            TYPE(eigensector_type), INTENT(IN)    :: es
            INTEGER, INTENT(IN)                   :: rankr
            LOGICAL, INTENT(IN)                   :: MASK(:)

         end subroutine
      end interface

      TYPE(eigensectorlist_type), OPTIONAL   :: GS_out
      LOGICAL, OPTIONAL, INTENT(IN)          :: COMPUTE_DYN
      integer, optional                      :: keldysh_level
      TYPE(sector_type)                      :: Asec
      TYPE(eigensector_type)                 :: Apm_es(2)
      TYPE(eigensector_type), POINTER        :: es    => NULL()
      TYPE(eigen_type),       POINTER        :: eigen => NULL(), OPeigen => &
                                                NULL()
      REAL(DBL)                              :: E0
      INTEGER                                :: applyA_timer, &
                                                compute_dyn_timer, &
                                                compute_green_timer, i1, j1, &
                                                isec_back, i
      INTEGER                                :: ipm, jpm, kpm, iph_max, iph, &
                                                iorb, jorb, iind, jjj1, jjj2, j
      REAL(DBL)                              :: Zpart, boltz
      LOGICAL                                :: already_computed_mean(green%N, 2)
      LOGICAL                                :: compute_dyn_correl
      LOGICAL                                :: ORBMASKvec(2, 2, green%N)
      LOGICAL, SAVE                          :: first_call = .true.
      TYPE(mask_type)                        :: MASK(2, 2)
      INTEGER                                :: isector, ieigen, ii, uup, ddn, &
                                                itot, i_, v(2), issz
      TYPE(eigen_type)                       :: Apmsym
      COMPLEX(DBL)                           :: dyn(green%Nw)
      TYPE(eigen_type), POINTER              :: Ai => NULL()
      LOGICAL                                :: NOT_COMMENSURATE
      INTEGER, allocatable                   :: indices_state_sector(:, :)
      INTEGER                                :: ktot, iiorb, jjorb, kk, jj, &
                                                iorb_f, jorb_f

#ifdef DEBUG
      write(log_unit, '(a)') "DEBUG: entering green_class_computeaa_compute_greenaa"
#endif
      call init_data()
     
      !--------------------------! 
      ! PARSE LIST OF GS SECTORS !
      !--------------------------!
  
      if(FLAG_MPI_GREENS < 2)then
         DO isector = 1, GS%nsector
            call init_sector()
            if(.not.(USE_TRANSPOSE_TRICK_MPI.and.NOT_COMMENSURATE))then
               DO ipm = 1, 2
                  DO jpm = 1, 2
                     IF(GREEN%compute(ipm, jpm))then
                        IF(associated(green%correl(ipm, jpm)%MM%MASK%mat))then
                           IF(ANY(green%correl(ipm, jpm)%MM%MASK%mat))then
                              DO iph = 1, iph_max 
                                 call init_iph()
                                 DO ieigen = 1, GS%es(isector)%lowest%neigen
                                    call loop_over_states          
                                 ENDDO
                                 CALL delete_H()
                              ENDDO 
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            endif
         ENDDO 
      else
         call scan_indices()
         write(*, *) ' === START MAIN LOOP === '
         do kk = rank+1, ktot, size2
            call fix_indices
            if(messages3) write(*, *) 'INIT SECTOR', rank
            call init_sector
            call fix_indices
            if(messages3) write(*, *) 'INIT PH', rank
            call init_iph
            call fix_indices
            if(messages3) write(*, *) 'LOOP OVER STATES', rank
            call loop_over_states
            if(messages3) write(*, *) 'DELETE H', rank
            CALL delete_H()
         enddo
      endif 
  
      if(size2 > 1)then
         if(FLAG_MPI_GREENS > 0)then
            write(*, *) ' === SYNCHRONIZE GREEN FUNCTIONS === '
            do ii = 1, 2
               do jj = 1, 2
                  IF(          GREEN%compute(ii, jj)            )then
                     IF(associated(green%correl(ii, jj)%MM%MASK%mat))then
                        IF(       ANY(green%correl(ii, jj)%MM%MASK%mat))then
                           DO i_ = 1, size(green%correlstat(ii, jj)%rc%vec(:))
                              call mpisum(green%correlstat(ii, jj)%rc%vec(i_))
                              call mpisum(green%correl(ii, jj)%vec(i_, :))
                           enddo
                        endif
                     endif
                  endif
               enddo
            enddo
         endif
      endif
  
      if(donot_compute_holepart)then
         write(*, *) ' --- ASSUME DO NOT COMPUTE HOLE PART --- '
         if(maxval(abs(green%correlstat(1, 2)%rc%vec(:))) < 1.d-15) &
              green%correlstat(1, 2)%rc%vec(:) = 1.d0-green%correlstat(2, 1)%rc%vec(:)
         if(maxval(abs(green%correlstat(2, 1)%rc%vec(:))) < 1.d-15) &
              green%correlstat(2, 1)%rc%vec(:) = 1.d0-green%correlstat(1, 2)%rc%vec(:)
      endif
  
      CALL reshuffle_vecAA(green, compute_dyn_correl)
      call clean_everything()
      write(log_unit, *) 'The partition function is : ', Zpart
  
#ifdef DEBUG
      write(log_unit, '(a)') "DEBUG: leaving green_class_computeaa_compute_greenaa"
#endif

      contains
     
#include "green_class_computeAA.h"

   end subroutine

end module
