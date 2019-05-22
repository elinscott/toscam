module correlations

   use correl_class, only: correl_type
   use genvar, only: DBL
   use green_class, only: green_type
   use readable_vec_class, only: readable_vec_type

   implicit none

   private

   !---------!
   ! PUBLIC  !
   !---------!

   TYPE(correl_type), SAVE, public ::  SNAMBU
   TYPE(correl_type), SAVE, public ::  SNAMBUret
   TYPE(green_type), SAVE, public  ::  GNAMBU, GNAMBUN
   TYPE(green_type), SAVE, public  ::  GNAMBUret
   TYPE(green_type), SAVE, public  ::  G(2), GN(2), GKret(2)
   TYPE(green_type), SAVE, public  ::  GF(2), GNF(2)
   TYPE(green_type), SAVE, public  ::  Spm

   !----------!
   ! PRIVATE  !
   !----------!

   TYPE(green_type), SAVE           :: Gret(2), GFret(2)
   TYPE(green_type), SAVE           :: Sz, N
   TYPE(green_type), SAVE           :: Szret, Spmret, Nret
   TYPE(green_type), SAVE           :: P3, P4
   TYPE(green_type), SAVE           :: P3ret, P4ret
   TYPE(correl_type), SAVE           :: CHI, CHIret

   LOGICAL                            :: SUPER
   REAL(DBL)                          :: Nc, beta, width
   INTEGER, ALLOCATABLE               :: triplets(:, :), quadruplets(:, :)
   TYPE(readable_vec_type), POINTER   :: vec_list(:) => NULL()
   LOGICAL                            :: compute_Sz, compute_Spm, &
                                         compute_N, compute_P3, &
                                         compute_P4
   LOGICAL                            :: write_cor = .true., dump_stat &
                                         = .false.

   public :: compute_correlations
   public :: init_correlations
   public :: write_info_correlations
   public :: write_static

contains

#include "correlations_write.h"
#include "correlations_init.h"

   subroutine compute_correlations(AIM, GS, retard, COMPUTE_ALL_CORREL, &
                                   only_dens)

      use aim_class, only: aim_type
      use apply_c, only: init_apply_c
      use apply_ns, only: init_apply_ns
      use apply_p, only: init_apply_p
      use common_def, only: dump_message, reset_timer, timer_fortran
      use eigen_sector_class, only: eigensectorlist_type
      use genvar, only: bosonic, log_unit, retarded
      use globalvar_ed_solver, only: ALL_FIRST_CALL
      use rcmatrix_class, only: rcmatrix_type
      use rcvector_class, only: rcvector_type

      implicit none

      TYPE(AIM_type), INTENT(IN)             :: AIM
      TYPE(eigensectorlist_type), INTENT(IN) :: GS
      LOGICAL, OPTIONAL, INTENT(IN)          :: COMPUTE_ALL_CORREL
      LOGICAL             :: retard
      INTEGER             :: ii, Nc, start_correl, bndsup(2), bndsdo(2)
      TYPE(rcmatrix_type) :: dmat
      LOGICAL             :: do_all_correl, only_dens
      LOGICAL, SAVE       :: first_call = .true.
      TYPE(rcvector_type) :: vec1, vec2

      do_all_correl = .false.
      IF (PRESENT(COMPUTE_ALL_CORREL)) do_all_correl = COMPUTE_ALL_CORREL

      CALL reset_timer(start_correl)

      Nc = AIM%Nc
      bndsup = (/1, Nc/)
      bndsdo = bndsup + Nc

      CALL dump_message(TEXT="#")
      CALL dump_message(TEXT="#########################################")
      CALL dump_message(TEXT="### BEGIN COMPUTATION OF CORRELATIONS ###")
      CALL dump_message(TEXT="#########################################")
      CALL dump_message(TEXT="#")

      IF (first_call) THEN
         if (.not. ALL_FIRST_CALL) first_call = .false.
         CALL init_apply_C(AIM)
         CALL init_apply_NS(AIM)
         CALL init_apply_P(triplets, quadruplets, AIM)
      ENDIF

      !-----------------------------------------------------------------!
      !-----------------------------------------------------------------!
      write (log_unit, *) ' ==== compute super green ==== '
      call matsubara_super_green

      if (retard .and. .not. only_dens) then
         write (log_unit, *) ' ==== compute retarded functions ==== '
         call compute_retarded_functions
      endif

      IF (do_all_correl .and. .not. only_dens) THEN
         write (log_unit, *) ' ==== compute dynamic bosonic functions ==== '
         call bosonic_dynamic
         write (log_unit, *) ' ==== compute reduced density matrix ==== '
         call reduced_density_matrix_toolbox
      ELSE
         write (log_unit, *) ' ==== compute equal time bosonic ==== '
         if (.not. only_dens) call bosonic_equal_time
      ENDIF

      CALL write_static(AIM, GS)
      !-----------------------------------------------------------------!
      !-----------------------------------------------------------------!

      CALL timer_fortran(start_correl, "### COMPUTING CORRELATIONS TOOK ###")

   contains

#include "correlations_dyn_stat.h"

   end subroutine

   subroutine compute_self_energy(S, G, AIM, retarded)

      use impurity_class, only: nambu_ec
      use globalvar_ed_solver, only: average_g, cutoff_rvb, mask_average
      use masked_matrix_class, only: delete_masked_matrix, masked_matrix_type
      use correl_class, only: average_correlations, correl_type
      use bath_class_hybrid, only: bath2hybrid
      use matrix, only: id, invmat
      use aim_class, only: aim_type

      implicit none

      TYPE(correl_type), INTENT(INOUT) :: S
      TYPE(correl_type), INTENT(IN)    :: G
      TYPE(AIM_type), INTENT(IN)       :: AIM
      TYPE(masked_matrix_type) :: E
      COMPLEX(DBL)             :: Gm1(G%N, G%N, G%Nw)
      COMPLEX(DBL)             :: bath1(G%N, G%N, G%Nw)
      INTEGER                  :: iw, nnn
      logical                  :: retarded, block

      CALL bath2hybrid(AIM%bath, S%freq%title)
      CALL Nambu_Ec(E, AIM%impurity%Ec)

      if (retarded) then
         if (S%Nw /= size(AIM%bath%hybridret%fctn, 3)) then
            write (*, *) 'S%Nw: ', S%Nw
            write (*, *) 'number of frequencies in hybrid : ', &
               size(AIM%bath%hybrid%fctn, 3)
            stop 'error compute_self_energy ED solver'
         endif
      else
         if (S%Nw /= size(AIM%bath%hybrid%fctn, 3)) then
            write (*, *) 'S%Nw: ', S%Nw
            write (*, *) 'number of frequencies in hybrid : ', &
               size(AIM%bath%hybrid%fctn, 3)
            stop 'error compute_self_energy ED solver'
         endif
      endif

      nnn = size(Gm1, 1)
      block = maxval(abs(AIM%BATH%Pb%rc%mat)) < cutoff_rvb .and. &
              maxval(abs(AIM%BATH%PVbc(1)%rc%mat)) < cutoff_rvb .and. &
              maxval(abs(AIM%BATH%PVbc(2)%rc%mat)) < cutoff_rvb

      Gm1 = G%fctn
      if (.not. retarded) then
         bath1 = AIM%bath%hybrid%fctn
      else
         bath1 = AIM%bath%hybridret%fctn
      endif

      if (.not. retarded) then
         if (abs(average_G) >= 3) call average_correlations(AIM%bath%hybrid, &
                                                            bath1, average_G >= 0, MASK_AVERAGE)
      else
         if (abs(average_G) >= 3) call average_correlations(AIM%bath%hybridret, &
                                                            bath1, average_G >= 0, MASK_AVERAGE)
      endif

      if (abs(average_G) >= 2) call average_correlations(G, Gm1, average_G >= &
                                                         0, MASK_AVERAGE)

      DO iw = 1, S%Nw
         CALL invmat(size(Gm1, 1), Gm1(:, :, iw))
      ENDDO

      DO iw = 1, S%Nw
         if (.not. retarded) then
            S%fctn(:, :, iw) = SNAMBU%freq%vec(iw)*Id(nnn) - (E%rc%mat + &
                                                              bath1(:, :, iw) + Gm1(:, :, iw))
         else
            if (aimag(SNAMBUret%freq%vec(iw)) < 0.00001) then
               write (*, *) 'ED, i*delta too small in compute self energy : ', &
                  aimag(SNAMBUret%freq%vec(iw))
               stop
            endif
            S%fctn(:, :, iw) = SNAMBUret%freq%vec(iw)*Id(nnn) - (E%rc%mat + &
                                                                 bath1(:, :, iw) + Gm1(:, :, iw))
         endif
      ENDDO

      CALL delete_masked_matrix(E)

   end subroutine

   subroutine Nambu_green(GNAMBU, G, GF)

      use mask_class, only: delete_mask, mask_type, new_mask
      use correl_class, only: correl2vec
      use genvar, only: dbl
      use masked_matrix_class, only: masked_matrix2vec
      use green_class, only: green_type

      implicit none

      TYPE(green_type), INTENT(INOUT)        :: GNAMBU
      TYPE(green_type), INTENT(IN), OPTIONAL :: G(2)
      TYPE(green_type), INTENT(IN), OPTIONAL :: GF(2)
      INTEGER         :: Nc, bndsup(2), bndsdo(2), iw, iorb, jorb, ipm, mipm
      TYPE(mask_type) :: MASK
      COMPLEX(DBL)    :: swapvec(GNAMBU%Nw)
#ifdef _complex
      COMPLEX(DBL)    :: dswapv
#else
      REAL(DBL)       :: dswapv
#endif

      Nc = G(1)%N
      bndsup = (/1, Nc/)
      bndsdo = bndsup + Nc

      IF (PRESENT(G)) THEN

         !-----------------------------------------------------------------------------!
         ! CAST GREEN'S FUNCTIONS INTO NAMBU FORM I.E.                                 !
         !               ( G1[-,+]  F1[-,-] )                                          !
         ! GNAMBU[-,+] = ( F2[+,+]  G2[+,-] )                                          !
         !                                                                             !
         ! WITH                                                                        !
         ! G1[-,+](i,j) = < c[i,up](z) * C[j,up] >                                     !
         ! G2[+,-](i,j) = < C[i,do](z) * c[j,do] >                                     !
         ! F1[-,-](i,j) = < c[i,up](z) * c[j,do] >                                     !
         ! F2[+,+](i,j) = < C[i,do](z) * C[j,up] > = TRANSPOSE(F1[-,-]) IF H IS REAL   !
         ! AND                                                                         !
         !               ( G1[+,-]  F1[+,+] )                                          !
         ! GNAMBU[+,-] = ( F2[-,-]  G2[-,+] )                                          !
         !                                                                             !
         ! WITH                                                                        !
         ! G1[+,-](i,j) = < C[i,up](z) * c[j,up] >                                     !
         ! G2[-,+](i,j) = < c[i,do](z) * C[j,do] >                                     !
         ! F1[+,+](i,j) = < C[i,up](z) * C[j,do] >                                     !
         ! F2[-,-](i,j) = < c[i,do](z) * c[j,up] > = TRANSPOSE(F1[+,+]) IF H IS REAL   !
         !-----------------------------------------------------------------------------!

         DO ipm = 1, 2
            mipm = 3 - ipm

            GNAMBU%correl(ipm, mipm)%fctn = 0.0_DBL
            GNAMBU%correlstat(ipm, mipm)%rc%mat = 0.0_DBL

            ! NORMAL PART

            ! SPIN UP
            GNAMBU%correl(ipm, mipm)%fctn(bndsup(1):bndsup(2), &
                                          bndsup(1):bndsup(2), :) = G(1)%correl(ipm, mipm)%fctn
            GNAMBU%correlstat(ipm, mipm)%rc%mat(bndsup(1):bndsup(2), &
                                                bndsup(1):bndsup(2)) = G(1)%correlstat(ipm, mipm)%rc%mat

            ! SPIN DOWN
            GNAMBU%correl(ipm, mipm)%fctn(bndsdo(1):bndsdo(2), &
                                          bndsdo(1):bndsdo(2), :) = G(2)%correl(mipm, ipm)%fctn
            GNAMBU%correlstat(ipm, mipm)%rc%mat(bndsdo(1):bndsdo(2), &
                                                bndsdo(1):bndsdo(2)) = G(2)%correlstat(mipm, ipm)%rc%mat

            ! ANOMALOUS PART

            IF (PRESENT(GF)) THEN
               ! UP/DOWN PART
               GNAMBU%correl(ipm, mipm)%fctn(bndsup(1):bndsup(2), &
                                             bndsdo(1):bndsdo(2), :) = GF(1)%correl(ipm, ipm)%fctn
               GNAMBU%correlstat(ipm, mipm)%rc%mat(bndsup(1):bndsup(2), &
                                                   bndsdo(1):bndsdo(2)) = GF(1)%correlstat(ipm, ipm)%rc%mat

               ! DOWN/UP PART

#ifdef _complex
               GNAMBU%correl(ipm, mipm)%fctn(bndsdo(1):bndsdo(2), &
                                             bndsup(1):bndsup(2), :) = GF(2)%correl(mipm, mipm)%fctn
               GNAMBU%correlstat(ipm, mipm)%rc%mat(bndsdo(1):bndsdo(2), &
                                                   bndsup(1):bndsup(2)) = GF(2)%correlstat(mipm, mipm)%rc%mat
#else
               DO iw = 1, GF(1)%Nw
                  GNAMBU%correl(ipm, mipm)%fctn(bndsdo(1):bndsdo(2), &
                                                bndsup(1):bndsup(2), iw) = TRANSPOSE(GF(1)%correl(ipm, &
                                                                                                  ipm)%fctn(:, :, iw))
               ENDDO
               GNAMBU%correlstat(ipm, mipm)%rc%mat(bndsdo(1):bndsdo(2), &
                                                   bndsup(1):bndsup(2)) = TRANSPOSE(GF(1)%correlstat(ipm, &
                                                                                                     ipm)%rc%mat)
#endif

            ENDIF
         ENDDO

      ELSE

         ! WE JUST NEED TO PARTICLE-HOLE TRANSFORM THE MATRIX ELEMENTS IN THE
         ! LOWER-RIGHT CORNER OF GNAMBU THAT WERE NOT DIRECTLY COMPUTED !

         CALL new_mask(MASK, GNAMBU%correl(2, 1)%MM%MASK)

         !-------------------!
         ! NORMAL SPIN DOWN  !
         !-------------------!

         DO iorb = Nc + 1, Nc*2
            DO jorb = Nc + 1, Nc*2
               IF (.NOT. MASK%mat(iorb, jorb)) THEN
                  swapvec = GNAMBU%correl(2, 1)%fctn(iorb, jorb, :)
                  GNAMBU%correl(2, 1)%fctn(iorb, jorb, :) = GNAMBU%correl(1, &
                                                                          2)%fctn(iorb, jorb, :)
                  GNAMBU%correl(1, 2)%fctn(iorb, jorb, :) = swapvec
                  dswapv = GNAMBU%correlstat(2, 1)%rc%mat(iorb, jorb)
                  GNAMBU%correlstat(2, 1)%rc%mat(iorb, jorb) = &
                     GNAMBU%correlstat(1, 2)%rc%mat(iorb, jorb)
                  GNAMBU%correlstat(1, 2)%rc%mat(iorb, jorb) = dswapv
               ENDIF
            ENDDO
         ENDDO

         !--------------------!
         ! ANOMALOUS UP-RIGHT !
         !--------------------!

         DO iorb = 1, Nc
            DO jorb = Nc + 1, Nc*2
               IF (.NOT. MASK%mat(iorb, jorb)) THEN
                  swapvec = GNAMBU%correl(2, 2)%fctn(iorb, jorb, :)
                  GNAMBU%correl(2, 2)%fctn(iorb, jorb, :) = GNAMBU%correl(1, &
                                                                          1)%fctn(iorb, jorb, :)
                  GNAMBU%correl(1, 1)%fctn(iorb, jorb, :) = swapvec
                  dswapv = GNAMBU%correlstat(2, 2)%rc%mat(iorb, jorb)
                  GNAMBU%correlstat(2, 2)%rc%mat(iorb, jorb) = &
                     GNAMBU%correlstat(1, 1)%rc%mat(iorb, jorb)
                  GNAMBU%correlstat(1, 1)%rc%mat(iorb, jorb) = dswapv
               ENDIF
            ENDDO
         ENDDO

         CALL correl2vec(GNAMBU%correl(1, 2))
         CALL correl2vec(GNAMBU%correl(2, 1))
         CALL masked_matrix2vec(GNAMBU%correlstat(1, 2))
         CALL masked_matrix2vec(GNAMBU%correlstat(2, 1))
         CALL delete_mask(MASK)

      ENDIF

   end subroutine

end module
