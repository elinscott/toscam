module dmft_solver_ed

   use AIM_class, only: AIM_type
   use bath_class, only: bath_type
   use impurity_class, only: impurity_type
   use eigen_sector_class, only: eigensectorlist_type

   implicit none

   private

   CHARACTER(LEN=100) :: BATHfile
   CHARACTER(LEN=100) :: CORRELfile
   CHARACTER(LEN=100) :: EDfile

   LOGICAL              :: start_from_old_gs
   CHARACTER(LEN=100) :: OLDGSFILE

   TYPE(impurity_type), SAVE :: impurity
   TYPE(bath_type), SAVE :: bath
   TYPE(AIM_type), SAVE :: AIM
   TYPE(eigensectorlist_type), SAVE :: GS

   TYPE h1occup
      LOGICAL :: ifdiag, run
      INTEGER :: ndeg
      INTEGER :: n
      COMPLEX(KIND=8), pointer :: Hn(:, :)
      REAL(KIND=8), pointer :: En(:)
      INTEGER, pointer  :: st_n(:)
      INTEGER, pointer  :: narr(:)
   ENDTYPE h1occup

   public :: stand_alone_ed

contains

#include "dmft_ed_solver_read_write.h"
#include "dmft_ed_solver_hubbardI.h"

   subroutine solver_ED_interface(nbathparam_in, iter_, mmu, g_out, self_out, &
                                  hybrid_in, Eimp, sigw, rdens, impurity_, Himp, gw, frequ, &
                                  flip_input_output, para_state_, retarded, restarted_, compute_all, &
                                  tot_rep, spm_, corhop_, only_density, skip_fit, &
                                  use_specific_set_parameters_, param_input_, param_output_, &
                                  nbathparam_, freeze_poles_delta_, freeze_poles_delta_iter_, &
                                  freeze_poles_delta_niter_, imp_causality, Jhund_, &
                                  use_input_delta_instead_of_fit, Jhund_matrix)

      use mesh, only: mirror_array
      use splines, only: resampleit
      use globalvar_ed_solver, only: average_g, energy_global_shift, &
         energy_global_shift2, first_iter_use_edinput, flag_slater_int, &
         fmos, fmos_hub1, freeze_poles_delta, freeze_poles_delta_iter, &
         freeze_poles_delta_niter, iterdmft, jhund, jjmatrix, mask_average, &
         mask_average_sigma, ncpt_approx, ncpt_approx_tot, &
         ncpt_flag_two_step_fit, ncpt_para, ncpt_tot, only_compute_density, &
         para_state, param_input, param_output, restarted, skip_fit_, &
         tot_repulsion, use_specific_set_parameters
      use impurity_class, only: hamiltonian, update_impurity, web
      use genvar, only: fermionic, log_unit, no_mpi, pi, rank, &
         size2, dp
      use correl_class, only: average_correlations
      use matrix, only: average_vec, diag, write_array
      use bath_class_hybrid, only: bath2hybrid, hybrid2bath
      use mpi_mod, only: mpibarrier
      use solver, only: solve_aim
      use correlations, only: compute_correlations, G, GF, GNAMBU, &
         GNAMBUN, GNAMBUret, SNAMBU, SNAMBUret, Spm
      use linalg, only: modi

      implicit none

      complex(8)        :: g_out(:, :, :), gw(:, :, :), self_out(:, :, :), &
                           hybrid_in(:, :, :), sigw(:, :, :), sigww(size(sigw, 1), &
                                                                    size(sigw, 2), size(sigw, 3))
      complex(8)        :: gw_(size(gw, 1), size(gw, 2), size(gw, 3)), &
                           sigw_(size(sigw, 1), size(sigw, 2), size(sigw, 3))
      integer           :: i, j, a1, a2, isector, iter_, nbathparam_in, &
                           ncpt_tot_back, ncpt_approx_tot_back
      integer, optional :: nbathparam_
      real(8), optional :: param_input_(nbathparam_in), &
                           param_output_(nbathparam_in), Jhund_
      complex(8)        :: Eimp(:, :), spm_(:, :, :), corhop_(:, :, :)
      real(8), optional :: Jhund_matrix(size(Eimp, 1)/2, size(Eimp, 1)/2)
      real(8)           :: rdens(:), mmu, frequ(:), tot_rep
      type(hamiltonian) :: Himp
      type(web)         :: impurity_
      logical           :: flip_input_output, para_state_, restarted_, &
                           retarded, compute_all
      logical, optional :: only_density, skip_fit, &
                           use_specific_set_parameters_, freeze_poles_delta_, &
                           imp_causality
      logical, optional :: use_input_delta_instead_of_fit
      integer, optional :: freeze_poles_delta_iter_, freeze_poles_delta_niter_
      real(8)           :: bandes(size(gw, 1), size(frequ))

      if (2*impurity_%N /= size(g_out, 1)) then
         write (*, *) 'size of impurity does not match size of impurity green &
              &functions'
         write (*, *) 'error, critical'
         stop
      endif

      Jhund = 0.d0
      if (present(Jhund_)) Jhund = Jhund_

      if (.not. allocated(JJmatrix)) then
         write (*, *) 'JJ matrix not allocated, it should be if ED solver &
              &initialized'
         stop
      endif

      JJmatrix = Jhund

      if (present(Jhund_matrix)) then

         JJmatrix = Jhund_matrix

         Jhund = maxval(abs(Jhund_matrix))
      endif

      if (flag_slater_int) then
         JJmatrix = 0.d0
         Jhund = 0.d0
      endif

      if (fmos_hub1) then
         call hub1_run
         return
      endif

      freeze_poles_delta_niter = 1
      if (present(freeze_poles_delta_niter_)) freeze_poles_delta_niter = &
         freeze_poles_delta_niter_

      freeze_poles_delta_iter = 1
      if (present(freeze_poles_delta_iter_)) freeze_poles_delta_iter = &
         freeze_poles_delta_iter_

      freeze_poles_delta = .false.
      if (present(freeze_poles_delta_)) freeze_poles_delta = &
         freeze_poles_delta_

      only_compute_density = .false.
      if (present(only_density)) only_compute_density = only_density

      use_specific_set_parameters = .false.
      if (present(use_specific_set_parameters_)) use_specific_set_parameters = &
         use_specific_set_parameters_

      skip_fit_ = .false.
      if (present(skip_fit)) skip_fit_ = skip_fit

      if (bath%nparam == 0) stop 'error ed solver 0 parameters in the bath'
      if (present(nbathparam_)) nbathparam_ = bath%nparam + ncpt_para*ncpt_tot
      if (.not. allocated(param_output)) allocate (param_output(bath%nparam + &
                                                                ncpt_para*ncpt_tot))
      if (.not. allocated(param_input)) allocate (param_input(bath%nparam + &
                                                              ncpt_para*ncpt_tot))
      if (present(param_input_)) then
         param_input = param_input_
         where (abs(param_input) > 1.d10) param_input = 0.
         if (use_specific_set_parameters) then
            write (*, *) &
               '---------------------------------------------------------'
            write (*, *) 'ED SOLVER, USE INPUT PARAMETER SPECIFIED : ', &
               param_input
            write (*, *) &
               '---------------------------------------------------------'
         endif
      endif

      iterdmft = iter_
      restarted = restarted_

      ! --- > FLIP INPUT

      Eimp = flip_matrix_Eimp(Eimp, flip_input_output)
      do i = 1, size(hybrid_in, 3)
         hybrid_in(:, :, i) = flip_matrix(hybrid_in(:, :, i), flip_input_output)
      enddo

      para_state = para_state_
      energy_global_shift2 = 0.d0
      energy_global_shift = 0.d0

      a1 = size(sigw, 1)
      a2 = size(sigw, 2)

      if (size2 > 1 .and. .not. no_mpi) then
         write (*, *) 'UPDATING IMPURITY, RANK = ', rank
         call mpibarrier
      endif

#ifdef _complex
      call update_impurity(impurity_%N, impurity, mmu, impurity_, Himp, Eimp= &
                           Eimp)
#else
      call update_impurity(impurity_%N, impurity, mmu, impurity_, Himp, Eimp= &
                           real(Eimp))
#endif

      if (size2 > 1 .and. .not. no_mpi) then
         write (*, *) 'FITTING BATH,    RANK = ', rank
         write (*, *) 'DMFT ITERATION / RANK = ', iter_, rank
         call mpibarrier
      endif

      if (.not. (iter_ == 1 .and. (first_iter_use_edinput .or. &
                                   start_from_old_gs))) then
         bath%hybrid%fctn(:, :, :) = hybrid_in
         ncpt_tot_back = ncpt_tot
         if (ncpt_tot_back /= 0 .and. ncpt_flag_two_step_fit) then
            ncpt_tot = 0
            ncpt_approx_tot_back = ncpt_approx_tot
            ncpt_approx_tot = 0
         endif
         CALL hybrid2bath(bath)
         if (ncpt_tot_back /= 0 .and. ncpt_flag_two_step_fit) then
            bath%hybrid%fctn(:, :, :) = hybrid_in
            ncpt_tot = ncpt_tot_back
            ncpt_approx_tot = ncpt_approx_tot_back
            CALL hybrid2bath(bath)
         endif
      else
         CALL bath2hybrid(bath, FERMIONIC)
         if (rank == 0) then
            !call PGSUBP(4, 4)
            do j = 1, size(bath%hybrid%fctn(:, 1, 1))
               !call plotarray(aimag(bath%hybrid%freq%vec), real
               ! (bath%hybrid%fctn(j, j, :)), 'ED re bath start iter1')
               !call plotarray(aimag(bath%hybrid%freq%vec),
               ! aimag(bath%hybrid%fctn(j, j, :)), 'ED im bath start iter1')
            enddo
            !call PGSUBP(1, 1)
         endif
      endif

      if (fmos) then
         call fmos_solver
         if (present(param_output_)) param_output_ = param_output
         return
      endif

      if (size2 > 1 .and. .not. no_mpi) then
         write (*, *) 'SOLVING AIM, RANK = ', rank
         call mpibarrier
      endif

      CALL solve_AIM(GS, AIM, OLDGSFILE=OLDGSFILE, COMPUTE_ALL_CORREL= &
                     compute_all)
      write (log_unit, *) &
         '##########################################################'
      do isector = 1, GS%nsector
         write (log_unit, *) '------ > PARSE EIGEN SECTORS : ', isector, &
            GS%nsector
         if (associated(GS%es(isector)%lowest%eigen)) write (log_unit, *) &
            'eigenvalues ...... ', GS%es(isector)%lowest%eigen(:)%val
      enddo

      if (size2 > 1 .and. .not. no_mpi) then
         write (*, *) 'COMPUTING CORRELATIONS, RANK = ', rank
         call mpibarrier
      endif

      write (log_unit, *) ' ================================ '
      write (log_unit, *) '..... compute correlations .....'
      CALL compute_correlations(AIM, GS, retarded, COMPUTE_ALL_CORREL= &
                                compute_all, only_dens=only_compute_density)
      write (log_unit, *) '..... correlations done ........'
      write (log_unit, *) ' ================================ '

      write (log_unit, *) ' ==== compute keldysh ==== '
      call keldysh
      write (log_unit, *) ' ==== keldysh done ==== '

      if (only_compute_density) goto 144

      if (present(use_input_delta_instead_of_fit)) then
         if (use_input_delta_instead_of_fit) then
            SNAMBU%fctn = SNAMBU%fctn + hybrid_in - bath%hybrid%fctn
         endif
      endif

      if (present(imp_causality)) then
         if (imp_causality) then
            do j = 1, size(self_out, 1)
               do i = 1, size(self_out, 3)
                  if (aimag(SNAMBU%fctn(j, j, i)) > 0.) SNAMBU%fctn(j, j, i) = &
                     CMPLX(real(SNAMBU%fctn(j, j, i)), 0.d0)
               enddo
               do i = 1, size(sigww, 3)
                  if (aimag(SNAMBUret%fctn(j, j, i)) > 0.) SNAMBUret%fctn(j, j, &
                                                                          i) = CMPLX(real(SNAMBUret%fctn(j, j, i)), 0.d0)
               enddo
            enddo
         endif
      endif

      if (size2 > 1 .and. .not. no_mpi) then
         write (*, *) 'PREPARING OUTPUTS, RANK = ', rank
         call mpibarrier
      endif

      g_out = GNAMBU%correl(2, 1)%fctn
      self_out = SNAMBU%fctn
      sigww = SNAMBUret%fctn
      gw = GNAMBUret%correl(2, 1)%fctn
      if (ncpt_approx /= 0) then
         call cpt_extract_G(GNAMBUN%correl(2, 1)%freq%vec, self_out, g_out, &
                            imaginary=.true.)
         call cpt_extract_G(GNAMBUret%correl(2, 1)%freq%vec, sigww, gw, &
                            imaginary=.false.)
      endif

      write (*, *) 'shape Spm cor', shape(Spm%correl(1, 2)%fctn)
      write (*, *) 'shape corhop', shape(GNAMBUN%correl(2, 1)%fctn)
      spm_(1, 1, :) = Spm%correl(1, 2)%fctn(1, 1, :)
      corhop_ = GNAMBUN%correl(2, 1)%fctn

      !call plotarray(aimag(GNAMBUN%correl(2, 1)%freq%vec),
      ! real(GNAMBUN%correl(2, 1)%fctn(1, 1, :)),
      ! 'inside_ed_solver_GnambuN_real')
      !call plotarray(aimag(GNAMBUN%correl(2, 1)%freq%vec),
      ! aimag(GNAMBUN%correl(2, 1)%fctn(1, 1, :)),
      ! 'inside_ed_solver_GnambuN_aimag')

      do j = 1, size(sigww, 1)
         bandes(j, :) = aimag(sigww(j, j, :))
      enddo
      !call plotarray(real(GNAMBUret%correl(2, 1)%freq%vec), bandes,
      ! 'inside_ed_solver_self_energy_bef_average')

      if (abs(average_G) == 1) call average_correlations(SNAMBU, self_out, &
                                                         average_G >= 0, MASK_AVERAGE_SIGMA)
      if (abs(average_G) == 1) call average_correlations(SNAMBUret, sigww, &
                                                         average_G >= 0, MASK_AVERAGE_SIGMA)

      do j = 1, size(sigww, 1)
         bandes(j, :) = aimag(sigww(j, j, :))
      enddo

      !call plotarray(real(GNAMBUret%correl(2, 1)%freq%vec), bandes,
      ! 'inside_ed_solver_self_energy_aft_average' //
      ! ADJUSTL(TRIM(SNAMBU%stat)))
      if (abs(average_G) == 1) call average_correlations(GNAMBU%correl(2, 1), &
                                                         g_out, average_G >= 0, MASK_AVERAGE_SIGMA)
      if (abs(average_G) == 1) call average_correlations(GNAMBUN%correl(2, 1), &
                                                         corhop_, average_G >= 0, MASK_AVERAGE_SIGMA)

      ! --- > FLIP OUTPUT

      do i = 1, size(g_out, 3)
         self_out(:, :, i) = flip_matrix(self_out(:, :, i), flip_input_output)
         g_out(:, :, i) = flip_matrix(g_out(:, :, i), flip_input_output)
         corhop_(:, :, i) = flip_matrix(corhop_(:, :, i), flip_input_output)
      enddo

      do i = 1, a1
         do j = 1, a2
            if (.not. flip_input_output .or. .not. (i > a1/2 .and. j > a2/2) &
                ) then
               call resampleit(real(GNAMBUret%correl(2, 1)%freq%vec), gw(i, j, &
                                                                         :), frequ, gw_(i, j, :), 0.d0)
               call resampleit(real(GNAMBUret%correl(2, 1)%freq%vec), sigww(i, &
                                                                            j, :), frequ, sigw_(i, j, :), 0.d0)
               gw(i, j, :) = gw_(i, j, :)
               sigw(i, j, :) = sigw_(i, j, :)
            else
               call mirror_array(real(GNAMBUret%correl(2, 1)%freq%vec), gw(i, &
                                                                           j, :), frequ, gw_(i, j, :))
               call mirror_array(real(GNAMBUret%correl(2, 1)%freq%vec), &
                                 sigww(i, j, :), frequ, sigw_(i, j, :))
               gw(i, j, :) = -conjg(gw_(i, j, :))
               sigw(i, j, :) = -conjg(sigw_(i, j, :))
            endif
         enddo
      enddo

144   continue
      DO i = 1, size(rdens)
         if (i <= size(rdens)/2) then
            rdens(i) = real(G(1)%correlstat(1, 2)%rc%mat(i, i), kind=DP)
         else
            j = modi(i, size(rdens)/2)
            rdens(i) = real(G(2)%correlstat(1, 2)%rc%mat(j, j), kind=DP)
         endif
      ENDDO

      if (abs(average_G) > 0) then
         write (145 + rank, *) 'DENSITY BEFORE AVERAGE OVER MASK : ', rdens
         call average_vec(rdens, diag(MASK_AVERAGE))
         write (145 + rank, *) 'AVERAGE DENSITY : ', rdens
         write (145 + rank, *) 'MASK FOR AVERAGE : ', diag(MASK_AVERAGE)
      endif

      call write_array(AIM%BATH%hybrid%fctn(:, :, i), ' hybrid ', unit= &
                       30000)

      !call dump_to_fortb

      tot_rep = tot_repulsion

      if (present(param_output_)) param_output_ = param_output

      do j = 1, size(bath%hybrid%fctn(:, 1, 1))
         bandes(j, :) = -aimag(gw(j, j, :))/pi
      enddo
      if (.not. flip_input_output) then
         do j = a1/2 + 1, a1
            call mirror_array(frequ, bandes(j, :))
         enddo
      endif
      !call plotarray(frequ, bandes, 'ED_IMP_DOS')

      open (unit=87120, file='chiloc_vertex_sigma_ret_spin_up')
      open (unit=87121, file='chiloc_vertex_sigma_matsu_spin_up')
      open (unit=87122, file='chiloc_vertex_green_ret_spin_up')
      open (unit=87123, file='chiloc_vertex_green_matsu_spin_up')
      do i = 1, size(frequ)
         write (87122, *) frequ(i), real(gw(1, 1, i)), aimag(gw(1, 1, i))
         write (87120, *) frequ(i), real(sigw(1, 1, i)), aimag(gw(1, 1, i))
      enddo
      do i = 1, size(bath%hybrid%freq%vec)
         write (87123, *) aimag(bath%hybrid%freq%vec(i)), real(g_out(1, 1, i)), &
            aimag(g_out(1, 1, i))
         write (87121, *) aimag(bath%hybrid%freq%vec(i)), real(self_out(1, 1, &
                                                                        i)), aimag(self_out(1, 1, i))
      enddo
      close (87120)
      close (87121)
      close (87122)
      close (87123)

   contains

#include "dmft_ed_solver_fmos.h"
#include "dmft_ed_solver_dump_fort.h"
#include "correlations_keldysh.h"
#include "dmft_ed_solver_hubbardI_run.h"

   end subroutine

   subroutine cpt_extract_G(vec, sig, g_out, imaginary)

      use globalvar_ed_solver, only: cpt_correct_green_out, t_cpt
      use matrix, only: id, invmat

      implicit none

      complex(8) :: vec(:), sig(:, :, :), green(size(T_cpt, 1), size(T_cpt, 2)), &
                    gg(size(vec), size(sig, 1), size(sig, 2))
      complex(8) :: g_out(:, :, :)
      integer    :: i, j, k, l, nn
      logical    :: imaginary

      nn = size(T_cpt, 1)
      do i = 1, size(vec)
         green = vec(i)*Id(nn) - T_cpt
         green(1:size(sig, 1), 1:size(sig, 2)) = green(1:size(sig, 1), &
                                                       1:size(sig, 2)) - sig(:, :, i)
         call invmat(nn, green)
         gg(i, 1:size(sig, 1), 1:size(sig, 2)) = green(1:size(sig, 1), &
                                                       1:size(sig, 2))
      enddo

      do i = 1, size(sig, 1)
         if (imaginary) then
            !call plotarray( aimag(vec), aimag(gg(:, i, i)),
            ! 'CPT_GREEN_MATSU_APPROX_diag' // trim(adjustl(tostring(i))) //
            ! "_" )
            !call plotarray( aimag(vec), aimag(g_out(i, i, :)),
            ! 'CPT_GREEN_MATSU_ORIG_diag' // trim(adjustl(tostring(i))) // "_"
            ! )
         else
            !call plotarray( real(vec), aimag(gg(:, i, i)),
            ! 'CPT_GREEN_REAL_APPROX_diag' // trim(adjustl(tostring(i))) //
            ! "_" )
            !call plotarray( real(vec), aimag(g_out(i, i, :)),
            ! 'CPT_GREEN_REAL_ORIG_diag' // trim(adjustl(tostring(i))) // "_"
            ! )
         endif
      enddo

      if (cpt_correct_green_out) then
         do i = 1, size(vec)
            g_out(:, :, i) = gg(i, 1:size(sig, 1), 1:size(sig, 2))
         enddo
      endif

   end subroutine

   function flip_matrix(mat, flip_input_output)

      implicit none

      complex(8) :: mat(:, :), flip_matrix(size(mat, 1), size(mat, 2))
      integer    :: i, j, a1, a2
      logical    :: flip_input_output

      a1 = size(mat, 1)
      a2 = size(mat, 2)
      do i = 1, a1
         do j = 1, a2
            if (.not. flip_input_output .or. .not. (i > a1/2 .and. j > a2/2) &
                ) then
               flip_matrix(i, j) = mat(i, j)
            else
               flip_matrix(i, j) = -conjg(mat(i, j))
            endif
         enddo
      enddo
   end function

   function flip_matrix_Eimp(mat, flip_input_output)

      implicit none

      complex(8) :: mat(:, :), flip_matrix_Eimp(size(mat, 1), size(mat, 2))
      integer    :: i, j, a1, a2, Nc
      logical    :: flip_input_output

      flip_matrix_Eimp = 0.d0
      a1 = size(mat, 1)
      a2 = size(mat, 2)
      Nc = a1/2
      if (.not. flip_input_output) then
         flip_matrix_Eimp = mat
         return
      endif
      flip_matrix_Eimp(1:Nc, 1:Nc) = mat(1:Nc, 1:Nc)
      do i = 1, Nc
         do j = 1, Nc
            flip_matrix_Eimp(Nc + i, Nc + j) = -mat(Nc + j, Nc + i)
         enddo
      enddo
   end function

   subroutine stand_alone_ed()

      use common_def, only: utils_system_call
      use random, only: init_rantab, rand_init
      use genvar, only: imi, matsubara, no_mpi, pi, ran_tab, rank, size2
      use stringmanip, only: tostring
      use impurity_class, only: hamiltonian, web
      use globalvar_ed_solver, only: jhund, jjmatrix, Jhund_Slater_type, uumatrix
      use mesh, only: build1dmesh
      use matrix, only: diag, write_array
      use mpi_mod, only: mpibarrier
      use solver, only: remove_filef

      implicit none

      real(8)                 :: mmu, bbeta, rdelta_width, wwmin, wmax, &
                                 rdelta_frequ_eta1_, rdelta_frequ_T_, &
                                 rdelta_frequ_w0_
      type(hamiltonian)       :: Himp
      type(web)               :: impurity_
      integer                 :: i, nn, ii, jj, kk, nmatsu_frequ, nmatsu_long, &
                                 FLAGIMP, nstep, Nww, bath_param_ed, &
                                 FLAG_AVERAGE_G, FLAG_ORDER_TYPE
      logical                 :: retarded, para_state_, supersc_state_, &
                                 compute_all
      real(8)                 :: tot_rep, frequ_min, frequ_max, dU1, ttemp, &
                                 JJhund
      real(8), allocatable    :: UUmatrix_loc(:, :), JJmatrix_loc(:, :), rdens(:)
      integer                 :: nreal_frequ
      integer                 :: iter_, cluster_problem_size, j, k
      real(8), allocatable    :: frequr(:), frequi(:)
      complex(8), allocatable :: tempmat(:, :)
      complex(8), allocatable :: spm_(:, :, :), corhop_(:, :, :), &
                                 hybrid_in(:, :, :), Eimp(:, :), &
                                 hybrid_in_long(:, :, :)
      complex(8), allocatable :: self_out(:, :, :), sigw(:, :, :), gw(:, :, :), &
                                 g_out(:, :, :)
      real(8), allocatable    :: temp(:), Eimp_(:, :), bathparams(:), &
                                 bathparams_output(:)
      integer(4)              :: my_seed, my_iter_dmft, l
      complex(8)              :: ccc
      logical                 :: nohole, fit_green, checkit, cluster_full, &
                                 check, use_input_delta, no_cdw
      real(8)                 :: ww, ddd, aa
      integer(4)              :: iw
      real(8)                 :: reg, a(3), get, sol(2)
      logical, parameter      :: init_every_time = .true.
      logical                 :: first_time
      real(8)                 :: dummy, dummy_tot

      open (unit=10001, file='PARAMS')
      read (10001, *) cluster_problem_size
      allocate (temp(1 + 2*2*cluster_problem_size), Eimp_(cluster_problem_size, &
                                                          cluster_problem_size), Eimp(2*cluster_problem_size, &
                                                                             2*cluster_problem_size), rdens(2*cluster_problem_size))
      allocate (UUmatrix_loc(cluster_problem_size, cluster_problem_size))
      allocate (JJmatrix_loc(cluster_problem_size, cluster_problem_size))
      JJmatrix_loc = 0.d0
      UUmatrix_loc = 0.d0
      Eimp_ = 0.d0
      Eimp = 0.d0

      do ii = 1, cluster_problem_size
         read (10001, *) (Eimp_(ii, jj), jj=1, cluster_problem_size)
      enddo
      Eimp = 0.d0
      Eimp(1:cluster_problem_size, 1:cluster_problem_size) = Eimp_
      do ii = 1, cluster_problem_size
         read (10001, *) (Eimp_(ii, jj), jj=1, cluster_problem_size)
      enddo
      Eimp(cluster_problem_size + 1:2*cluster_problem_size, &
           cluster_problem_size + 1:2*cluster_problem_size) = Eimp_
      do ii = 1, cluster_problem_size
         read (10001, *) (UUmatrix_loc(ii, jj), jj=1, cluster_problem_size)
      enddo
      read (10001, *) nmatsu_frequ
      read (10001, *) nmatsu_long
      read (10001, *) nreal_frequ
      allocate (spm_(2*cluster_problem_size, 2*cluster_problem_size, &
                     nmatsu_long), corhop_(2*cluster_problem_size, &
                                           2*cluster_problem_size, nmatsu_long))
      allocate (g_out(2*cluster_problem_size, 2*cluster_problem_size, &
                      nmatsu_long))
      allocate (tempmat(cluster_problem_size, cluster_problem_size))
      allocate (hybrid_in(2*cluster_problem_size, 2*cluster_problem_size, &
                          nmatsu_frequ), hybrid_in_long(2*cluster_problem_size, &
                                                        2*cluster_problem_size, nmatsu_long))
      allocate (self_out(2*cluster_problem_size, 2*cluster_problem_size, &
                         nmatsu_long))
      allocate (sigw(2*cluster_problem_size, 2*cluster_problem_size, &
                     nreal_frequ))
      allocate (gw(2*cluster_problem_size, 2*cluster_problem_size, nreal_frequ))
      allocate (frequr(nreal_frequ), frequi(nmatsu_frequ))
      read (10001, *) compute_all
      read (10001, *) FLAG_ORDER_TYPE
      read (10001, *) retarded
      read (10001, *) frequ_min
      read (10001, *) frequ_max
      call build1Dmesh(frequr, nreal_frequ, frequ_min, frequ_max)
      read (10001, *) mmu
      read (10001, *) bbeta
      read (10001, *) dU1
      read (10001, *) rdelta_width
      read (10001, *) rdelta_frequ_eta1_
      read (10001, *) rdelta_frequ_T_
      read (10001, *) rdelta_frequ_w0_
      read (10001, *) my_iter_dmft
      read (10001, *) JJhund
      read (10001, *) use_input_delta
      read (10001, *) no_cdw
      read (10001, *) fit_green
      read (10001, *) nohole
      cluster_full = .true.
      !optional arguments....
      JJmatrix_loc = Jhund
      do ii = 1, cluster_problem_size
         read (10001, *, end=61) (JJmatrix_loc(ii, jj), jj=1, &
                                  cluster_problem_size)
      enddo
      read (10001, *, end=61) cluster_full

61    continue
      close (10001)

      first_time = my_iter_dmft < 0
      my_iter_dmft = abs(my_iter_dmft)

      if (rank == 0 .or. no_mpi) then
         inquire (file='ed.sector_bound_file', exist=checkit)
         if (.not. checkit) then
            open (unit=1451, file='ed.sector_bound_file')
            write (1451, *) - 1000, 1000
            write (1451, *) - 1000, 1000
            write (1451, *) - 1000, 1000
            close (1451)
         endif
      endif

      do ii = 1, cluster_problem_size
         do jj = 1, cluster_problem_size
            if (abs(UUmatrix_loc(ii, jj)) < 1.d-6 .and. abs(JJmatrix_loc(ii, &
                                                                         jj)) > 1.d-6) then
               write (*, *) 'ERROR zero repulsion (so mask repulsion will be &
                    &false)'
               write (*, *) 'and finite hunds coupling, i, j : ', i, j
               write (*, *) 'U, J : ', UUmatrix_loc(ii, jj), JJmatrix_loc(ii, jj)
               write (*, *) 'they share the same mask within the code, this &
                    &should match, critical'
               stop
            endif
         enddo
      enddo

      do ii = 1, cluster_problem_size
         do jj = 1, cluster_problem_size
            if (ii /= jj .and. abs(UUmatrix_loc(ii, jj)) > 1.d-6) then
               if (.not. Jhund_Slater_type) then
                  UUmatrix_loc(ii, jj) = UUmatrix_loc(ii, jj) - &
                                         2.5d0*JJmatrix_loc(ii, jj)
               else
                  UUmatrix_loc(ii, jj) = UUmatrix_loc(ii, jj) - &
                                         2.0d0*JJmatrix_loc(ii, jj)
               endif
            endif
         enddo
      enddo

      if (rank == 0 .or. no_mpi) then
         write (*, *) ' ==================================== '
         write (*, *) '  my_iter_dmft :  ', my_iter_dmft
         write (*, *) ' ==================================== '
      endif

      iter_ = 1

      if (first_time) then
         write (*, *) ' ====================================&
              &======================= '
         write (*, *) 'THIS INTERFACE WAS NOT YET CALLED - INITIALIZATION '
         write (*, *) ' RANK = ', rank
         write (*, *) ' ====================================&
              &======================= '
      endif
      if (.not. allocated(Himp%dU)) allocate (Himp%dU(1), Himp%eps(1), &
                                              Himp%teta(1, 1))
      if (.not. allocated(impurity_%nneigh)) allocate (impurity_%nneigh(1), &
                                                       impurity_%site(cluster_problem_size))
      if (.not. allocated(Himp%Vrep)) allocate (Himp%Vrep(1, 1))

      impurity_%N = cluster_problem_size
      impurity_%nneigh = 0
      if (.not. allocated(impurity_%cadran)) &
         allocate (impurity_%cadran(impurity_%N, 1))
      impurity_%cadran = 0
      impurity_%site = 1
      impurity_%teta = 0
      impurity_%nneigh = 1
      Himp%dU = dU1
      Himp%eps = 0.
      Himp%teta = 0.
      Himp%Vrep = 0.

      my_seed = 4781
#ifndef NO_SYS_CALL
      call utils_system_call("rm -r ED_out")
      call utils_system_call("rm -r AGR")
      call utils_system_call("mkdir AGR")
      call utils_system_call("mkdir ED_out")
      call utils_system_call("rm ./seed")
#else
      write (*, *) 'WARNING not erasing directories AGR and ED_out'
      write (*, *) '        because NO_SYS_CALL was set'
#endif

      call init_rantab(iseed_=my_seed)
      call rand_init(my_seed)

      hybrid_in = 0.

      INQUIRE (file='delta_input_full_1', EXIST=check)
      if (check) then
         open (unit=10001, file='delta_input_full_1', form='unformatted')
         j = 0
         do
            k = cluster_problem_size
            read (10001, end=66) ttemp, tempmat
            j = j + 1
            frequi(j) = ttemp
            hybrid_in(1:k, 1:k, j) = tempmat
         enddo
66       continue
         if (j /= nmatsu_frequ) then
            write (*, *) 'not the right number of matsubara frequencies in &
                 &delta_input_full'
            stop
         endif
         close (10001)
      endif

      INQUIRE (file='delta_input_full_2', EXIST=check)
      if (check) then
         open (unit=10001, file='delta_input_full_2', form='unformatted')
         j = 0
         do
            k = cluster_problem_size
            read (10001, end=67) ttemp, tempmat
            j = j + 1
            hybrid_in(k + 1:2*k, k + 1:2*k, j) = tempmat
         enddo
67       continue
         if (j /= nmatsu_frequ) then
            write (*, *) 'not the right number of matsubara frequencies in &
                 &delta_input_full'
            stop
         endif
         close (10001)
      endif

      if (use_input_delta) then
         write (*, *) 'WARNING WARNING WARNING using input delta to compute &
              &Self, instead of fitted delta'
         write (*, *) 'this will break causality'
      endif

      if (.not. cluster_full) then
         write (*, *) 'WARNING WARNING WARNING removing off diag elements in ED &
              &solver'
         do j = 1, nmatsu_frequ
            do i = 1, 2*cluster_problem_size
               do k = 1, 2*cluster_problem_size
                  if (i /= k) hybrid_in(i, k, j) = 0.d0
               enddo
            enddo
         enddo
      else
         if (rank == 0) then
            do j = 1, nmatsu_frequ
               if (maxval(abs(hybrid_in(:, :, j) - transpose(hybrid_in(:, :, &
                                                                       j)))) > 1.d-3) then
                  write (*, *) 'WARNING hybridization not symmetric !!!'
                  write (*, *) 'frequ : ', j
                  write (*, *)
                  do k = 1, cluster_problem_size
                     write (*, '(200f8.3)') (abs(hybrid_in(k, l, j) - &
                                                 hybrid_in(l, k, j)), l=1, cluster_problem_size)
                  enddo
                  write (*, *) '---------real--------------'
                  do k = 1, cluster_problem_size
                     write (*, '(200f8.3)') (real(hybrid_in(k, l, j)), l=1, &
                                             cluster_problem_size)
                  enddo
                  write (*, *) '---------imag--------------'
                  do k = 1, cluster_problem_size
                     write (*, '(200f8.3)') (aimag(hybrid_in(k, l, j)), l= &
                                             1, cluster_problem_size)
                  enddo
                  !stop
               endif
            enddo
         endif
      endif

      hybrid_in_long = 0.d0
      do k = 1, nmatsu_frequ
         hybrid_in_long(:, :, k) = hybrid_in(:, :, k)
      enddo

      do i = 1, 2*cluster_problem_size
         iw = nmatsu_frequ - 10
         ww = dble(2*iw - 1)*pi/bbeta
         ddd = -Aimag(hybrid_in(i, i, iw))*ww
         aa = Real(hybrid_in(i, i, iw))*(ww**2)
         do k = nmatsu_frequ + 1, nmatsu_long
            ww = dble(2*k - 1)*pi/bbeta
            hybrid_in_long(i, i, k) = ddd/(imi*ww) + aa/(ww**2)
         enddo
         !call plotarray( (/( dble(2*iw-1)*pi/bbeta, iw = 1, nmatsu_long )/),
         ! real(hybrid_in_long(i, i, :)), 'hybridlong_re_' //
         ! trim(adjustl(tostring(i))))
         !call plotarray( (/( dble(2*iw-1)*pi/bbeta, iw = 1, nmatsu_long )/),
         ! aimag(hybrid_in_long(i, i, :)), 'hybridlong_im_' //
         ! trim(adjustl(tostring(i))))
      enddo

      do k = 1, 2*cluster_problem_size
         !call plotarray(frequi, aimag(hybrid_in(k, k, :)), 'diag' //
         ! trim(adjustl(tostring(k))))
      enddo

      if (size(hybrid_in, 2) >= 4) then
         !call plotarray(frequi, aimag(hybrid_in(1, 3, :)), 'offdiag13')
         !call plotarray(frequi, aimag(hybrid_in(1, 4, :)), 'offdiag14')
      endif

      if (rank == 0) write (*, *) 'CALLING ED SOLVER WITH [x] real frequencies &
           &and [y] nmatsu frequencies : ', nreal_frequ, nmatsu_long

      inquire (file='ed.average', exist=checkit)
      if (checkit) then
         open (unit=6613, file='ed.average')
         read (6613, *) FLAG_AVERAGE_G
         close (6613)
      else
         FLAG_AVERAGE_G = 0
      endif

      para_state_ = (FLAG_ORDER_TYPE == 1)
      nstep = 1
      FLAGIMP = 1
      if (rank == 0) write (*, *) 'initialization ed solver'

      call initialize_ed_solver(mmu, bbeta, impurity_, Himp, rdelta_width, &
                                nmatsu_long, FLAGIMP, nstep, nreal_frequ, frequ_min, frequ_max, &
                                para_state_, FLAG_ORDER_TYPE == 3, bath_param_ed, &
                                rdelta_frequ_eta1_=rdelta_frequ_eta1_, rdelta_frequ_T_= &
                                rdelta_frequ_T_, rdelta_frequ_w0_=rdelta_frequ_w0_, init= &
                                init_every_time .or. first_time, average_G_=FLAG_AVERAGE_G, &
                                UUmatrix_in_size=cluster_problem_size, UUmatrix_in= &
                                UUmatrix_loc, fit_green_func_and_not_hybrid_=fit_green, &
                                donot_compute_holepart_=nohole)

      if (size2 > 1 .and. .not. no_mpi) then
         write (*, *) 'ED SOLVER INITIALIZED, RANK = ', rank
         call mpibarrier
      endif

      INQUIRE (file='ed.fit.param', EXIST=check)

      if (allocated(bathparams)) deallocate (bathparams, bathparams_output)
      allocate (bathparams(bath_param_ed), bathparams_output(bath_param_ed))

      if (check) then
         open (unit=1212, file='ed.fit.param')
         read (1212, *, end=78) (bathparams(j), j=1, bath_param_ed)
         write (*, *) 'READING PREVIOUS FIT : ', bathparams
         close (1212)
         if (.false.) then
78          continue
            write (*, *) 'the size of the ED solver bath has changed, erasing &
                 &fitting parameters file'
            close (1212)
#ifndef NO_SYS_CALL
            if (rank == 0 .or. no_mpi) then
               call utils_system_call("rm ed.fit.param")
            endif
#else
            if (rank == 0 .or. no_mpi) then
               call remove_filef("ed.fit.param")
            endif
#endif
            !call system("rm ed.sector_bound_file")
            bathparams = (/(ran_tab(j), j=1, bath_param_ed)/)
         endif
      else
         bathparams = (/((-1.+2*ran_tab(j))/30., j=1, bath_param_ed)/)
         write (*, *) 'STARTING WITH RANDOM PARAM FOR THE FIT : ', bathparams
      endif

      if (bath_param_ed == 0) then
         write (*, *) 'error initialization problem, 0 bath param : ', &
            bath_param_ed
      endif
      if (rank == 0) then
         write (*, *) '############################################'
         write (*, *) 'calling ED solver, with [x] Bath parameters : ', &
            bath_param_ed
         write (*, *) '############################################'
      endif
      self_out = 0.
      g_out = 0.
      gw = 0.
      sigw = 0.

      if (size2 > 1 .and. .not. no_mpi) then
         write (*, *) 'CALLING ED SOLVER, RANK = ', rank
         call mpibarrier
      endif

      write (*, *) 'allocated JJmatrix? :', allocated(JJmatrix), &
         allocated(UUmatrix)

      call solver_ED_interface(bath_param_ed, iter_, mmu, g_out, self_out, &
                               hybrid_in_long, Eimp, sigw, rdens, impurity_, Himp, gw, frequr, &
                               flip_input_output=FLAG_ORDER_TYPE /= 3, para_state_= &
                               FLAG_ORDER_TYPE == 1, retarded=retarded, restarted_=.false., &
                               compute_all=compute_all, tot_rep=tot_rep, spm_=spm_, corhop_ &
                               =corhop_, imp_causality=.true., param_output_= &
                               bathparams_output, param_input_=bathparams, &
                               use_specific_set_parameters_=.true., Jhund_=JJhund, &
                               use_input_delta_instead_of_fit=use_input_delta, Jhund_matrix= &
                               JJmatrix_loc)

      if (size2 > 1 .and. .not. no_mpi) then
         write (*, *) 'ED SOLVER TERMINATED, RANK = ', rank
         call mpibarrier
      endif

#ifndef NO_SYS_CALL
      if (rank == 0 .or. no_mpi) then
         call utils_system_call("rm ed.fit.param")
      endif
#else
      if (rank == 0 .or. no_mpi) then
         call remove_filef("ed.fit.param")
      endif
#endif
      if (rank == 0 .or. no_mpi) then
         open (unit=1212, file='ed.fit.param')
         write (1212, *) (bathparams_output(j), j=1, bath_param_ed)
         write (*, *) 'WRITING NEW FIT : ', bathparams_output
         close (1212)
      endif

      if (rank == 0) then
         write (*, *) ' impurity densities  : ', rdens
         write (*, *) ' total density       : ', sum(rdens)
         write (*, *) ' total double occ    : ', tot_rep
      endif

      do k = 1, cluster_problem_size
         !call plotarray((/( dble(j), j = 1, nmatsu_long)/), real(self_out(k,
         ! k, :)), 'self_out_diag_real_' // trim(adjustl(toString(k))))
         !call plotarray((/( dble(j), j = 1, nmatsu_long)/), aimag(self_out(k,
         ! k, :)), 'self_out_diag_imag_' // trim(adjustl(toString(k))))
      enddo

      if (rank == 0) then
         call write_array(real(self_out(:, :, nmatsu_long)), 'self at matsu &
              &long RE', short=.true., unit=6)
         call write_array(real(g_out(:, :, nmatsu_long)), 'g at matsu long &
              &RE', short=.true., unit=6)
         call write_array(abs(real(self_out(:, :, nmatsu_frequ - 1) - self_out(:, &
              :, nmatsu_long))), ' correction of Eimp nmatsu_long-nmatsu_frequ &
              &RE', short=.true., unit=6)
         call write_array(aimag(self_out(:, :, nmatsu_long)), 'self at matsu &
              &long IM', short=.true., unit=6)
         call write_array(aimag(g_out(:, :, nmatsu_long)), 'g at matsu long &
              &IM', short=.true., unit=6)
         call write_array(abs(aimag(self_out(:, :, nmatsu_frequ - 1) - self_out(:, &
              :, nmatsu_long))), ' correction of Eimp nmatsu_long-nmatsu_frequ &
              &IM', short=.true., unit=6)
      endif

      self_out(:, :, nmatsu_frequ - 1) = self_out(:, :, nmatsu_long)
      self_out(:, :, nmatsu_frequ - 5) = self_out(:, :, nmatsu_long)

      if (no_cdw) then
         k = cluster_problem_size
         nn = 2*k
         do j = 1, nmatsu_frequ
            do ii = 1, nn
               do jj = 1, nn
                  if (jj > nn/4 .and. ii > nn/4 .and. ii <= nn/2 .and. jj <= &
                      nn/2) then
                     self_out(ii - nn/4, jj - nn/4, j) = self_out(ii, jj, j)
                     sigw(ii - nn/4, jj - nn/4, j) = sigw(ii, jj, j)
                  endif
                  if (jj > nn/2 + nn/4 .and. ii > nn/2 + nn/4) then
                     self_out(ii - nn/4, jj - nn/4, j) = self_out(ii, jj, j)
                     sigw(ii - nn/4, jj - nn/4, j) = sigw(ii, jj, j)
                  endif
               enddo
            enddo
         enddo
      endif

      if (rank == 0 .or. no_mpi) then
         open (unit=10001, file='_sigma_output_full_1', form='unformatted')
         k = cluster_problem_size
         do j = 1, nmatsu_frequ
            write (10001) self_out(1:k, 1:k, j)
         enddo
         close (10001)
         open (unit=10001, file='_sigma_output_full_2', form='unformatted')
         k = cluster_problem_size
         do j = 1, nmatsu_frequ
            write (10001) self_out(k + 1:2*k, k + 1:2*k, j)
         enddo
         close (10001)

         open (unit=10001, file='_sigma_output_full_real_1', form= &
               'unformatted')
         k = cluster_problem_size
         do j = 1, nreal_frequ
            write (10001) sigw(1:k, 1:k, j)
         enddo
         close (10001)
         open (unit=10001, file='_sigma_output_full_real_2', form= &
               'unformatted')
         k = cluster_problem_size
         do j = 1, nreal_frequ
            write (10001) sigw(k + 1:2*k, k + 1:2*k, j)
         enddo
         close (10001)

         ! ============= with more frequencies =========== !
         open (unit=10001, file='_sigma_output_full_long_1', form= &
               'unformatted')
         k = cluster_problem_size
         do j = 1, nmatsu_long
            write (10001) self_out(1:k, 1:k, j)
         enddo
         close (10001)
         open (unit=10001, file='_sigma_output_full_long_2', form= &
               'unformatted')
         k = cluster_problem_size
         do j = 1, nmatsu_long
            write (10001) self_out(k + 1:2*k, k + 1:2*k, j)
         enddo
         close (10001)
         ! ============= with more frequencies =========== !

         open (unit=10001, file='sigma_output')
         do j = 1, nmatsu_frequ
            temp(1) = frequi(j)
            do k = 1, 2*cluster_problem_size
               temp(1 + 2*k - 1) = real(self_out(k, k, j))
               temp(1 + 2*k) = aimag(self_out(k, k, j))
            enddo
            write (10001, '(200f15.8)') (temp(k), k=1, &
                                         2*2*cluster_problem_size + 1)
         enddo
         close (10001)

         open (unit=10001, file='sigma_output_real')
         do j = 1, nreal_frequ
            temp(1) = frequr(j)
            do k = 1, 2*cluster_problem_size
               temp(1 + 2*k - 1) = real(sigw(k, k, j))
               temp(1 + 2*k) = aimag(sigw(k, k, j))
            enddo
            write (10001, '(200f15.8)') (temp(k), k=1, &
                                         2*2*cluster_problem_size + 1)
         enddo
         close (10001)

         open (unit=10001, file='green_output_real')
         do j = 1, nreal_frequ
            temp(1) = frequr(j)
            do k = 1, 2*cluster_problem_size
               temp(1 + 2*k - 1) = real(gw(k, k, j))
               temp(1 + 2*k) = aimag(gw(k, k, j))
            enddo
            write (10001, '(200f15.8)') (temp(k), k=1, &
                                         2*2*cluster_problem_size + 1)
         enddo
         close (10001)

         open (unit=10001, file='green_output_matsu')
         do j = 1, nmatsu_frequ
            temp(1) = frequi(j)
            do k = 1, 2*cluster_problem_size
               temp(1 + 2*k - 1) = real(g_out(k, k, j))
               temp(1 + 2*k) = aimag(g_out(k, k, j))
            enddo
            write (10001, '(200f15.8)') (temp(k), k=1, &
                                         2*2*cluster_problem_size + 1)
         enddo
         close (10001)

         do k = 1, 2
            open (unit=10001, file='green_output_matsu_full'// &
                  TRIM(ADJUSTL(toString(k))), form='unformatted')
            do j = 1, nmatsu_frequ
               i = cluster_problem_size
               if (k == 1) then
                  write (10001) g_out(1:i, 1:i, j)
               else
                  write (10001) g_out(i + 1:2*i, i + 1:2*i, j)
               endif
            enddo
            close (10001)
         enddo

         ! ============= with more frequencies =========== !
         do k = 1, 2
            open (unit=10001, file='green_output_matsu_full_long_'// &
                  TRIM(ADJUSTL(toString(k))), form='unformatted')
            do j = 1, nmatsu_long
               i = cluster_problem_size
               if (k == 1) then
                  write (10001) g_out(1:i, 1:i, j)
               else
                  write (10001) g_out(i + 1:2*i, i + 1:2*i, j)
               endif
            enddo
            close (10001)
         enddo
         ! ============= with more frequencies =========== !
      endif

      write (*, *) 'ED CLEANING UP, RANK = ', rank

      if (init_every_time) then
         deallocate (UUmatrix_loc, JJmatrix_loc)
         deallocate (UUmatrix, JJmatrix, rdens, frequr, frequi, tempmat, &
                     spm_, corhop_, hybrid_in, Eimp, hybrid_in_long, self_out, sigw, &
                     gw, g_out, temp, Eimp_, bathparams, bathparams_output)
         deallocate (Himp%dU, Himp%eps, Himp%teta)
         deallocate (impurity_%nneigh, impurity_%site)
         deallocate (Himp%Vrep)
         deallocate (impurity_%cadran)
         call finalize_ed_solver(finalize_iter=.true., finalize=.true.)
      endif

      write (*, *) 'ED DONE RETURN, RANK = ', rank

      return
   end subroutine

   subroutine initialize_ed_solver(mmu, bbeta, impurity_, Himp, rdelta_width, &
                                   nmatsu_frequ, FLAGIMP, nstep, Nww, wwmin, wmax, para_state_, &
                                   supersc_state_, bath_param_ed, rdelta_frequ_eta1_, rdelta_frequ_T_, &
                                   rdelta_frequ_w0_, init, average_G_, UUmatrix_in_size, UUmatrix_in, &
                                   fit_green_func_and_not_hybrid_, donot_compute_holepart_)

      use globalvar_ed_solver, only: average_g, beta_ed, demax, demax0, &
         donot_compute_holepart, fast_fit, fit_green_func_and_not_hybrid, &
         flag_slater_int, fmos, force_no_pairing, force_nupdn_basis, force_para_state, gen_cpt, &
         jjmatrix, min_all_bath_param, ncpt_approx, ncpt_approx_tot, &
         ncpt_chain_coup, ncpt_para, ncpt_tot, para_state, &
         rdelta_frequ_eta1, rdelta_frequ_t, rdelta_frequ_w0, &
         slater_coulomb_c, slater_coulomb_r, supersc_state, ucc, uccc, uccr, &
         use_input_eimp_matrices, use_precomputed_slater_matrix, uumatrix
      use impurity_class, only: define_impurity, hamiltonian, web
      use genvar, only: fermionic, log_unit, rank, size2
      use aim_class, only: new_aim, update_aim_pointer
      use bath_class_hybrid, only: bath2hybrid
      use common_def, only: utils_assert
      use correl_class, only: new_correl
      use stringmanip, only: tostring
      use strings, only: remove
      use bath_class, only: read_bath
      use matrix, only: write_array
      use correlations, only: init_correlations, SNAMBU, SNAMBUret
      use solver, only: init_solver

      implicit none

      real(8)           :: mmu, bbeta, rdelta_width, wmin, wmax, wwmin
      type(hamiltonian) :: Himp
      type(web)         :: impurity_
      integer           :: nmatsu_frequ, FLAGIMP, nstep, Nww, bath_param_ed
      logical           :: check, para_state_, supersc_state_
      real(8)           :: rdelta_frequ_eta1_, rdelta_frequ_T_, rdelta_frequ_w0_
      logical, optional :: init, fit_green_func_and_not_hybrid_, &
                           donot_compute_holepart_
      integer, optional :: average_G_
      integer, optional :: UUmatrix_in_size
      real(8), optional :: UUmatrix_in(impurity_%N, impurity_%N)
      integer           :: chan, g1, g2, g3, g4

      use_precomputed_slater_matrix = .false.

      inquire (file='./Uc.dat.ED', exist=flag_slater_int)
      if (flag_slater_int) Then
         write (log_unit, *) '#############################################'
         write (log_unit, *) '       USING SLATER INTERACTIONS             '
         open (unit=11122, file='./Uc.dat.ED', form='unformatted')
         read (11122) chan
         if (chan /= impurity_%N) then
            write (log_unit, *) 'ERROR, slater file contains [x] orbital sets : &
                 &', chan
            write (log_unit, *) ' and impurity model contains [x] orb : ', &
               impurity_%N
            stop
         endif
         if (chan > 7) then
            write (log_unit, *) 'ERROR, slater file contains too many orbitals &
                 &(max is 6 for memory reason): ', chan
            write (log_unit, *) ' remove this failsafe in case you want to do &
                 &more...'
            stop
         endif
         if (allocated(Slater_Coulomb_c)) deallocate (Slater_Coulomb_r, &
                                                      Slater_Coulomb_c)
         allocate (Slater_Coulomb_c(chan, chan, chan, chan), &
                   Slater_Coulomb_r(chan, chan, chan, chan))
         read (11122) Slater_Coulomb_c
         close (11122)
         Slater_Coulomb_r = Slater_Coulomb_c
#ifndef _complex
         if (maxval(abs(Slater_Coulomb_c - Slater_Coulomb_r)) > 1.d-2) then
            write (*, *) 'Slater Coulomb interaction is complex, and running &
                 &Lanczos in real mode'
            write (*, *) 'please recompile with _complex flag, to run the &
                 &complex Lanczos'
            write (*, *) 'maxval Slater_Coulomb_c : ', &
               maxval(abs(Slater_Coulomb_c))
            write (*, *) 'maxval complex-real : ', &
               maxval(abs(Slater_Coulomb_c - Slater_Coulomb_r))
            write (*, *) 'critical'
            stop
         endif
         Slater_Coulomb_c = conjg(Slater_Coulomb_c) !we always use the
         ! conjugate of Slater_Coulomb_c
#endif
         inquire (file='Uc.dat.ED.matrix', exist= &
                  use_precomputed_slater_matrix)
         if (use_precomputed_slater_matrix) then
            if (allocated(UCC)) deallocate (UCCr, UCCc, UCC)
            open (unit=71123, file='Uc.dat.ED.matrix', form='unformatted')
            read (71123) g1, g2, g3, g4
            allocate (UCC(0:g1 - 1, 0:g2 - 1, g3, g4))
            read (71123) UCC(0:g1 - 1, 0:g2 - 1, 1:g3, 1:g4)
            read (71123) g4, g1, g2, g3
            if (g4 /= chan) then
               write (*, *) 'OUPS, the file Uc.dat.ED.matrix was for [x] &
                    &orbitals : ', g1
               stop
            endif
            allocate (UCCr(0:g1 - 1, 0:g2 - 1, 1:g3), UCCc(0:g1 - 1, 0:g2 - 1, 1:g3))
            read (71123) UCCc(0:g1 - 1, 0:g2 - 1, 1:g3)
            UCCr = real(UCCc, 8)
            close (71123)
         endif
         write (log_unit, *) '#############################################'
      endif

      donot_compute_holepart = .false.
      if (present(donot_compute_holepart_)) then
         donot_compute_holepart = donot_compute_holepart_
      endif

      use_input_Eimp_matrices = .false.
      if (allocated(UUmatrix)) deallocate (UUmatrix)
      if (allocated(JJmatrix)) deallocate (JJmatrix)
      allocate (JJmatrix(impurity_%N, impurity_%N))

      if (present(UUmatrix_in)) then
         allocate (UUmatrix(UUmatrix_in_size, UUmatrix_in_size))
         UUmatrix = UUmatrix_in
         use_input_Eimp_matrices = .true.
         call utils_assert(impurity_%N == UUmatrix_in_size, 'Error in dmft_solver_ed: &
               & size-mismatch between cluster problem (' // trim(adjustl(tostring(impurity_%N))) &
               & // ') and correlation matrix (' // trim(adjustl(tostring(UUmatrix_in_size))) // ')')
         if (flag_slater_int) UUmatrix = 0.d0
         call write_array(UUmatrix, ' repulsion matrix input, please check and &
              &comment out ', short=.true., unit=6)
         if (UUmatrix_in_size /= impurity_%N) then
            write (*, *) 'error in the inputs, UUmatrix_in'
            stop
         endif
      endif

      if (flag_slater_int) Himp%dU = 0.d0

      rdelta_frequ_eta1 = rdelta_frequ_eta1_
      rdelta_frequ_T = rdelta_frequ_T_
      rdelta_frequ_w0 = rdelta_frequ_w0_

      fit_green_func_and_not_hybrid = .false.
      if (present(fit_green_func_and_not_hybrid_)) &
         fit_green_func_and_not_hybrid = fit_green_func_and_not_hybrid_

      average_G = 0
      if (present(average_G_)) average_G = average_G_

      para_state = para_state_
      supersc_state = supersc_state_

      beta_ED = bbeta
      ! =============================================== !
      if (present(init)) then
         if (init) CALL read_DMFT_parameters(nmatsu_frequ, FLAGIMP, nstep, &
                                             impurity_%N)
      endif
#ifdef OPENMP_MPI_SAFE
      if (MAXT > 1 .and. size2 > 1 .and. .not. no_mpi) then
         fast_fit = .false.
      endif
#endif
      ! =============================================== !
      dEmax = 1.d0/beta_ED*dEmax0

      if (supersc_state .and. .not. force_nupdn_basis) min_all_bath_param = &
         -abs(min_all_bath_param)
      if (force_no_pairing) min_all_bath_param = abs(min_all_bath_param)

      write (145 + rank, *) ' .... dEmax .... : ', dEmax

      wmin = -max(abs(wwmin), abs(wmax))
      wmax = max(abs(wwmin), abs(wmax))

      write (log_unit, *) ' ==================================&
           &============== '
      write (log_unit, *) '.... DEFINE IMPURITY ....'
      CALL define_impurity(impurity, mmu, impurity_, Himp)
      write (log_unit, *) '......INITIALIZE CORRELATIONS MATRIX......'
      CALL init_correlations(impurity, bbeta, CORRELfile, rdelta_width, wmin, &
                             wmax, nmatsu_frequ, Nww)
      write (log_unit, *) ' ==================================&
           &============== '

      ! =============================================== !
      if (present(init)) then
         if (init) CALL read_bath(bath, BATHfile, impurity%Nc)
      endif
      ! =============================================== !

      write (log_unit, *) '.... copy SNAMBU correl to hybrid bath ....'
      CALL new_correl(bath%hybrid, SNAMBU)
      bath%hybrid%title = 'hybrid'
      write (log_unit, *) 'copy SNAMBUret to hybrid retarded....'
      CALL new_correl(bath%hybridret, SNAMBUret)
      bath%hybridret%title = 'hybridret'
      write (log_unit, *) '.... bath2hybrid ....'
      CALL bath2hybrid(bath, FERMIONIC)
      write (log_unit, *) '.... new AIM ....'
      call update_AIM_pointer(AIM, impurity, bath)

      ! =============================================== !
      if (present(init)) then
         if (init) CALL new_AIM(AIM, impurity, bath)
      endif
      ! =============================================== !

      if (fmos) then
         bath_param_ed = bath%nparam
         return
      endif

      IF (start_from_old_gs) THEN
         write (log_unit, *) '... init solver from OLD file ...'
         CALL init_solver(EDfile, AIM, OLDGSFILE=OLDGSFILE)
      ELSE
         write (log_unit, *) '... init solver from scratch ...'
         CALL init_solver(EDfile, AIM)
      ENDIF

      CALL write_header(AIM, log_unit)

      if (use_precomputed_slater_matrix .and. allocated(UCC)) then
         UCC(:, :, :, 1) = UCC(:, :, :, 1)*(2**AIM%Nb)
         UCC(:, :, :, 2) = UCC(:, :, :, 2)*(2**AIM%Nb)
      endif

      !For cluster perturbation theory:
      if (ncpt_approx /= 0) then
         if (gen_cpt) then
            ncpt_chain_coup = ncpt_approx
         else
            ncpt_chain_coup = 1
         endif
         ncpt_approx_tot = AIM%Nb*ncpt_approx
         ncpt_tot = AIM%Nb*(ncpt_approx*(ncpt_approx + 1)/2) + &
                    AIM%Nb*ncpt_chain_coup
      else
         ncpt_approx_tot = 0
         ncpt_tot = 0
         ncpt_chain_coup = 0
      endif

      if (para_state .or. force_para_state) then
         if (supersc_state .and. .not. force_nupdn_basis .and. &
             .not. force_no_pairing) then
            ncpt_para = 2
         else
            ncpt_para = 1
         endif
      else
         if (supersc_state .and. .not. force_nupdn_basis .and. &
             .not. force_no_pairing) then
            ncpt_para = 3
         else
            ncpt_para = 2
         endif
      endif

      bath_param_ed = bath%nparam + ncpt_para*ncpt_tot
      write (log_unit, *) ' NUMBER OF BATH PARAMETERS : ', bath_param_ed + &
         ncpt_para*ncpt_tot

   end subroutine

   subroutine finalize_ed_solver(finalize_iter, finalize)

      use genvar, only: log_unit
      use impurity_class, only: delete_impurity
      use aim_class, only: delete_aim

      implicit none

      logical, optional :: finalize, finalize_iter

      if (present(finalize_iter)) then
         CALL delete_impurity(impurity)
         CALL write_footer(AIM, GS, log_unit)
      endif
      if (present(finalize)) then
         call delete_AIM(AIM)
      endif
      return
   end subroutine

end module
