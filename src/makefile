# Second-level makefile for TOSCAM master scripts
#
# Do not edit; system-specific changes should be made to the appropriate
#  conf.arch files
#

ifndef ROOTDIR
   ROOTDIR=..
endif
ifdef ARCH
   include $(ROOTDIR)/config/conf.$(ARCH)
else
   include $(ROOTDIR)/config/conf.default
endif

prep = -fpp -cxxlib

.SUFFIXES:
.SUFFIXES: .mod .f90 .o .a

#
# Libraries
#

LIBUTILS   = -L$(ROOTDIR)/lib/utils/obj -lutils
LIBSPLINES = -L$(ROOTDIR)/lib/splines/obj -lsplines

LIBS       = $(MKL) $(LIBUTILS) $(LIBSPLINES)

#
# Include
#

VPATH = $(ROOTDIR)/src:$(ROOTDIR)/ed_solver/obj:$(ROOTDIR)/lib/utils/obj

INCUTILS    = -I$(ROOTDIR)/lib/utils/obj
INCSOLVER   = -I$(ROOTDIR)/ed_solver/obj
INCSCRIPTS  = -I$(ROOTDIR)/src/obj
INC         = $(INCSOLVER) $(INCUTILS) $(INCSCRIPTS)

#
# Module redirection
#

ifeq ($(F90),$(filter $(F90),ifort mpifort))
   REDIRECT = -module $(ROOTDIR)/src/obj
endif
ifeq ($(F90),$(filter $(F90),gfortran, mpif90))
   REDIRECT = -J$(ROOTDIR)/src/obj
endif

#
# Scripts
#

SCRIPTS = ed_solver.out \
          onetep_dmft_one_iter.out \
          dmft_gen_correl_output_file.out \
          onetep_dmft_embedding_gen_sigma.out \
          dmft_collect_sigma_embed.out \
          onetep_dmft_stitch.out \
          onetep_dmft_build_heterostructure.out \
          onetep_dmft_sc.out \
          onetep_dmft_calculate_average_e.out \
          onetep_dmft_build_optical_tensor.out \
          onetep_dmft_read_n.out \
          onetep_dmft_local_rotation.out \
          onetep_dmft_read_spin.out \
          onetep_dmft_trace_energy.out \
          check_convergence_full.out \
          check_convergence_single_renorm.out \
          check_convergence_renorm.out \
          onetep_dmft_histogram.out \
          check_convergence.out \
          onetep_dmft_read_overlap_matrices.out \
          ed_solver_serial.out \
          onetep_dmft_movie.out \
          onetep_dmft_local_rotation_find.out \
          onetep_dmft_local_rotation_one_atom.out \
          onetep_dmft_proj_green_func_matrix.out \
          onetep_dmft_local_dimers.out \
          dmft_all_iterations_onetep.out \
          onetep_dmft_compute_optics.out \
          onetep_dmft_summarize_sim.out \
          onetep_dmft_density_matsu.out \
          onetep_dmft_local_uniform.out \
          onetep_dmft_duplicate_unitcell.out \
          onetep_dmft.out \
          onetep_split.out \
          onetep_split_serial.out \
          dmft_collect.out \
          dmft_project.out \
          dmft_project_back.out \
          dmft_group_dimers.out \
          dmft_gen_correl.out \
          dmft_project_sigma.out \
          dmft_collect_dos.out

#
# Specific compilation rules
#

default: build

build: objdir $(SCRIPTS)

clean:
	$(RM) obj/* *.o *.mod

objdir:
	@( if [ ! -d $(ROOTDIR)/src/obj ] ; \
           then mkdir $(ROOTDIR)/src/obj ; \
   fi )

#
# General rules
#

%.out: %.f90
	$(F90) $(prep) -o $(ROOTDIR)/bin/$@ $(FFLAGS) $^ $(LIBS) $(INC) $(REDIRECT)

#
# Dependencies
#

check_convergence.out: check_convergence.f90 aim2_class.o aim_class.o apply_c.o apply_ns.o \
    apply_p.o bath_class.o bath_class_hybrid.o bath_class_vec.o block_lanczos.o common_def.o \
    conj_grad.o correl_class.o correlations.o density_matrix.o derivative_noise.o \
    dmft_solver_ed.o ed_arpack.o eigen_class.o eigen_sector_class.o fermion_hilbert_class.o \
    fermion_ket_class.o fermion_sector2_class.o fourier_transform_mod.o frequency_class.o genvar.o \
    geometry.o geomlib.o globalvar_ed_solver.o green_class.o green_class_compute_dynamic.o \
    green_class_compute_symmetric.o green_class_computeaa.o green_class_computeab.o h_class.o \
    haim2_class.o haimsz_class.o haimupdo_class.o impurity_class.o init_and_close_my_sim.o \
    lanczos_fast.o linalg.o lockmod.o mask_class.o masked_matrix_class.o masked_matrix_class_mod.o \
    matrix.o mesh.o minimization_wrapping.o mpirout.o namelistmod.o openmpmod.o overlap_module.o \
    quantum_algebra.o random.o rcmatrix_class.o rcvector_class.o readable_vec_class.o \
    sector_class.o smooth_data.o solver.o sorting.o specialfunction.o splines.o splines2.o \
    string4.o string5.o stringmanip.o strings.o tools_algebra.o tools_fit.o tridiag_class.o \
    vertex.o

check_convergence_full.out: check_convergence_full.f90 aim2_class.o aim_class.o apply_c.o \
    apply_ns.o apply_p.o bath_class.o bath_class_hybrid.o bath_class_vec.o block_lanczos.o \
    common_def.o conj_grad.o correl_class.o correlations.o density_matrix.o \
    derivative_noise.o dmft_solver_ed.o ed_arpack.o eigen_class.o eigen_sector_class.o \
    fermion_hilbert_class.o fermion_ket_class.o fermion_sector2_class.o fourier_transform_mod.o \
    frequency_class.o genvar.o geometry.o geomlib.o globalvar_ed_solver.o green_class.o \
    green_class_compute_dynamic.o green_class_compute_symmetric.o green_class_computeaa.o \
    green_class_computeab.o h_class.o haim2_class.o haimsz_class.o haimupdo_class.o \
    impurity_class.o init_and_close_my_sim.o lanczos_fast.o linalg.o lockmod.o mask_class.o \
    masked_matrix_class.o masked_matrix_class_mod.o matrix.o mesh.o minimization_wrapping.o \
    mpirout.o namelistmod.o openmpmod.o overlap_module.o quantum_algebra.o random.o \
    rcmatrix_class.o rcvector_class.o readable_vec_class.o sector_class.o smooth_data.o solver.o \
    sorting.o specialfunction.o splines.o splines2.o string4.o string5.o stringmanip.o strings.o \
    tools_algebra.o tools_fit.o tridiag_class.o vertex.o

check_convergence_renorm.out: check_convergence_renorm.f90 aim2_class.o aim_class.o apply_c.o \
    apply_ns.o apply_p.o bath_class.o bath_class_hybrid.o bath_class_vec.o block_lanczos.o \
    common_def.o conj_grad.o correl_class.o correlations.o density_matrix.o \
    derivative_noise.o dmft_solver_ed.o ed_arpack.o eigen_class.o eigen_sector_class.o \
    fermion_hilbert_class.o fermion_ket_class.o fermion_sector2_class.o fourier_transform_mod.o \
    frequency_class.o genvar.o geometry.o geomlib.o globalvar_ed_solver.o green_class.o \
    green_class_compute_dynamic.o green_class_compute_symmetric.o green_class_computeaa.o \
    green_class_computeab.o h_class.o haim2_class.o haimsz_class.o haimupdo_class.o \
    impurity_class.o init_and_close_my_sim.o lanczos_fast.o linalg.o lockmod.o mask_class.o \
    masked_matrix_class.o masked_matrix_class_mod.o matrix.o mesh.o minimization_wrapping.o \
    mpirout.o namelistmod.o openmpmod.o overlap_module.o quantum_algebra.o random.o \
    rcmatrix_class.o rcvector_class.o readable_vec_class.o sector_class.o smooth_data.o solver.o \
    sorting.o specialfunction.o splines.o splines2.o string4.o string5.o stringmanip.o strings.o \
    tools_algebra.o tools_fit.o tridiag_class.o vertex.o

check_convergence_single_renorm.out: check_convergence_single_renorm.f90 aim2_class.o aim_class.o \
    apply_c.o apply_ns.o apply_p.o bath_class.o bath_class_hybrid.o bath_class_vec.o \
    block_lanczos.o common_def.o conj_grad.o correl_class.o correlations.o \
    density_matrix.o derivative_noise.o dmft_solver_ed.o ed_arpack.o eigen_class.o \
    eigen_sector_class.o fermion_hilbert_class.o fermion_ket_class.o fermion_sector2_class.o \
    fourier_transform_mod.o frequency_class.o genvar.o geometry.o geomlib.o globalvar_ed_solver.o \
    green_class.o green_class_compute_dynamic.o green_class_compute_symmetric.o \
    green_class_computeaa.o green_class_computeab.o h_class.o haim2_class.o haimsz_class.o \
    haimupdo_class.o impurity_class.o init_and_close_my_sim.o lanczos_fast.o linalg.o lockmod.o \
    mask_class.o masked_matrix_class.o masked_matrix_class_mod.o matrix.o mesh.o \
    minimization_wrapping.o mpirout.o namelistmod.o openmpmod.o overlap_module.o quantum_algebra.o \
    random.o rcmatrix_class.o rcvector_class.o readable_vec_class.o sector_class.o smooth_data.o \
    solver.o sorting.o specialfunction.o splines.o splines2.o string4.o string5.o stringmanip.o \
    strings.o tools_algebra.o tools_fit.o tridiag_class.o vertex.o

dmft_all_iterations_onetep.out: dmft_all_iterations_onetep.f90 stringmanip.o

dmft_check_cluster_convergence.out: dmft_check_cluster_convergence.f90 aim2_class.o aim_class.o \
    apply_c.o apply_ns.o apply_p.o bath_class.o bath_class_hybrid.o bath_class_vec.o \
    block_lanczos.o common_def.o conj_grad.o correl_class.o correlations.o \
    density_matrix.o derivative_noise.o dmft_solver_ed.o ed_arpack.o eigen_class.o \
    eigen_sector_class.o fermion_hilbert_class.o fermion_ket_class.o fermion_sector2_class.o \
    fourier_transform_mod.o frequency_class.o genvar.o geometry.o geomlib.o globalvar_ed_solver.o \
    green_class.o green_class_compute_dynamic.o green_class_compute_symmetric.o \
    green_class_computeaa.o green_class_computeab.o h_class.o haim2_class.o haimsz_class.o \
    haimupdo_class.o impurity_class.o init_and_close_my_sim.o lanczos_fast.o linalg.o lockmod.o \
    mask_class.o masked_matrix_class.o masked_matrix_class_mod.o matrix.o mesh.o \
    minimization_wrapping.o mpirout.o namelistmod.o openmpmod.o overlap_module.o quantum_algebra.o \
    random.o rcmatrix_class.o rcvector_class.o readable_vec_class.o sector_class.o smooth_data.o \
    solver.o sorting.o specialfunction.o splines.o splines2.o string4.o string5.o stringmanip.o \
    strings.o tools_algebra.o tools_fit.o tridiag_class.o vertex.o

dmft_collect.out: dmft_collect.f90 strings.o

dmft_collect_dos.out: dmft_collect_dos.f90 strings.o

dmft_collect_sigma_embed.out: dmft_collect_sigma_embed.f90 strings.o

dmft_gen_correl.out: dmft_gen_correl.f90 stringmanip.o

dmft_gen_correl_output_file.out: dmft_gen_correl_output_file.f90 stringmanip.o

dmft_group_dimers.out: dmft_group_dimers.f90 stringmanip.o

dmft_project.out: dmft_project.f90 common_def.o genvar.o geomlib.o linalg.o lockmod.o matrix.o \
    openmpmod.o random.o sorting.o stringmanip.o strings.o tools_algebra.o

dmft_project_back.out: dmft_project_back.f90 common_def.o genvar.o geomlib.o linalg.o lockmod.o \
    matrix.o openmpmod.o random.o sorting.o stringmanip.o strings.o tools_algebra.o

dmft_project_orbitals_back.out: dmft_project_orbitals_back.f90 common_def.o genvar.o geomlib.o \
    linalg.o lockmod.o matrix.o openmpmod.o random.o sorting.o stringmanip.o strings.o \
    tools_algebra.o

dmft_project_sigma.out: dmft_project_sigma.f90 common_def.o genvar.o geomlib.o linalg.o lockmod.o \
    matrix.o openmpmod.o random.o sorting.o stringmanip.o strings.o tools_algebra.o

ed_solver.out: ed_solver.f90 aim2_class.o aim_class.o apply_c.o apply_ns.o apply_p.o bath_class.o \
    bath_class_hybrid.o bath_class_vec.o block_lanczos.o common_def.o conj_grad.o \
    correl_class.o correlations.o density_matrix.o derivative_noise.o dmft_solver_ed.o ed_arpack.o \
    eigen_class.o eigen_sector_class.o fermion_hilbert_class.o fermion_ket_class.o \
    fermion_sector2_class.o fourier_transform_mod.o frequency_class.o genvar.o geometry.o \
    geomlib.o globalvar_ed_solver.o green_class.o green_class_compute_dynamic.o \
    green_class_compute_symmetric.o green_class_computeaa.o green_class_computeab.o h_class.o \
    haim2_class.o haimsz_class.o haimupdo_class.o impurity_class.o init_and_close_my_sim.o \
    lanczos_fast.o linalg.o lockmod.o mask_class.o masked_matrix_class.o masked_matrix_class_mod.o \
    matrix.o mesh.o minimization_wrapping.o mpirout.o namelistmod.o openmpmod.o overlap_module.o \
    quantum_algebra.o random.o rcmatrix_class.o rcvector_class.o readable_vec_class.o \
    sector_class.o smooth_data.o solver.o sorting.o specialfunction.o splines.o splines2.o \
    string4.o string5.o stringmanip.o strings.o tools_algebra.o tools_fit.o tridiag_class.o \
    vertex.o

ed_solver_complex.out: ed_solver_complex.f90 aim2_class.o aim_class.o apply_c.o apply_ns.o \
    apply_p.o bath_class.o bath_class_hybrid.o bath_class_vec.o block_lanczos.o common_def.o \
    conj_grad.o correl_class.o correlations.o density_matrix.o derivative_noise.o \
    dmft_solver_ed.o ed_arpack.o eigen_class.o eigen_sector_class.o fermion_hilbert_class.o \
    fermion_ket_class.o fermion_sector2_class.o fourier_transform_mod.o frequency_class.o genvar.o \
    geometry.o geomlib.o globalvar_ed_solver.o green_class.o green_class_compute_dynamic.o \
    green_class_compute_symmetric.o green_class_computeaa.o green_class_computeab.o h_class.o \
    haim2_class.o haimsz_class.o haimupdo_class.o impurity_class.o init_and_close_my_sim.o \
    lanczos_fast.o linalg.o lockmod.o mask_class.o masked_matrix_class.o masked_matrix_class_mod.o \
    matrix.o mesh.o minimization_wrapping.o mpirout.o namelistmod.o openmpmod.o overlap_module.o \
    quantum_algebra.o random.o rcmatrix_class.o rcvector_class.o readable_vec_class.o \
    sector_class.o smooth_data.o solver.o sorting.o specialfunction.o splines.o splines2.o \
    string4.o string5.o stringmanip.o strings.o tools_algebra.o tools_fit.o tridiag_class.o \
    vertex.o

ed_solver_serial.out: ed_solver_serial.f90 aim2_class.o aim_class.o apply_c.o apply_ns.o apply_p.o \
    bath_class.o bath_class_hybrid.o bath_class_vec.o block_lanczos.o common_def.o \
    conj_grad.o correl_class.o correlations.o density_matrix.o derivative_noise.o \
    dmft_solver_ed.o ed_arpack.o eigen_class.o eigen_sector_class.o fermion_hilbert_class.o \
    fermion_ket_class.o fermion_sector2_class.o fourier_transform_mod.o frequency_class.o genvar.o \
    geometry.o geomlib.o globalvar_ed_solver.o green_class.o green_class_compute_dynamic.o \
    green_class_compute_symmetric.o green_class_computeaa.o green_class_computeab.o h_class.o \
    haim2_class.o haimsz_class.o haimupdo_class.o impurity_class.o init_and_close_my_sim.o \
    lanczos_fast.o linalg.o lockmod.o mask_class.o masked_matrix_class.o masked_matrix_class_mod.o \
    matrix.o mesh.o minimization_wrapping.o mpirout.o namelistmod.o openmpmod.o overlap_module.o \
    quantum_algebra.o random.o rcmatrix_class.o rcvector_class.o readable_vec_class.o \
    sector_class.o smooth_data.o solver.o sorting.o specialfunction.o splines.o splines2.o \
    string4.o string5.o stringmanip.o strings.o tools_algebra.o tools_fit.o tridiag_class.o \
    vertex.o

gensig.out: gensig.f90

onetep_dmft.out: onetep_dmft.f90 common_def.o genvar.o geomlib.o init_and_close_my_sim.o linalg.o \
    lockmod.o namelistmod.o openmpmod.o random.o string4.o string5.o stringmanip.o strings.o

onetep_dmft_build_heterostructure.out: onetep_dmft_build_heterostructure.f90 common_def.o genvar.o \
    lockmod.o namelistmod.o openmpmod.o string4.o string5.o stringmanip.o

onetep_dmft_build_optical_tensor.out: onetep_dmft_build_optical_tensor.f90 common_def.o \
    derivative_noise.o fourier_transform_mod.o genvar.o geometry.o geomlib.o \
    init_and_close_my_sim.o linalg.o lockmod.o matrix.o mesh.o mpirout.o namelistmod.o openmpmod.o \
    random.o smooth_data.o sorting.o splines.o splines2.o string4.o string5.o stringmanip.o \
    tools_algebra.o tools_fit.o

onetep_dmft_calculate_average_e.out: onetep_dmft_calculate_average_e.f90

onetep_dmft_compute_optics.out: onetep_dmft_compute_optics.f90 common_def.o derivative_noise.o \
    fourier_transform_mod.o genvar.o geometry.o geomlib.o linalg.o lockmod.o matrix.o mesh.o \
    mpirout.o openmpmod.o random.o smooth_data.o sorting.o splines.o splines2.o stringmanip.o \
    tools_algebra.o tools_fit.o

onetep_dmft_density_matsu.out: onetep_dmft_density_matsu.f90 genvar.o geomlib.o linalg.o lockmod.o \
    openmpmod.o stringmanip.o

onetep_dmft_duplicate_unitcell.out: onetep_dmft_duplicate_unitcell.f90

onetep_dmft_embedding_gen_sigma.out: onetep_dmft_embedding_gen_sigma.f90 common_def.o genvar.o \
    lockmod.o namelistmod.o openmpmod.o string4.o string5.o stringmanip.o

onetep_dmft_fmos.out: onetep_dmft_fmos.f90 routinesfmos.o

onetep_dmft_hf.out: onetep_dmft_hf.f90 routineshf.o

onetep_dmft_histogram.out: onetep_dmft_histogram.f90 common_def.o genvar.o geomlib.o linalg.o \
    lockmod.o openmpmod.o random.o sorting.o stringmanip.o

onetep_dmft_local_dimers.out: onetep_dmft_local_dimers.f90 common_def.o genvar.o geomlib.o \
    linalg.o lockmod.o openmpmod.o random.o sorting.o

onetep_dmft_local_rotation.out: onetep_dmft_local_rotation.f90 genvar.o geometry2.o geomlib.o \
    linalg.o lockmod.o openmpmod.o

onetep_dmft_local_rotation_find.out: onetep_dmft_local_rotation_find.f90 common_def.o \
    conj_grad.o genvar.o geomlib.o init_and_close_my_sim.o linalg.o lockmod.o \
    minimization_wrapping.o mpirout.o namelistmod.o openmpmod.o random.o string4.o string5.o \
    stringmanip.o

onetep_dmft_local_rotation_one_atom.out: onetep_dmft_local_rotation_one_atom.f90 genvar.o \
    geomlib.o linalg.o lockmod.o openmpmod.o

onetep_dmft_local_uniform.out: onetep_dmft_local_uniform.f90 genvar.o geomlib.o linalg.o lockmod.o \
    openmpmod.o

onetep_dmft_movie.out: onetep_dmft_movie.f90 genvar.o geomlib.o linalg.o lockmod.o openmpmod.o \
    stringmanip.o

onetep_dmft_one_iter.out: onetep_dmft_one_iter.f90 common_def.o genvar.o geomlib.o linalg.o \
    lockmod.o matrix.o namelistmod.o openmpmod.o random.o sorting.o specialfunction.o string4.o \
    string5.o stringmanip.o strings.o tools_algebra.o

onetep_dmft_proj_green_func_matrix.out: onetep_dmft_proj_green_func_matrix.f90 common_def.o \
    derivative_noise.o fourier_transform_mod.o genvar.o geometry.o geomlib.o linalg.o lockmod.o \
    matrix.o mesh.o mpirout.o openmpmod.o random.o smooth_data.o sorting.o splines.o splines2.o \
    stringmanip.o tools_algebra.o tools_fit.o

onetep_dmft_read_n.out: onetep_dmft_read_n.f90 common_def.o genvar.o geomlib.o linalg.o lockmod.o \
    matrix.o openmpmod.o random.o sorting.o tools_algebra.o

onetep_dmft_read_overlap_matrices.out: onetep_dmft_read_overlap_matrices.f90 common_def.o genvar.o \
    geomlib.o linalg.o lockmod.o matrix.o openmpmod.o random.o sorting.o tools_algebra.o

onetep_dmft_read_spin.out: onetep_dmft_read_spin.f90

onetep_dmft_sc.out: onetep_dmft_sc.f90 common_def.o genvar.o geomlib.o init_and_close_my_sim.o \
    linalg.o lockmod.o namelistmod.o openmpmod.o random.o string4.o string5.o stringmanip.o \
    strings.o

onetep_dmft_stitch.out: onetep_dmft_stitch.f90 common_def.o genvar.o geomlib.o linalg.o lockmod.o \
    namelistmod.o openmpmod.o string4.o string5.o stringmanip.o

onetep_dmft_summarize_sim.out: onetep_dmft_summarize_sim.f90 common_def.o derivative_noise.o \
    fourier_transform_mod.o genvar.o geometry.o geomlib.o init_and_close_my_sim.o linalg.o \
    lockmod.o matrix.o mesh.o mpirout.o namelistmod.o openmpmod.o random.o smooth_data.o sorting.o \
    splines.o splines2.o string4.o string5.o stringmanip.o tools_algebra.o tools_fit.o

onetep_dmft_trace_energy.out: onetep_dmft_trace_energy.f90

onetep_split.out: onetep_split.f90 common_def.o genvar.o geomlib.o init_and_close_my_sim.o \
    linalg.o lockmod.o mpirout.o namelistmod.o openmpmod.o random.o string4.o string5.o \
    stringmanip.o strings.o

onetep_split_serial.out: onetep_split_serial.f90 common_def.o genvar.o geomlib.o \
    init_and_close_my_sim.o linalg.o lockmod.o mpirout.o namelistmod.o openmpmod.o random.o \
    string4.o string5.o stringmanip.o strings.o

optical_conductivity_units.out: optical_conductivity_units.f90

plot_sigma.out: plot_sigma.f90

read_sigma.out: read_sigma.f90

self_energy_onetep.out: self_energy_onetep.f90

test_units_for_optics.out: test_units_for_optics.f90

variables_heterostructure.out: variables_heterostructure.f90 common_def.o genvar.o lockmod.o \
    namelistmod.o openmpmod.o string4.o string5.o stringmanip.o

