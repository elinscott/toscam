# Second-level makefile for TOSCAM exact diagonalisation solver
#
# Do not edit; system-specific changes should be made to the appropriate
#  conf.arch files
#

ifndef ROOTDIR
   ROOTDIR=../..
endif
ifdef ARCH
   include $(ROOTDIR)/config/conf.$(ARCH)
else
   include $(ROOTDIR)/config/conf.default
endif

#
# Facilitating execution from a different directory
#

VPATH=$(ROOTDIR)/src/ed_solver:$(ROOTDIR)/src/ed_solver/obj

#
# Compiler-dependent flags
#

ifeq ($(COMPILER),$(filter $(COMPILER),intel))
   REDIRECT = -module $(ROOTDIR)/src/ed_solver/obj
   prep = -fpp -D__complex
endif
ifeq ($(COMPILER),$(filter $(COMPILER),gnu))
   REDIRECT = -J$(ROOTDIR)/src/ed_solver/obj
   prep = -cpp -D__complex
endif
ifeq ($(COMPILER),$(filter $(COMPILER),nag))
   REDIRECT = -mdir $(ROOTDIR)/src/ed_solver/obj
endif
ifneq ($(COMPILER), $(filter $(COMPILER),intel gnu nag))
   $(error Please specify compiler type using COMPILER entry in config file (intel/gnu/nag))
endif

#
# List of source files
#

ROUTINES = aim2_class.f90 \
           aim_class.f90 \
           apply_c.f90 \
           apply_ns.f90 \
           apply_p.f90 \
           bath_class.f90 \
           bath_class_hybrid.f90 \
           bath_class_vec.f90 \
           block_lanczos.f90 \
           correl_class.f90 \
           correlations.f90 \
           density_matrix.f90 \
           dmft_solver_ed.f90 \
           ed_arpack.f90 \
           eigen_class.f90 \
           eigen_sector_class.f90 \
           fermion_hilbert_class.f90 \
           fermion_ket_class.f90 \
           fermion_sector2_class.f90 \
           frequency_class.f90 \
           globalvar_ed_solver.f90 \
           green_class.f90 \
           green_class_compute_dynamic.f90 \
           green_class_compute_symmetric.f90 \
           green_class_computeaa.f90 \
           green_class_computeab.f90 \
           h_class.f90 \
           haim2_class.f90 \
           haimsz_class.f90 \
           haimupdo_class.f90 \
           impurity_class.f90 \
           lanczos_fast.f90 \
           mask_class.f90 \
           masked_matrix_class.f90 \
           masked_matrix_class_mod.f90 \
           overlap_module.f90 \
           quantum_algebra.f90 \
           rcmatrix_class.f90 \
           rcvector_class.f90 \
           readable_vec_class.f90 \
           sector_class.f90 \
           solver.f90 \
           vertex.f90

MODULE_SOURCES= $(filter %.f90,$(ROUTINES))
MODULES = $(MODULE_SOURCES:.f90=.mod)

#
# Libraries
#

LIBSPLINES = -L$(ROOTDIR)/src/splines/obj -lsplines
LIBSLATEC  = -L$(ROOTDIR)/src/slatec/static -lslatec
LIBS       = $(LIBSPLINES)

#
# Include
#

INCUTILS = -I$(ROOTDIR)/src/utils/obj
INC      = $(INCUTILS)

#  
# Specific compilation rules
#

default: build

build: $(MODULES)

lib: 

	ar ru libdmftsolverED.a `ls *.o | sort | uniq`
	${RANLIB} libdmftsolverED.a

clean:
	$(RM) obj/*

#
# General rules
#

%.mod: %.f90
	$(F90) -c $(prof) $(FFLAGS) $(prep) $< $(MKL) $(LIBS) $(INC) -o $(ROOTDIR)/src/ed_solver/obj/$(@:.mod=.o) $(REDIRECT)

#
# Dependencies
#

aim2_class.mod: aim2_class.f90 aim_class.mod bath_class.mod impurity_class.mod masked_matrix_class.mod

aim_class.mod: aim_class.f90 bath_class.mod globalvar_ed_solver.mod impurity_class.mod

apply_c.mod: apply_c.f90 aim_class.mod eigen_class.mod eigen_sector_class.mod fermion_hilbert_class.mod fermion_ket_class.mod fermion_sector2_class.mod globalvar_ed_solver.mod sector_class.mod

apply_ns.mod: apply_ns.f90 aim_class.mod eigen_class.mod eigen_sector_class.mod fermion_hilbert_class.mod fermion_ket_class.mod fermion_sector2_class.mod sector_class.mod

apply_p.mod: apply_p.f90 aim_class.mod eigen_class.mod eigen_sector_class.mod fermion_ket_class.mod fermion_sector2_class.mod sector_class.mod

bath_class.mod: bath_class.f90 correl_class.mod globalvar_ed_solver.mod masked_matrix_class.mod

bath_class_hybrid.mod: bath_class_hybrid.f90 bath_class.mod bath_class_vec.mod correl_class.mod globalvar_ed_solver.mod impurity_class.mod masked_matrix_class.mod

bath_class_vec.mod: bath_class_vec.f90 bath_class.mod globalvar_ed_solver.mod masked_matrix_class.mod masked_matrix_class_mod.mod

block_lanczos.mod: block_lanczos.f90 eigen_class.mod globalvar_ed_solver.mod h_class.mod

correl_class.mod: correl_class.f90 frequency_class.mod globalvar_ed_solver.mod masked_matrix_class_mod.mod

correlations.mod: correlations.f90 aim_class.mod apply_c.mod apply_ns.mod apply_p.mod bath_class_hybrid.mod correl_class.mod density_matrix.mod eigen_sector_class.mod globalvar_ed_solver.mod green_class.mod green_class_computeaa.mod green_class_computeab.mod impurity_class.mod mask_class.mod masked_matrix_class.mod overlap_module.mod rcmatrix_class.mod rcvector_class.mod readable_vec_class.mod

density_matrix.mod: density_matrix.f90 aim_class.mod eigen_class.mod eigen_sector_class.mod globalvar_ed_solver.mod rcmatrix_class.mod readable_vec_class.mod sector_class.mod

dmft_solver_ed.mod: dmft_solver_ed.f90 aim_class.mod bath_class.mod bath_class_hybrid.mod correl_class.mod correlations.mod eigen_sector_class.mod globalvar_ed_solver.mod impurity_class.mod solver.mod

ed_arpack.mod: ed_arpack.f90 eigen_class.mod globalvar_ed_solver.mod h_class.mod rcvector_class.mod

eigen_class.mod: eigen_class.f90 globalvar_ed_solver.mod rcvector_class.mod

eigen_sector_class.mod: eigen_sector_class.f90 eigen_class.mod globalvar_ed_solver.mod sector_class.mod

fermion_hilbert_class.mod: fermion_hilbert_class.f90 fermion_ket_class.mod quantum_algebra.mod

fermion_ket_class.mod: fermion_ket_class.f90

fermion_sector2_class.mod: fermion_sector2_class.f90 fermion_hilbert_class.mod quantum_algebra.mod

frequency_class.mod: frequency_class.f90

globalvar_ed_solver.mod: globalvar_ed_solver.f90

green_class.mod: green_class.f90 correl_class.mod masked_matrix_class.mod

green_class_compute_dynamic.mod: green_class_compute_dynamic.f90 eigen_class.mod eigen_sector_class.mod frequency_class.mod globalvar_ed_solver.mod h_class.mod lanczos_fast.mod rcvector_class.mod

green_class_compute_symmetric.mod: green_class_compute_symmetric.f90 eigen_class.mod eigen_sector_class.mod sector_class.mod

green_class_computeaa.mod: green_class_computeaa.f90 aim_class.mod eigen_class.mod eigen_sector_class.mod globalvar_ed_solver.mod green_class.mod green_class_compute_dynamic.mod green_class_compute_symmetric.mod h_class.mod mask_class.mod sector_class.mod

green_class_computeab.mod: green_class_computeab.f90 aim_class.mod eigen_class.mod eigen_sector_class.mod globalvar_ed_solver.mod green_class.mod green_class_compute_symmetric.mod mask_class.mod sector_class.mod

h_class.mod: h_class.f90 aim_class.mod globalvar_ed_solver.mod haimsz_class.mod haimupdo_class.mod sector_class.mod

haim2_class.mod: haim2_class.f90 aim2_class.mod aim_class.mod fermion_hilbert_class.mod fermion_ket_class.mod globalvar_ed_solver.mod impurity_class.mod masked_matrix_class.mod

haimsz_class.mod: haimsz_class.f90 aim_class.mod fermion_hilbert_class.mod fermion_ket_class.mod globalvar_ed_solver.mod haim2_class.mod masked_matrix_class_mod.mod

haimupdo_class.mod: haimupdo_class.f90 aim_class.mod fermion_hilbert_class.mod fermion_ket_class.mod fermion_sector2_class.mod globalvar_ed_solver.mod haim2_class.mod masked_matrix_class_mod.mod

impurity_class.mod: impurity_class.f90 globalvar_ed_solver.mod masked_matrix_class.mod masked_matrix_class_mod.mod

lanczos_fast.mod: lanczos_fast.f90 eigen_class.mod globalvar_ed_solver.mod h_class.mod rcvector_class.mod

mask_class.mod: mask_class.f90

masked_matrix_class.mod: masked_matrix_class.f90 globalvar_ed_solver.mod masked_matrix_class_mod.mod

masked_matrix_class_mod.mod: masked_matrix_class_mod.f90 mask_class.mod

overlap_module.mod: overlap_module.f90 aim_class.mod eigen_class.mod eigen_sector_class.mod readable_vec_class.mod sector_class.mod

quantum_algebra.mod: quantum_algebra.f90

rcmatrix_class.mod: rcmatrix_class.f90

rcvector_class.mod: rcvector_class.f90 globalvar_ed_solver.mod

readable_vec_class.mod: readable_vec_class.f90 fermion_ket_class.mod

sector_class.mod: sector_class.f90 fermion_hilbert_class.mod fermion_sector2_class.mod

solver.mod: solver.f90 aim_class.mod eigen_sector_class.mod ed_arpack.mod globalvar_ed_solver.mod vertex.mod

vertex.mod: vertex.f90 aim_class.mod apply_c.mod eigen_class.mod eigen_sector_class.mod globalvar_ed_solver.mod sector_class.mod

