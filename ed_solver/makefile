  HOME=/scratch/ebl27/toscam_full
  SHELL  = /bin/bash
  FCOMPL = mpifort  
  
  prof  = 
#-p -g

  prep   = -fpp -D__complex   
 #FFLAG  = -O0 -traceback -g -std03 >& error_output  
  FFLAG  = -O1 -openmp -std03 >& error_output

# FFLAG +=  -DNO_SYS_CALL -DOPENMP_MPI_SAFE -g -traceback   
#-DOPENMP_MPI -DNO_SYS_CALL -DOPENMP_MPI_SAFE for TOSCAM

  FFLAG2 = -O0   
  RANLIB = ranlib

ROUTINES =  mask_class.f90 masked_mat_class.f90 quantum_algebra.f90 quantum_fermion_ket_class.f90 quantum_readable_vec.f90 quantum_fermion_hilbert_class.f90 quantum_fermion_sector_class.f90 quantum_sector_class.f90 frequency_axis.f90 Global_Variables.f90 masked_matrix_real_comp.f90 rcvector_class.f90 rcmatrix_class.f90 eigen_class.f90 eigen_sector_class.f90 impurity_class.f90 correl_class.f90 bath_class.f90 bath_class_vec.f90 bath_class_hybrid.f90 quantum_algebra_Anderson_Impurity_class.f90 quantum_algebra_Anderson_Impurity_class2.f90 HAIM_offdiag_arrays.f90 HAIMsz.f90 HAIMupdo.f90 density_matrix.f90 overlap.f90 apply_C.f90 apply_NS.f90 apply_P.f90 H_class.f90 Lanczos_gpu.f90 Lanczos_plain.f90 ED_ARPACK.f90 Block_Lanczos.f90 green_class.f90 green_class_compute_dynamic.f90 green_class_compute_symmetric_combination.f90 green_class_computeAA.f90 green_class_computeAB.f90 correlations.f90 Chiloc_four_leg.f90 solver.f90 dmft_ed_solver.f90

ROUTINES2 = eigen_sector_class.f90

MODFILES=-I$(HOME)/genlib/LIBUTILS/MODFILES -I./MODFILES -I$(HOME)/genlib/LIB_GPU_NUMERICS/MODFILES 

comp:  ${ROUTINES} 
	@rm *.o *.mod ./MODFILES/*.mod || echo 'files erased...'
	@rm libdmftsolverED.a || echo 'old lib erased...'
	@echo "compile"
	$(FCOMPL) -c ${prof} $(FFLAG) ${prep} ${ROUTINES} ${MODFILES} 
	@echo "set of files..."
	@mv *.mod MODFILES/

comp2:  ${ROUTINES2}
	@echo "compile2"
	$(FCOMPL) -c $(FFLAG2) ${prof} ${prep} ${ROUTINES2} ${MODFILES}
	@echo "set of files..."
	@mv *.mod MODFILES/

build: comp comp2 libdmftsolverED.a

libdmftsolverED.a: 
	ar ru $@ `ls *.o | sort | uniq`
	${RANLIB} $@
	@echo 'LIBRARY IS DONE.........'
	@rm *.o || echo 'files erased....'
	@echo ' '

