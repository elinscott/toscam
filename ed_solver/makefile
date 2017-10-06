
  SHELL  = /bin/bash
  FCOMPL = mpif90   

  prof  = 
#-p -g

  prep   = -cpp -D_complex   
  FFLAG  = -heap-arrays 30000 -O1 -ip -arch i386 -parallel -DOPENMP_MPI -openmp -std03 >& error_output   
# FFLAG  = -arch i386 -openmp -debug all -g -fpe0 -O0  -static -assume nobuffered_io -traceback -check all -check noarg_temp_created -std03 >& error_output  
  FFLAG  = -O1 -fopenmp -cpp -DOPENMP_MPI -ffree-line-length-4000 >& error_output
  
#-DOPENMP_MPI -DNO_SYS_CALL -DOPENMP_MPI_SAFE for TOSCAM

  FFLAG2 = -O0 -arch i386 -heap-arrays  
  FFLAG2 = -O0  -ffree-line-length-4000 
  RANLIB = ranlib

ROUTINES = Global_Variables.f90 masked_matrix_real_comp.f90 rcvector_class.f90 rcmatrix_class.f90 eigen_class.f90 eigen_sector_class.f90 impurity_class.f90 correl_class.f90 bath_class.f90 bath_class_vec.f90 bath_class_hybrid.f90 quantum_algebra_Anderson_Impurity_class.f90 quantum_algebra_Anderson_Impurity_class2.f90 HAIM_offdiag_arrays.f90 HAIMsz.f90 HAIMupdo.f90 density_matrix.f90 overlap.f90 apply_C.f90 apply_NS.f90 apply_P.f90 H_class.f90 Lanczos_gpu.f90 Lanczos_plain.f90 ED_ARPACK.f90 Block_Lanczos.f90 green_class.f90 green_class_compute_dynamic.f90 green_class_compute_symmetric_combination.f90 green_class_computeAA.f90 green_class_computeAB.f90 correlations.f90 Lanczos_Cullum_mult.f90 Lanczos_Cullum_lesub.f90 Lanczos_Cullum_hleval.f90 Lanczos_Cullum_leval.f90 Lanczos_Cullum_hlevec.f90 Lanczos_Cullum_levec.f90 Lanczos_Cullum_wrapper.f90 Chiloc_four_leg.f90 solver.f90 dmft_ed_solver.f90

ROUTINES2 = eigen_sector_class.f90

MODFILES=-I$(HOME)/MYLIBRARY/LIBUTILS/MODFILES -I./MODFILES -I$(HOME)/MYLIBRARY/LIB_GPU_NUMERICS/MODFILES -I$(HOME)/MYLIBRARY/LIB_GPU_NUMERICS_MAGMA/include -I$(HOME)/MYLIBRARY/LIB_GPU_NUMERICS_MAGMA/control   

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






