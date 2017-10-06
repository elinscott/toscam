
   HOME=/rscratch/ebl27/MYLIBRARY_ONETEP/
   SHELL    =  /bin/bash 
   RANLIB   =  ranlib 

   prep     = -fpp -D__DUMP_DISLIN  

   FCOMPL   =  mpifort  
   FFLAG    =  -O1 -openmp -DNOSVD     

   CCOMPL   =  mpicc     
   CFLAG    =  -O1 

   C++      =  mpicpc   
   CXXFLAG  =  -O1  

ROUTINES=core_functions_fortran.f90 test_memory_limit.f90 openmpmod.f90 GeneralVariables.f90 tools_arbitrary_precision_module.f90 tools_arbitrary_precision_interfaces.f90 tools_arbitrary_precision.f90 fourier_transform.f90 core_functions_def.f90 conj_grad_civelli.f90 geometry_lib.f90 tools_nms_lib.f90 read_scs_file.f90 read_structure_constant.f90 rotate_overlap_unity.f90 geometry2.f90 tetraedron_method_tool.f90 tetraedron_method.f90 spheric_harmonics.f90 Strings2.f90 Strings3.f90 Strings4.f90 Strings5.f90 tools_fit.f90 FuncAnalysis.f90 mixing_tool.f90 linalg.f90 units.f90 fast_fourier_mkl_dfti.f90 fast_fourier1.f90 fast_fourier2.f90 fast_fourier.f90 tools_algebra.f90 random.f90 MPI_routines.f90 metropolis_tool.f90 band_structure_complex.f90 sorting.f90 matrices.f90 matrices_tridiag.f90 matrices_analytical.f90 derivative_noisy_data.f90 Combinatory.f90 Strings.f90 liste.f90 liste3.f90 specFunc.f90 hartree_fock.f90 rotation_matrix.f90 geometry.f90 smooth_data.f90 splines2.f90 splines.f90 files_handling.f90 symmetrize_tool.f90 tools_fit.f90 mesh.f90 integrals.f90 namelist_module.f90 pointers_smart.f90 kramers_kronig.f90 dft_lib.f90 tools_minimization_wrapping.f90 core_functions_initsim_closesim.f90  drive_quantity_toward_value.f90

ROUTINESC=IEEE.c 

ROUTINESCXX=

MODFILES=-I$(HOME)/MYLIBRARY/LIB_GPU_NUMERICS/MODFILES 

###########################################################################################################

all:    
	@echo '*********** Finished compilation **************'
	@echo '       use make deltemp to remove bad files    '
	@echo '***********************************************'
rout:  
	@rm ./MODFILES/*.mod  || echo "no such files..."
	@rm *.mod || echo "no such files..."
	@rm *.o || echo "no such files..."
	@rm *.a || echo "no such files..."
	$(FCOMPL) -c ${prep} $(FFLAG) ${ROUTINES} ${MODFILES}
	${CCOMPL} -c ${CFLAG} ${prof} ${ROUTINESC} 
	@echo "done...."
	@echo "now compile C++ files..." 
	$(C++) -c $(CXXFLAG) memory_address_functions.cc 
	@mv *.mod MODFILES
	@rm *.linkinfo || echo "no such files..."

build: rout libdo

libdo: 
	@echo "build library now...."
	ar ru libutils.a `ls *.o | sort | uniq `
	${RANLIB} libutils.a
	@rm *.o || echo "no such files...."
	@echo 'LIBRARY IS DONE.........'

.SUFFIXES : .cc
.cc.o:
	$(C++) $(CFLAGS) -c $<

