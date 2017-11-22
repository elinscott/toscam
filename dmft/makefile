
#########################################################################
  HOME    =  /u/fs1/ebl27/scratch/toscam_full
  SHELL   =  /bin/bash
  FCOMPL  =  mpifort 
  FFLAG   =  -O1 -qopenmp -cxxlib -fpp   
#########################################################################
MKL              = -lmkl 
DIERICK          = -L$(HOME)/genlib/LIB_DIERICK/ -ldierick -L$(HOME)/genlib/LIB_FIELD_SPLINE -lfield  
SLATEC           = -L$(HOME)/genlib/LIBSLATEC/ -lslatec 
GPU              = -L${HOME}/genlib/LIB_GPU_NUMERICS -lgpu 
CUDA             =  
#########################################################################
ALL_LIB          =  ${MKL} $(DIERICK) $(SLATEC) #${X11}
#########################################################################

 LIBSOLVER   = -L/${HOME}/ed_solver -ldmftsolverED   
 LIBSOLVER_C = -L/${HOME}/ed_solver -ldmftsolverED_complex   
 LIBUTILS    = -L/${HOME}/genlib/LIBUTILS -lutils  

#########################################################################
 
 MODFILES    = -I$(HOME)/ed_solver/MODFILES -I$(HOME)/genlib/LIBUTILS/MODFILES  -I$(HOME)/genlib/LIB_GPU_NUMERICS/MODFILES  -I${HOME}/genlib/LIB_GPU_NUMERICS/MODFILES    

 MODFILES_C  = -I$(HOME)/ed_solver/MODFILES_COMPLEX -I$(HOME)/genlib/LIBUTILS/MODFILES  -I$(HOME)/genlib/LIB_GPU_NUMERICS/MODFILES -I${HOME}/LIB_GPU_NUMERICS/MODFILES

#########################################################################

 STATIC          =  ${LIBDMFT}$(LIBSOLVER)$(LIBPRO)${LIBUTILS}${GPU}${ALL_LIB} 
 STATIC_C        =  ${LIBDMFT}$(LIBSOLVER_C)$(LIBPRO)${LIBUTILS}${GPU}${ALL_LIB}
 FILES           =  dmft_one_iteration.f90
 FILES2          =  dmft_all_iterations.f90
 FILES3          =  dmft_standalone_ed.f90
 FILES4          =  dmft_one_iteration_split.f90
 FILES5          =  dmft_check_convergence.f90
 FILES6          =  dmft_check_cluster_convergence.f90 
 FILES7          =  dmft_collect_greens.f90 
 FILES8          =  dmft_project_orbitals.f90  
 FILES9          =  dmft_project_orbitals_back.f90
 FILES10         =  dmft_group_dimers.f90
 FILES11         =  dmft_generate_edcorrel.f90
 FILES12         =  dmft_project_orbitals_sigma.f90
 FILES13         =  duplicate_unitcell.f90
 FILES14         =  build_local_rotations_along_octahedra.f90  
 FILES15         =  dmft_collect_dos.f90
 FILES16         =  build_dimers.f90
 FILES17         =  build_uniform_mask.f90
 FILES18         =  dmft_one_iteration_split_serial.f90
 FILES19         =  onetep.dmft.density.matsu.f90
 FILES20         =  onetep.dmft.summarize.sim.f90
 FILES21         =  onetep.dmft.compute.optics.f90
 FILES22         =  dmft_all_iterations.onetep.f90
 FILES23         =  onetep.dmft.proj.green.func.matrix.f90
 FILES24         =  dmft_standalone_ed_serial.f90
 FILES25         =  build_local_rotations_for_one_atom.f90
 FILES26         =  onetep.dmft.find.rot.for.d.orbital.transform.f90
 FILES27         =  onetep.dmft.do.time.evolution.f90
 FILES28         =  onetep.dmft.read.overlap.matrices.f90
 FILES29         =  onetep.dmft.correct.spin.f90
 FILES30         =  onetep.dmft.histogram.f90
 FILES31         =  dmft_check_cluster_convergence_renom.f90
 FILES32         =  dmft_check_convergence_renorm.f90
 FILES33         =  onetep.dmft.trace.energy.f90
 FILES34         =  onetep.dmft.correct.n.f90
 FILES35         =  onetep.dmft.build.optical.tensor.f90
 FILES36         =  onetep.dmft.calculate.average.e.f90
 FILES37         =  dmft_onetep_sc.f90
 FILES38         =  onetep.dmft.hf.f90
 FILES39         =  onetep.dmft.build.heterostructure.f90
 FILES40         =  onetep.dmft.stitch.f90
 FILES41         =  dmft_collect_sigma_embed.f90
 FILES42         =  onetep.dmft.embedding.gen.sigma.f90
 FILES43         =  dmft_generate_edcorrel.outputfile.f90
#########################################################################

prog: ${FILES}
	@echo " ........... compile : static / dynamic linking ........... "
	${FCOMPL} -o ed_solver.out                           ${FFLAG} ${FILES3}  ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep_dmft_one_iter.out                ${FFLAG} ${FILES}   ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o dmft_gen_correl.output.file.out         ${FFLAG} ${FILES43} ${MODFILES} ${STATIC} ${CUDA}       
	${FCOMPL} -o onetep.dmft.embedding.gen.sigma.out     ${FFLAG} ${FILES42} ${MODFILES} ${STATIC} ${CUDA}  
	${FCOMPL} -o dmft_collect_sigma_embed.out            ${FFLAG} ${FILES41} ${MODFILES} ${STATIC} ${CUDA}  
	${FCOMPL} -o onetep.dmft.stitch.out                  ${FFLAG} ${FILES40} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.build.heterostructure.out   ${FFLAG} ${FILES39} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.sc                          ${FFLAG} ${FILES37} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.calculate.average.e.out     ${FFLAG} ${FILES36} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.build.optical.tensor.out    ${FFLAG} ${FILES35} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.read.n.out                  ${FFLAG} ${FILES34} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.local.rotation              ${FFLAG} ${FILES14} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.read.spin.out               ${FFLAG} ${FILES29} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.trace.energy.out            ${FFLAG} ${FILES33} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o check_convergence_full.out              ${FFLAG} ${FILES6}  ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o check_convergence_single_renorm.out     ${FFLAG} ${FILES32} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o check_convergence_renorm.out            ${FFLAG} ${FILES31} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.histogram.out               ${FFLAG} ${FILES30} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o check_convergence.out                   ${FFLAG} ${FILES5}  ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.read.overlap.matrices.out   ${FFLAG} ${FILES28} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o ed_solver_serial.out                    ${FFLAG} ${FILES24} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.movie                       ${FFLAG} ${FILES27} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.local.rotation.find         ${FFLAG} ${FILES26} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.local.rotation.one.atom     ${FFLAG} ${FILES25} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.proj.green.func.matrix.out  ${FFLAG} ${FILES23} ${MODFILES} ${STATIC} ${CUDA}   
	${FCOMPL} -o onetep.dmft.local.dimers.out            ${FFLAG} ${FILES16} ${MODFILES} ${STATIC} ${CUDA}   
	${FCOMPL} -o dmft_all_iterations.onetep.out          ${FFLAG} ${FILES22} ${MODFILES} ${STATIC} ${CUDA}	
	${FCOMPL} -o onetep.dmft.compute.optics.out          ${FFLAG} ${FILES21} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.summarize.sim.out           ${FFLAG} ${FILES20} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.density.matsu.out           ${FFLAG} ${FILES19} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.local.uniform               ${FFLAG} ${FILES17} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep.dmft.duplicate.unitcell          ${FFLAG} ${FILES13} ${MODFILES} ${STATIC} ${CUDA}       
	${FCOMPL} -o onetep.dmft                             ${FFLAG} ${FILES2}  ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep_split.out                        ${FFLAG} ${FILES4}  ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o onetep_split_serial.out                 ${FFLAG} ${FILES18} ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o dmft_collect.out                        ${FFLAG} ${FILES7}  ${MODFILES} ${STATIC} ${CUDA}
	${FCOMPL} -o dmft_project.out                        ${FFLAG} ${FILES8}  ${MODFILES} ${STATIC} ${CUDA}            
	${FCOMPL} -o dmft_project_back.out                   ${FFLAG} ${FILES9}  ${MODFILES} ${STATIC} ${CUDA}            
	${FCOMPL} -o dmft_group_dimers.out                   ${FFLAG} ${FILES10} ${MODFILES} ${STATIC} ${CUDA}       
	${FCOMPL} -o dmft_gen_correl.out                     ${FFLAG} ${FILES11} ${MODFILES} ${STATIC} ${CUDA}       
	${FCOMPL} -o dmft_project_sigma.out                  ${FFLAG} ${FILES12} ${MODFILES} ${STATIC} ${CUDA}       
	${FCOMPL} -o dmft_collect_dos.out                    ${FFLAG} ${FILES15} ${MODFILES} ${STATIC} ${CUDA}

#	${FCOMPL} -o ed_solver_complex.out                   ${FFLAG} ${FILES3}  ${MODFILES_C} ${STATIC_C} ${CUDA}
#############################
#############################
#############################

	@echo " ...................... done .............................. "

#########################################################################

