
   SHELL              =    /bin/bash 
   RANLIB             =     ranlib 
   CUDA               =     
   CXX                =     icpc -c 
   NVCC               =     
   LINK               =     icpc  

   MAGMADEF=
   INCLUDES          := -I. 
   LIB               :=  

   C++        =  icpc -Wno-deprecated  
   FC         =  ifort -cxxlib -Wno-deprecated    

   CCFLAGS    =  -O1 
   CXXFLAGS   =  -O1 
   FCFLAGS    =  -O1 -fpp 

###############################
##     CUDA NOT INSTALLED    ##
##  UNCOMMENT THE FOLLOWING  ##
###############################

  FCFILES    =  fortran_cuda_routines_no_cuda_installed.f90 fortran_cuda_library.f90
  CXXFILES   =  no_cuda_installed.c

rout:  
	@rm ./MODFILES/*.mod  || echo "no such files..."
	@rm *.mod || echo "no such files..."
	@rm *.o || echo "no such files..."
	@rm *.a || echo "no such files..."
	${FC}                -c ${FCFLAGS}  ${FCFILES} || echo "no files.."
	${C++}  $(INCLUDES)  -c ${CXXFLAGS} ${CXXFILES} || echo "no files.."
	@mv *.mod MODFILES
	@rm *.linkinfo || echo "no such files..."
	@echo "done...."

build: rout libdo

libdo: 
	@echo "build library now...."
	ar ru libgpu.a `ls *.o | sort | uniq `
	${RANLIB} libgpu.a
	@rm *.o || echo "no such files...."
	@echo 'LIBRARY IS DONE.........'

clean :
	rm -f *.o *.linkinfo 
