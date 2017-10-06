
   SHELL   =  /bin/bash
   FCOMPL  =  ifort
#   FFLAG   =   -O1 -arch i386 -parallel -openmp -static
#   FCOMPL  =  mpif90 
   FFLAG   =   -O1 #-m64  
#  FFLAG   =  -debug all -g -fpe0  -O0 -static -assume nobuffered_io -traceback -check all -check noarg_temp_created  

   RANLIB  =  ranlib


##########################################################################################
 ROUTINES = bispev.f clocur.f cocosp.f concon.f concur.f cualde.f curev.f curfit.f dblint.f evapol.f fourco.f fpader.f fpadno.f fpadpo.f fpback.f fpbacp.f fpbfout.f fpbisp.f fpbspl.f fpchec.f fpched.f fpchep.f fpclos.f fpcoco.f fpcons.f fpcosp.f fpcsin.f fpcurf.f fpcuro.f fpcyt1.f fpcyt2.f fpdeno.f fpdisc.f fpfrno.f fpgivs.f fpgrdi.f fpgrpa.f fpgrre.f fpgrsp.f fpinst.f fpintb.f fpknot.f fpopdi.f fpopsp.f fporde.f fppara.f fppasu.f fpperi.f fppocu.f fppogr.f fppola.f fprank.f fprati.f fpregr.f fprota.f fprppo.f fprpsp.f fpseno.f fpspgr.f fpsphe.f fpsuev.f fpsurf.f fpsysy.f fptrnp.f fptrpe.f insert.f parcur.f parder.f parsur.f percur.f pogrid.f polar.f profil.f regrid.f spalde.f spgrid.f sphere.f splder.f splev.f splint.f sproot.f surev.f surfit.f 
##########################################################################################

#-----------------------------------------------------------------------
# Target "all" makes everything
#-----------------------------------------------------------------------
all:    setup lib prog
	@echo ' *********** Finished compilation **************'
	@echo '       use make deltemp to remove bad files     '
	@echo ' ********************************************** '
	
#-----------------------------------------------------------------------
# Rules for compiling the .o files
#-----------------------------------------------------------------------
rout:  ${ROUTINES}
	@rm *.mod *.o *.a || echo "file first erased..."
	$(FCOMPL) -c $(FFLAG) ${ROUTINES}
	@echo 'compiled' 
	@echo 'set of files...'

build: rout lib
#-----------------------------------------------------------------------
# Target "lib" is used to built the subroutine library.
#-----------------------------------------------------------------------
lib: libdo 

libdo : 
	ar ru libdierick.a `ls *.o | sort | uniq `
	${RANLIB} libdierick.a
	@echo 'LIBRARY IS DONE.........'
	@rm *.o
	@echo ' '
	@echo ' '

#-----------------------------------------------------------------------
# Target "clean" is used to remove all the intermediate files.
#-----------------------------------------------------------------------
delfile: 
	@rm -f *.o


