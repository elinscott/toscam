
   SHELL    =  /bin/bash 
   RANLIB   =  ranlib 
   FCOMPL   =  ifort  
   FFLAG    =  -O1 #-arch i386 -parallel -openmp -static   

   #FCOMPL  =  mpif90
   #FFLAG   =   -O1 -m64

###########################################################################################################
ROUTINES=css.f csstr.f cvrf.f dchold.f dcopy.f ddot.f dlv.f dmaket.f drdfun.f dsetup.f expbs.f expfn.f gaspbs.f gaspfn.f gcvcss.f gcvfc.f hsort.f ifind.f inpoly.f lscv.f m2deb.f mkpoly.f mltdrb.f mltdtd.f msort.f multeb.f multgb.f multrb.f radbas.f radfun.f rcss.f rcssr.f rcsswt.f sortm.f 
###########################################################################################################

default: build

all:    
	@echo '*********** Finished compilation **************'
	@echo '       use make deltemp to remove bad files    '
	@echo '***********************************************'
rout:  
	@rm *.a || echo "no such files..."
	$(FCOMPL) -c $(FFLAG) ${ROUTINES} 
	@echo "done...."

build: rout libdo

libdo: 
	@echo "build library now...."
	ar ru libfield.a `ls *.o | sort | uniq `
	${RANLIB} libfield.a
	#@rm *.o || echo "no such files...."
	@echo 'LIBRARY IS DONE.........'

clean:
	$(RM) *.o *.a


.SUFFIXES : .cc
.cc.o:
	$(C++) $(CFLAGS) -c $<

