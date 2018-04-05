# Second-level makefile for TOSCAM splines
#
# System-specific changes should be made to the appropriate
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

VPATH = $(ROOTDIR)/lib/splines:$(ROOTDIR)/lib/splines/obj

#
# Module redirection
#

ifeq ($(F90),$(filter $(F90),ifort mpifort))
   REDIRECT = -module $(ROOTDIR)/lib/splines/obj
endif
ifeq ($(F90),$(filter $(F90),gfortran, mpif90))
   REDIRECT = -J$(ROOTDIR)/lib/splines/obj
endif

#
# List of source files
#

ROUTINES = bispev.f clocur.f cocosp.f concon.f concur.f cualde.f curev.f \
           curfit.f dblint.f evapol.f fourco.f fpader.f fpadno.f fpadpo.f \
           fpback.f fpbacp.f fpbfout.f fpbisp.f fpbspl.f fpchec.f fpched.f \
           fpchep.f fpclos.f fpcoco.f fpcons.f fpcosp.f fpcsin.f fpcurf.f \
           fpcuro.f fpcyt1.f fpcyt2.f fpdeno.f fpdisc.f fpfrno.f fpgivs.f \
           fpgrdi.f fpgrpa.f fpgrre.f fpgrsp.f fpinst.f fpintb.f fpknot.f \
           fpopdi.f fpopsp.f fporde.f fppara.f fppasu.f fpperi.f fppocu.f \
           fppogr.f fppola.f fprank.f fprati.f fpregr.f fprota.f fprppo.f \
           fprpsp.f fpseno.f fpspgr.f fpsphe.f fpsuev.f fpsurf.f fpsysy.f \
           fptrnp.f fptrpe.f insert.f parcur.f parder.f parsur.f percur.f \
           pogrid.f polar.f profil.f regrid.f spalde.f spgrid.f sphere.f \
           splder.f splev.f splint.f sproot.f surev.f surfit.f

OUTPUT_SOURCES= $(filter %.f,$(ROUTINES))
OUTPUTS = $(OUTPUT_SOURCES:.f=.o) \

#
# Specific compilation rules
#

default: build

build: routines lib

routines: $(OUTPUTS)

lib: 
	@( cd $(ROOTDIR)/lib/splines/obj ; \
	ar r libsplines.a `ls *.o | sort | uniq `; \
	$(RANLIB) libsplines.a)

clean: 
	$(RM) obj/*

#
# General rules
#

%.o: %.f
	$(F90) -c $(FFLAGS) $< -o $(ROOTDIR)/lib/splines/obj/$@ $(REDIRECT)

#
# Dependencies
#

bispev.o: bispev.f
clocur.o: clocur.f
cocosp.o: cocosp.f
concon.o: concon.f
concur.o: concur.f
cualde.o: cualde.f
curev.o: curev.f
curfit.o: curfit.f
dblint.o: dblint.f
evapol.o: evapol.f
fourco.o: fourco.f
fpader.o: fpader.f
fpadno.o: fpadno.f
fpadpo.o: fpadpo.f
fpback.o: fpback.f
fpbacp.o: fpbacp.f
fpbfout.o: fpbfout.f
fpbisp.o: fpbisp.f
fpbspl.o: fpbspl.f
fpchec.o: fpchec.f
fpched.o: fpched.f
fpchep.o: fpchep.f
fpclos.o: fpclos.f
fpcoco.o: fpcoco.f
fpcons.o: fpcons.f
fpcosp.o: fpcosp.f
fpcsin.o: fpcsin.f
fpcurf.o: fpcurf.f
fpcuro.o: fpcuro.f
fpcyt1.o: fpcyt1.f
fpcyt2.o: fpcyt2.f
fpdeno.o: fpdeno.f
fpdisc.o: fpdisc.f
fpfrno.o: fpfrno.f
fpgivs.o: fpgivs.f
fpgrdi.o: fpgrdi.f
fpgrpa.o: fpgrpa.f
fpgrre.o: fpgrre.f
fpgrsp.o: fpgrsp.f
fpinst.o: fpinst.f
fpintb.o: fpintb.f
fpknot.o: fpknot.f
fpopdi.o: fpopdi.f
fpopsp.o: fpopsp.f
fporde.o: fporde.f
fppara.o: fppara.f
fppasu.o: fppasu.f
fpperi.o: fpperi.f
fppocu.o: fppocu.f
fppogr.o: fppogr.f
fppola.o: fppola.f
fprank.o: fprank.f
fprati.o: fprati.f
fpregr.o: fpregr.f
fprota.o: fprota.f
fprppo.o: fprppo.f
fprpsp.o: fprpsp.f
fpseno.o: fpseno.f
fpspgr.o: fpspgr.f
fpsphe.o: fpsphe.f
fpsuev.o: fpsuev.f
fpsurf.o: fpsurf.f
fpsysy.o: fpsysy.f
fptrnp.o: fptrnp.f
fptrpe.o: fptrpe.f
insert.o: insert.f
parcur.o: parcur.f
parder.o: parder.f
parsur.o: parsur.f
percur.o: percur.f
pogrid.o: pogrid.f
polar.o: polar.f
profil.o: profil.f
regrid.o: regrid.f
spalde.o: spalde.f
spgrid.o: spgrid.f
sphere.o: sphere.f
splder.o: splder.f
splev.o: splev.f
splint.o: splint.f
sproot.o: sproot.f
surev.o: surev.f
surfit.o: surfit.f
