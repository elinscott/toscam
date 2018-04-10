# Second-level makefile for TOSCAM utils
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

VPATH=$(ROOTDIR)/lib/utils/obj:$(ROOTDIR)/lib/utils

prep     = -fpp

LIBS     = 

#
# Module redirection
#

ifeq ($(F90),$(filter $(F90),ifort mpifort))
   REDIRECT = -module $(ROOTDIR)/lib/utils/obj
endif
ifeq ($(F90),$(filter $(F90),gfortran, mpif90))
   REDIRECT = -J$(ROOTDIR)/lib/utils/obj
endif

#
# List of source files
#

ROUTINES= common_def.f90 \
          conj_grad.f90 \
          derivative_noise.f90 \
          fourier_transform_f90.mod \
          genvar.f90 \
          geometry.f90 \
          geometry2.f90 \
          geomlib.f90 \
          init_and_close_my_sim.f90 \
          linalg.f90 \
          liste3.f90 \
          lockf90.mod \
          matrix.f90 \
          mesh.f90 \
          minimization_wrapping.f90 \
          mpirout.f90 \
          namelistf90.mod \
          openmpf90.mod \
          random.f90 \
          smooth_data.f90 \
          sorting.f90 \
          specialfunction.f90 \
          splines2.f90 \
          splines.f90 \
          stringmanip.f90 \
          strings.f90 \
          string4.f90 \
          string5.f90 \
          tools_algebra.f90 \
          tools_fit.f90 \
          tridiag_class.f90
ROUTINESC=IEEE.c
ROUTINESCXX=memory_address_functions.cc

MODULE_SOURCES= $(filter %.f90,$(ROUTINES))
MODULE_SOURCES_C= $(filter %.c,$(ROUTINESC))
MODULE_SOURCES_CXX= $(filter %.cc,$(ROUTINESCXX))

MODULES_F  = $(MODULE_SOURCES:.f90=.mod)
MODULES_C  = $(MODULE_SOURCES_C:.c=.o)
MODULES_CC = $(MODULE_SOURCES_CXX:.cc=.o) 

MODULES = $(MODULES_F) $(MODULES_C) $(MODULES_CC)

#
# Specific compilation rules
#

default: build

routines: $(MODULES) $(ROUTINESC) $(ROUTINESCXX)

build: routines lib

lib: 
	@( cd $(ROOTDIR)/lib/utils/obj ; \
	ar ru libutils.a $(MODULES_C) $(MODULES_CC); \
	$(RANLIB) libutils.a)

clean:
	$(RM) obj/*

#
# General rules
#

%.mod: %.f90
	$(F90) -c $(prep) $(FFLAGS) $(LIBS) $< -o $(ROOTDIR)/lib/utils/obj/$(@:.mod=.o) $(REDIRECT)

%.o: %.c
	$(CCOMPL) -c $(CFLAGS) $< -o $(ROOTDIR)/lib/utils/obj/$(@:.c=.o)

%.o: %.cc
	$(C++) -c $(CXXFLAGS) $< -o $(ROOTDIR)/lib/utils/obj/$(@:.cc=.o)

#
# Dependencies
#

conj_grad.mod: conj_grad.f90 common_def.mod genvar.mod
 
common_def.mod: common_def.f90 genvar.mod
 
derivative_noise.mod: derivative_noise.f90 fourier_transform_mod.mod \
                      matrix.mod
 
fourier_transform_mod.mod: fourier_transform_mod.f90 genvar.mod
 
genvar.mod: genvar.f90 openmpmod.mod
 
geometry.mod: geometry.f90 genvar.mod geomlib.mod linalg.mod matrix.mod \
              random.mod sorting.mod

geometry2.mod: geometry2.f90
 
geomlib.mod: geomlib.f90 genvar.mod
 
init_and_close_my_sim.mod: init_and_close_my_sim.f90 common_def.mod \
                           genvar.mod linalg.mod namelistmod.mod random.mod
 
linalg.mod: linalg.f90 genvar.mod geomlib.mod
 
liste3.mod: liste3.f90 common_def.mod genvar.mod

lockmod.mod: lockmod.f90
 
matrix.mod: matrix.f90 common_def.mod genvar.mod linalg.mod random.mod \
            sorting.mod tools_algebra.mod
 
mesh.mod: mesh.f90 genvar.mod geometry.mod linalg.mod sorting.mod splines.mod \
          stringmanip.mod
 
minimization_wrapping.mod: minimization_wrapping.f90 \
                           conj_grad.mod genvar.mod linalg.mod \
                           mpirout.mod random.mod stringmanip.mod

mpirout.mod: mpirout.f90 genvar.mod linalg.mod random.mod
 
namelistmod.mod: namelistmod.f90 common_def.mod genvar.mod string4.mod \
                 string5.mod
 
openmpmod.mod: openmpmod.f90 lockmod.mod
 
random.mod: random.f90 common_def.mod genvar.mod linalg.mod
 
smooth_data.mod: smooth_data.f90 fourier_transform_mod.mod genvar.mod
 
sorting.mod: sorting.f90 linalg.mod random.mod
 
specialfunction.mod: specialfunction.f90 genvar.mod linalg.mod
 
splines.mod: splines.f90 derivative_noise.mod genvar.mod geometry.mod \
             mpirout.mod random.mod smooth_data.mod sorting.mod splines2.mod \
             tools_fit.mod

splines2.mod: splines2.f90

tridiag_class.mod: tridiag_class.f90 common_def.mod genvar.mod matrix.mod
 
stringmanip.mod: stringmanip.f90 
 
strings.mod: strings.f90 
 
string4.mod: string4.f90 string5.mod

string5.mod: string5.f90 stringmanip.mod
 
tools_algebra.mod: tools_algebra.f90 genvar.mod
 
tools_fit.mod: tools_fit.f90 

IEEE.o: IEEE.c 

memory_address_functions.o: memory_address_functions.cc
