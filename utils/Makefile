#
# Second-level makefile for utilities for ONES package
#
# This version for multiple architectures by Peter Haynes
#
# System-specific changes should be made to the appropriate
#  conf.arch files
#

EXECDIR = $(ROOTDIR)/bin
INCPATH = -I$(ROOTDIR)/include
LIBPATH = -L$(ROOTDIR)/lib/$(ARCH)

include $(ROOTDIR)/config/conf.$(ARCH)

UTILS = gen_coul_pot.x usp2upf.x

default: $(UTILS:.x=)

clean:
	@( $(RM) -f *~ ; \
	cd $(EXECDIR) ; \
	$(RM) -f $(UTILS:.x=.*) )

gen_coul_pot: generate_coulomb_potentials.f90
	$(F90) -o $(EXECDIR)/gen_coul_pot.$(ARCH) \
		$(FFLAGS) $(OPTFLAGS) $(INCPATH) \
		generate_coulomb_potentials.f90

usp2upf: usp2upf.F90
	$(F90) -o $(EXECDIR)/usp2upf.$(ARCH) \
        $(FFLAGS) $(OPTFLAGS) $(INCPATH) $(LIBS) \
        usp2upf.F90
