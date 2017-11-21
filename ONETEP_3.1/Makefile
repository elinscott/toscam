#
# Top-level makefile for ONETEP
#
# This version for multiple architectures by Peter Haynes
#
# System-specific changes should be made to the appropriate
#  conf.arch files
#

ROOTDIR = $(PWD)

ARCH = `$(ROOTDIR)/bin/arch`

TARGET ?= onetep

default: help

onetep: objdir libdir
	@( cd obj/$(ARCH) ; \
	$(MAKE) -f $(ROOTDIR)/src/Makefile $(TARGET) VPATH=$(ROOTDIR)/src:$(VPATH) \
		ARCH=$(ARCH) ROOTDIR=$(ROOTDIR) )

debug:
	$(MAKE) onetep TARGET=debug

timing:
	$(MAKE) onetep TARGET=timing

profile:
	$(MAKE) onetep TARGET=profile

warnings:
	$(MAKE) onetep TARGET=warnings

utils: utildir
	@( cd utils ; \
	$(MAKE) -f Makefile ARCH=$(ARCH) ROOTDIR=$(ROOTDIR) EXTFLAGS=$(EXTFLAGS) )

clean:
	@$(RM) -fr obj/$(ARCH)

cleanall: cleanutils
	@$(RM) -fr obj/* bin/onetep.* *~ src/*~

cleanutils:
	@(cd utils ; \
	$(MAKE) -f Makefile clean ARCH=$(ARCH) ROOTDIR=$(ROOTDIR) )

dist:
	@( date > Archived.txt ; \
	tar -cf onetep.tar Archived.txt Makefile bin config include lib \
        obj src utils INSTALL )

help:
	@echo
	@echo '        ####### #     # ####### ####### ####### ######   '
	@echo '        #     # ##    # #          #    #       #     #  '
	@echo '        #     # # #   # #          #    #       #     #  '
	@echo '        #     # #  #  # #####      #    #####   ######   '
	@echo '        #     # #   # # #          #    #       #        '
	@echo '        #     # #    ## #          #    #       #        '
	@echo '        ####### #     # #######    #    ####### #        '
	@echo
	@echo 'To install type:'
	@echo '> gmake onetep'
	@echo
	@echo 'Other options:'
	@echo 'To compile a version for debugging:'
	@echo '> gmake clean ; gmake debug'
	@echo 'To compile a version for profiling:'
	@echo '> gmake clean ; gmake profile'
	@echo 'To remove object files for this architecture:'
	@echo '> gmake clean'
	@echo 'To remove object files for all architectures:'
	@echo '> gmake cleanall'
	@echo 'To create an archive for distribution:'
	@echo '> gmake dist'
	@echo 'To make utilities associated with ONETEP:'
	@echo '> gmake utils'
	@echo

checkdeps:
	utils/check_dependencies

objdir:
	@( if [ ! -d obj ] ; \
	        then mkdir obj ; \
	fi )
	@( if [ ! -d obj/$(ARCH) ] ; \
		then mkdir obj/$(ARCH) ; \
	fi )

libdir:
	@( if [ ! -d lib ] ; \
	        then mkdir lib ; \
	fi )
	@( if [ ! -d lib/$(ARCH) ] ; \
	        then mkdir lib/$(ARCH) ; \
	fi )


utildir:
	@( if [ ! -d utils ] ; \
		then echo 'utils directory not found!' ; exit 1 ; \
	fi )
