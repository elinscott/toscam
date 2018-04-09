#
# Top-level makefile for TOSCAM
#
# Do not edit. System-specific changes should be made to the appropriate
# conf.arch files
#

ROOTDIR = $(PWD)

ARCH = `$(ROOTDIR)/bin/arch`

.PHONY: default toscam ed_solver utils clean cleanall help

default: help

toscam: src_obj lib ed_solver
	@echo Compiling toscam scripts
	@( cd bin ; \
	$(MAKE) -f $(ROOTDIR)/src/makefile ARCH=$(ARCH) ROOTDIR=$(ROOTDIR))
	@echo Compilation complete

ed_solver: solver_obj
	@echo Compiling exact diagonalisation solver
	@( cd ed_solver/obj ; \
	$(MAKE) -f $(ROOTDIR)/ed_solver/makefile ARCH=$(ARCH) ROOTDIR=$(ROOTDIR))

lib: lib_obj
	@echo Compiling libraries
	@( cd lib/slatec/static ; \
	$(MAKE) -f makefile ARCH=$(ARCH) ROOTDIR=$(ROOTDIR))
	@( cd lib/splines/obj ; \
	$(MAKE) -f $(ROOTDIR)/lib/splines/makefile ARCH=$(ARCH) ROOTDIR=$(ROOTDIR))
	@( cd lib/utils/obj ; \
	$(MAKE) -f $(ROOTDIR)/lib/utils/makefile ARCH=$(ARCH) ROOTDIR=$(ROOTDIR))

clean:
	@$(RM) -fr bin/*
	@( cd src ; \
	$(MAKE) -f makefile clean)
	@( cd ed_solver ; \
	$(MAKE) -f makefile clean)

cleanlib:
	@( cd lib/slatec/static ; \
	$(MAKE) -f makefile clean)
	@( cd lib/splines ; \
	$(MAKE) -f makefile clean)
	@( cd lib/utils ; \
	$(MAKE) -f makefile clean)

cleanall: clean cleanlib

help:
	@echo
	@echo '         _____ ___  ___  ___   _   __  __          '
	@echo '        |_   _/ _ \/ __|/ __| /_\ |  \/  |         '
	@echo '          | || (_) \__ \ (__ / _ \| |\/| |         '
	@echo '          |_| \___/|___/\___/_/ \_\_|  |_|         '
	@echo
	@echo
	@echo 'To install type'
	@echo '> make toscam'
	@echo
	@echo 'Other options'
	@echo 'To compile using a specific architecture:'
	@echo '> make toscam ARCH=<arch>'
	@echo 'To make the exact diagonalisation solver:'
	@echo '> make ed_solver'
	@echo 'To make libraries associated with TOSCAM:'
	@echo '> make lib'
	@echo 'To remove all compiled files:'
	@echo '> make cleanall'
	@echo

src_obj:
	@( if [ ! -d src/obj ] ; \
		then mkdir src/obj ; \
	fi )

solver_obj:
	@( if [ ! -d ed_solver/obj ] ; \
		then mkdir ed_solver/obj ; \
	fi )

lib_obj:
	@( if [ ! -d lib/utils/obj ] ; \
		then mkdir lib/utils/obj ; \
	fi )
	@( if [ ! -d lib/splines/obj ] ; \
		then mkdir lib/splines/obj ; \
	fi )

