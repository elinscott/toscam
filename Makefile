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
	$(MAKE) -f $(ROOTDIR)/src/scripts/makefile ARCH=$(ARCH) ROOTDIR=$(ROOTDIR))
	@echo Compilation complete

ed_solver: solver_obj
	@echo Compiling exact diagonalisation solver
	@( cd src/ed_solver/obj ; \
	$(MAKE) -f $(ROOTDIR)/src/ed_solver/makefile ARCH=$(ARCH) ROOTDIR=$(ROOTDIR))

lib: lib_obj
	@echo Compiling libraries
	@( cd src/splines/obj ; \
	$(MAKE) -f $(ROOTDIR)/src/splines/makefile ARCH=$(ARCH) ROOTDIR=$(ROOTDIR))
	@( cd src/utils/obj ; \
	$(MAKE) -f $(ROOTDIR)/src/utils/makefile ARCH=$(ARCH) ROOTDIR=$(ROOTDIR))

clean:
	@$(RM) -fr bin/*
	@( cd src/scripts ; \
	$(MAKE) -f makefile clean)
	@( cd src/ed_solver ; \
	$(MAKE) -f makefile clean)
	@( cd src/splines ; \
	$(MAKE) -f makefile clean)
	@( cd src/utils ; \
	$(MAKE) -f makefile clean)

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
	@echo '> make clean'
	@echo

src_obj:
	@( if [ ! -d src/obj ] ; \
		then mkdir src/obj ; \
	fi )

solver_obj:
	@( if [ ! -d src/ed_solver/obj ] ; \
		then mkdir src/ed_solver/obj ; \
	fi )

lib_obj:
	@( if [ ! -d src/utils/obj ] ; \
		then mkdir src/utils/obj ; \
	fi )
	@( if [ ! -d src/splines/obj ] ; \
		then mkdir src/splines/obj ; \
	fi )

