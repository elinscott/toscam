# List of available tasks
.PHONY: install tests clean help

install:
	@# install testcode
	git submodule init
	git submodule update

tests:
ifeq ($(wildcard ./testcode/bin/testcode.py),)
	@$(MAKE) install
endif
ifeq ($(ARCH),)
	$(warning Please set ARCH.)
	@$(MAKE) help
else
	ARCH=$(ARCH) testcode/bin/testcode.py $(TC) --verbose
endif

KEEP_FILES := bench* *.dat input_onetep_dmft.txt mask_projections
KEEP_TESTS := $(patsubst %,! -name %,${KEEP_FILES})
CLEAN_FILES := test.out.* test.err.*
CLEAN_TESTS := $(patsubst %,-o -name %,${CLEAN_FILES})

clean:
	find test_?? -type f ${KEEP_TESTS} ${CLEAN_TESTS} | xargs -i rm -v {}

help:
	@echo ''
	@echo 'TOSCAM QC tests'
	@echo ''
	@echo 'To run:'
	@echo '> make tests ARCH=<arch>'
	@echo 'To remove all test outputs in the test subdirectories:'
	@echo '> make clean'
	@echo ''
	@echo 'For information on testcode see http://testcode.readthedocs.org'
	@echo ''

