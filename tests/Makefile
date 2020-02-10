# List of available tasks
.PHONY: tests clean help

tests:
ifeq ($(ARCH),)
	$(warning Please set ARCH.)
	@$(MAKE) help
else
	ARCH=$(ARCH) testcode/bin/testcode.py $(TC) --verbose
endif

KEEP_FILES := *.dat input_onetep_dmft.txt mask_projections
KEEP_TESTS := $(patsubst %,! -name %,${KEEP_FILES})
CLEAN_FILES := test.out.* test.err.*
CLEAN_TESTS := $(patsubst %,-o -name %,${CLEAN_FILES})

clean:
	find test_?? -type f ${KEEP_TESTS} ${CLEAN_TESTS} | xargs -i rm -v {}

help:
	@echo ''
	@echo 'TOSCAM QC tests'
	@echo ''
	@echo 'To run type'
	@echo '> make tests ARCH=<arch>'
	@echo
	@echo 'Other options'
	@echo 'To remove all test outputs in the test subdirectories:'
	@echo '> make clean'
	@echo

