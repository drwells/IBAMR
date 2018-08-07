## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = CCLaplace CCPoisson PETScOps PhysBdryOps SCLaplace SCPoisson VCLaplace VCViscousSolver

all:
	if test "$(CONFIGURATION_BUILD_DIR)" != ""; then \
	  $(MAKE) examples ; \
	fi ;

COMPLETE_EXAMPLES = CCLaplace CCPoisson PETScOps PhysBdryOps SCLaplace SCPoisson \
                    VCLaplace VCViscousSolver

examples: $(COMPLETE_EXAMPLES)
.PHONY: examples $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples

# gtest-PhysBdryOps is currently disabled
GTESTS = gtest-CCLaplace gtest-CCPoisson gtest-PETScOps gtest-SCLaplace \
         gtest-SCPoisson gtest-VCLaplace gtest-VCViscousSolver

# gtest-long-PhysBdryOps is currently disabled
GTESTS_LONG = gtest-long-CCLaplace gtest-long-CCPoisson gtest-long-PETScOps \
              gtest-long-SCLaplace gtest-long-SCPoisson gtest-long-VCLaplace \
              gtest-long-VCViscousSolver
gtest: $(GTESTS)

$(GTESTS):
	$(eval DIRECTORY := $(subst gtest-,,$@))
	cd $(DIRECTORY) && $(MAKE) $(AM_MAKEFLAGS) gtest

gtest-long: $(GTESTS_LONG)

$(GTESTS_LONG):
	$(eval DIRECTORY := $(subst gtest-long-,,$@))
	cd $(DIRECTORY) && $(MAKE) $(AM_MAKEFLAGS) gtest-long

gtest-all: gtest-long
.PHONY: gtest gtest-long gtest-all $(GTESTS) $(GTESTS_LONG)
