## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = CCLaplace CCPoisson PETScOps PhysBdryOps SCLaplace SCPoisson VCLaplace VCViscousSolver

all:
	if test "$(CONFIGURATION_BUILD_DIR)" != ""; then \
	  $(MAKE) examples ; \
	fi ;

## Standard make targets.
examples:
	@(cd CCLaplace   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd CCPoisson   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd PETScOps    && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd PhysBdryOps && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd SCLaplace   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd SCPoisson   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd VCLaplace   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd VCViscousSolver   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;


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
