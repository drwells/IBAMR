## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = IB IBFE IMP adv_diff advect navier_stokes ConstraintIB CIB CIBFE level_set multiphase_flow vc_navier_stokes

all:
	if test "$(CONFIGURATION_BUILD_DIR)" != ""; then \
	  $(MAKE) examples ; \
	fi ;

COMPLETE_EXAMPLES = IB IBFE adv_diff advect navier_stokes ConstraintIB CIB \
                    level_set multiphase_flow vc_navier_stokes
# Not all examples can be compiled at the present time
INCOMPLETE_EXAMPLES = IMP CIBFE
examples: $(COMPLETE_EXAMPLES)
.PHONY: examples $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples

# gtest-IMP and gtest-CIBFE are currently disabled
GTESTS = gtest-IB gtest-IBFE gtest-adv_diff gtest-advect gtest-navier_stokes \
         gtest-ConstraintIB gtest-CIB gtest-level_set gtest-multiphase_flow \
         gtest-vc_navier_stokes

# gtest-long-IMP and gtest-long-CIBFE are currently disabled
GTESTS_LONG = gtest-long-IB gtest-long-IBFE gtest-long-adv_diff \
              gtest-long-advect gtest-long-navier_stokes gtest-long-ConstraintIB \
              gtest-long-CIB gtest-long-level_set gtest-long-multiphase_flow \
              gtest-long-vc_navier_stokes

if GTEST_ENABLED
gtest: $(GTESTS)

$(GTESTS):
	$(eval DIRECTORY := $(subst gtest-,,$@))
	cd $(DIRECTORY) && $(MAKE) $(AM_MAKEFLAGS) gtest

gtest-long: $(GTESTS_LONG)

$(GTESTS_LONG):
	$(eval DIRECTORY := $(subst gtest-long-,,$@))
	cd $(DIRECTORY) && $(MAKE) $(AM_MAKEFLAGS) gtest-long

gtest-all: gtest-long
else
gtest:
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	@echo " To run tests, Google Test Framework must be enabled.                "
	@echo " Configure IBAMR with additional options:                            "
	@echo "                      --enable-gtest --with-gtest=path               "
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
gtest-long: gtest
gtest-all: gtest
endif
.PHONY: gtest gtest-long gtest-all $(GTESTS) $(GTESTS_LONG)
