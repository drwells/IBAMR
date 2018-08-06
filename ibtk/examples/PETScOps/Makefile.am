## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = ProlongationMat

all:
	if test "$(CONFIGURATION_BUILD_DIR)" != ""; then \
	  $(MAKE) examples ; \
	fi ;

gtest:
	cd ProlongationMat && $(MAKE) $(AM_MAKEFLAGS) gtest

gtest-long:
	cd ProlongationMat && $(MAKE) $(AM_MAKEFLAGS) gtest-long
.PHONY: gtest gtest-long

## Standard make targets.
examples:
	@(cd ProlongationMat && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
