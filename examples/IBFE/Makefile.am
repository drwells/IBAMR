## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = explicit

COMPLETE_EXAMPLES = $(SUBDIRS)
INCOMPLETE_EXAMPLES =
.PHONY: examples $(COMPLETE_EXAMPLES)
examples: $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples

gtest:
	cd explicit && $(MAKE) $(AM_MAKEFLAGS) gtest

gtest-long:
	cd explicit && $(MAKE) $(AM_MAKEFLAGS) gtest-long

gtest-all: gtest-long
.PHONY: gtest gtest-long gtest-all
