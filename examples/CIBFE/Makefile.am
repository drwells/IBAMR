## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = ex0 ex1

COMPLETE_EXAMPLES = $(SUBDIRS)
INCOMPLETE_EXAMPLES =
.PHONY: examples $(COMPLETE_EXAMPLES)
examples: $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples

GTESTS = gtest-ex0 gtest-ex1
GTESTS_LONG = gtest-long-ex0 gtest-long-ex1

# tests are disabled for CIBFE
gtest: $(GTESTS)

# cd $(shell echo $@ | sed "s/^gtest-//") && $(MAKE) $(AM_MAKEFLAGS) gtest
$(GTESTS):

gtest-long: $(GTESTS_LONG)

# cd $(shell echo $@ | sed "s/^gtest-long-//") && $(MAKE) $(AM_MAKEFLAGS) gtest-long
$(GTESTS_LONG):

gtest-all: gtest-long
.PHONY: gtest gtest-long gtest-all $(GTESTS) $(GTESTS_LONG)
