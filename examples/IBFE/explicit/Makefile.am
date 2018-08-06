## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = ex0 ex1 ex2 ex3 ex4 ex5 ex6 ex7 ex8 ex9 ex10 ex11

COMPLETE_EXAMPLES = $(SUBDIRS)
INCOMPLETE_EXAMPLES =
.PHONY: examples $(COMPLETE_EXAMPLES)
examples: $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples

GTESTS = gtest-ex0 gtest-ex1 gtest-ex2 gtest-ex3 gtest-ex4 gtest-ex5 gtest-ex6 \
         gtest-ex7 gtest-ex8 gtest-ex9 gtest-ex10 gtest-ex11

GTESTS_LONG = gtest-long-ex0 gtest-long-ex1 gtest-long-ex2 gtest-long-ex3 \
              gtest-long-ex4 gtest-long-ex5 gtest-long-ex6 gtest-long-ex7 \
              gtest-long-ex8 gtest-long-ex9 gtest-long-ex10 gtest-long-ex11

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
