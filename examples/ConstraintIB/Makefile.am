## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = eel2d eel3d

COMPLETE_EXAMPLES = eel2d eel3d falling_sphere flow_past_cylinder \
                    flow_past_cylinder_HF knifefish impulsively_started_cylinder \
                    moving_plate oscillating_rigid_cylinder stokes_first_problem
INCOMPLETE_EXAMPLES =
.PHONY: examples $(COMPLETE_EXAMPLES)
examples: $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples

GTESTS = gtest-eel2d gtest-eel3d gtest-falling_sphere gtest-flow_past_cylinder \
         gtest-flow_past_cylinder_HF gtest-knifefish \
         gtest-impulsively_started_cylinder gtest-moving_plate \
         gtest-oscillating_rigid_cylinder gtest-stokes_first_problem
GTESTS_LONG = gtest-long-eel2d gtest-long-eel3d gtest-long-falling_sphere \
              gtest-long-flow_past_cylinder gtest-long-flow_past_cylinder_HF \
              gtest-long-knifefish gtest-long-impulsively_started_cylinder \
              gtest-long-moving_plate gtest-long-oscillating_rigid_cylinder \
              gtest-long-stokes_first_problem

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
