bin_PROGRAMS = phyltr-dp phyltr-fpt phyltr-gen-gtree phyltr-gen-stree \
               phyltr-event-combinations
AM_LDFLAGS = $(BOOST_LDFLAGS) $(NHPARSER_LDFLAGS)
AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(NHPARSER_CPPFLAGS)
LDADD = $(BOOST_PROGRAM_OPTIONS_LIB) $(NHPARSER_LIB)

phyltr_dp_SOURCES = phyltr-dp.cc common.cc common.hh common-impl.hh \
                    utils/binary-tree.hh utils/binary-tree-impl.hh \
		    utils/assert.hh

phyltr_fpt_SOURCES = phyltr-fpt.cc common.cc common.hh common-impl.hh \
                     utils/binary-tree.hh  utils/binary-tree-impl.hh \
		     utils/assert.hh

phyltr_gen_gtree_SOURCES = phyltr-gen-gtree.cc common.cc common.hh \
                           common-impl.hh

phyltr_gen_stree_SOURCES = phyltr-gen-stree.cc common.cc common.hh \
                           common-impl.hh

phyltr_event_combinations_SOURCES = phyltr-event-combinations.cc common.cc \
                             common.hh common-impl.hh utils/binary-tree.hh \
			     utils/binary-tree-impl.hh utils/assert.hh
