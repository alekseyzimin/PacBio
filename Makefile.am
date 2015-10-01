ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -Wall -I$(srcdir) -I$(srcdir)/include
AM_CXXFLAGS = $(ALL_CXXFLAGS) -g -O3 $(JELLYFISH2_0_CFLAGS)
AM_LDFLAGS = $(JELLYFISH2_0_LIBS)

noinst_HEADERS = $(YAGGO_SOURCES)
bin_PROGRAMS =
dist_bin_SCRIPTS =
data_DATA =
BUILT_SOURCES = $(YAGGO_SOURCES)
CLEANFILES =
DISTCLEANFILES = $(BUILT_SOURCES)

# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO   " $@;
.yaggo.hpp:
#	$(V_YAGGO)$(YAGGO) --license $(srcdir)/header-license -o $@ $<
	$(V_YAGGO)$(YAGGO) -o $@ $<

YAGGO_SOURCES = # Append all file to be built by yaggo

##################
# Build binaries #
##################
bin_PROGRAMS += bin/jf_aligner bin/create_mega_reads bin/longest_path	\
                bin/merge_coords bin/findGapsInCoverageOfPacbios

noinst_HEADERS += src_jf_aligner/superread_parser.hpp			\
                  src_jf_aligner/super_read_name.hpp			\
                  src_jf_aligner/fine_aligner.hpp			\
                  src_jf_aligner/coarse_aligner.hpp			\
                  src_jf_aligner/overlap_graph.hpp			\
                  src_jf_aligner/misc.hpp				\
                  src_jf_aligner/coords_parsing.hpp			\
                  src_jf_aligner/union_find.hpp				\
                  src_jf_aligner/pb_aligner.hpp				\
                  src_jf_aligner/output_file.hpp			\
                  src_jf_aligner/jf_aligner.hpp				\
                  src_jf_aligner/frag_info.hpp				\
                  src_jf_aligner/lf_forward_size_list.hpp		\
                  src_jf_aligner/lf_forward_list.hpp			\
                  src_jf_aligner/least_square_2d.hpp			\
                  src_jf_aligner/bounded_counter.hpp src_lis/lis.hpp	\
                  src_lis/lis_align.hpp src_psa/barrier.hpp		\
                  src_psa/global_timer.hpp src_psa/psa.hpp		\
                  src_psa/compact_dna.hpp src_psa/mer_sa.hpp		\
                  src_psa/said_traits.hpp src_psa/compact_index.hpp	\
                  src_psa/mer_sa_imp.hpp src_psa/slice.hpp		\
                  src_psa/compact_iterator.hpp				\
                  src_psa/parallel_iterator_traits.hpp			\
                  src_psa/const_iterator_traits.hpp			\
                  src_psa/prefetch_iterator_traits.hpp

noinst_HEADERS += include/charb.hpp include/exp_buffer.hpp		\
                  include/multiplexer.hpp include/reallocators.hpp	\
                  include/debug.hpp include/misc.hpp			\
                  include/jflib/multiplexed_parser.hpp			\
                  include/jflib/multiplexed_io.hpp			\
                  include/jflib/pool.hpp				\
                  include/boost/asio/coroutine.hpp			\
                  include/boost/asio/yield.hpp

bin_jf_aligner_SOURCES = src_jf_aligner/jf_aligner.cc			\
                         src_jf_aligner/pb_aligner.cc			\
                         src_jf_aligner/superread_parser.cc		\
                         src_jf_aligner/super_read_name.cc		\
                         src_jf_aligner/coarse_aligner.cc		\
                         src_jf_aligner/fine_aligner.cc			\
                         src_jf_aligner/misc.cc lib/global_timer.cc

bin_create_mega_reads_SOURCES = src_jf_aligner/create_mega_reads.cc	\
                                src_jf_aligner/pb_aligner.cc		\
                                src_jf_aligner/superread_parser.cc	\
                                src_jf_aligner/union_find.cc		\
                                src_jf_aligner/overlap_graph.cc		\
                                src_jf_aligner/super_read_name.cc	\
                                src_jf_aligner/coarse_aligner.cc	\
                                src_jf_aligner/fine_aligner.cc		\
                                src_jf_aligner/misc.cc			\
                                lib/global_timer.cc
bin_longest_path_SOURCES =							\
                           src_jf_aligner/longest_path_overlap_graph2.cc	\
                           src_jf_aligner/pb_aligner.cc				\
                           src_jf_aligner/superread_parser.cc			\
                           src_jf_aligner/union_find.cc				\
                           src_jf_aligner/overlap_graph.cc			\
                           src_jf_aligner/super_read_name.cc			\
                           src_jf_aligner/coords_parsing.cc			\
                           src_jf_aligner/misc.cc
bin_merge_coords_SOURCES = src_jf_aligner/merge_coords.cc

bin_findGapsInCoverageOfPacbios_SOURCES = src_mega_reads/findGapsInCoverageOfPacbios.cc

dist_bin_SCRIPTS += src_mega_reads/build_mega_read_graph.sh			\
                    src_mega_reads/mega_reads_assemble.sh			\
                    src_mega_reads/mega_reads_assemble_nomatch.sh		\
                    src_mega_reads/run_big_nucmer_job_parallel.sh		\
                    src_mega_reads/run_big_nucmer_job_parallel_blasr_out.sh     \
                    src_mega_reads/build_mega_reads.pl				\
                    src_mega_reads/make_mr_frg.pl				\
                    src_mega_reads/extract_best_match_coords.pl			\
                    src_mega_reads/merge_matches_coords_file.pl			\
                    src_mega_reads/filter_matches.pl				\
                    src_mega_reads/reconciliate_mega_reads.maximal.nucmer.pl	\
                    src_mega_reads/join_mega_reads_trim.onepass.nomatch.pl	\
                    src_mega_reads/split_contig_file.pl				\
                    src_mega_reads/join_mega_reads_trim.onepass.pl		\
                    src_mega_reads/split_long_unitigs.pl			\
                    src_mega_reads/make_mate_frg.pl				\
                    src_mega_reads/determineUnjoinablePacbioSubmegas.perl refine_alignments.pl



YAGGO_SOURCES += src_jf_aligner/jf_aligner_cmdline.hpp			\
	         src_jf_aligner/create_mega_reads_cmdline.hpp		\
	         src_jf_aligner/longest_path_overlap_graph2_cmdline.hpp	\
	         src_jf_aligner/merge_coords_cmdline.hpp		\
	         src_mega_reads/findGapsInCoverageOfPacbios_cmdline.hpp


# # Print the value of a variable
# print-%:
# 	@echo -n $($*)
