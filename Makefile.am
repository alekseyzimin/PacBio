ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -Wall -I$(srcdir) -I$(srcdir)/include $(BOOST_CPPFLAGS)
AM_CXXFLAGS = $(ALL_CXXFLAGS) -std=c++0x $(JELLYFISH2_0_CFLAGS)
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
                bin/merge_coords bin/findGapsInCoverageOfPacbios bin/get_super_read_sizes

noinst_HEADERS += src_jf_aligner/zstr.hpp src_jf_aligner/strict_fstream.hpp \
		  src_jf_aligner/superread_parser.hpp			\
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
                  src_psa/prefetch_iterator_traits.hpp			\
                  src_psa/48bit_index.hpp src_psa/48bit_iterator.hpp

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

bin_get_super_read_sizes_SOURCES = src_mega_reads/get_super_read_sizes.cc

dist_bin_SCRIPTS += src_mega_reads/mega_reads_assemble_cluster.sh src_mega_reads/mega_reads_assemble_cluster2.sh \
                    src_mega_reads/make_mr_frg.pl src_mega_reads/fastqToFasta.pl \
                    src_mega_reads/reconciliate_mega_reads.maximal.nucmer.pl	\
                    src_mega_reads/join_mega_reads_trim.onepass.nomatch.pl	\
                    src_mega_reads/join_mega_reads_trim.onepass.pl		\
                    src_mega_reads/make_mate_frg.pl				\
                    src_mega_reads/determineUnjoinablePacbioSubmegas.perl \
		    src_ref/join_mega_reads_trim.onepass.ref.pl src_ref/join_mega_reads_trim.onepass.rna.pl\
		    src_mega_reads/refine.sh src_mega_reads/refine_alignments.pl  \
		    src_mega_reads/split_matches_file.pl \
		    src_ref/mega_reads_assemble_ref.sh \
		    src_reconcile/reconcile_matches.pl \
		    src_reconcile/samba.sh src_reconcile/close_scaffold_gaps.sh \
		    src_reconcile/output_reconciled_scaffolds.pl \
		    src_reconcile/merge_matches_and_tile_coords_file.pl \
		    src_mega_reads/reconcile_alignments.sh \
		    src_mega_reads/split_long_unitigs.pl \
		    src_mega_reads/recover_scaffolds.pl  src_mega_reads/trim_by_kunitigs.pl  src_mega_reads/trim_mega_reads.pl \
		    src_reconcile/fill_unaligned_gaps.pl \
		    src_reconcile/introduce_errors_fasta_file.pl \
		    src_reconcile/compute_coverage.pl  \
		    src_reconcile/evaluate_splits.pl \
		    src_reconcile/polca.sh src_reconcile/sizesToScaff.pl \
		    src_reconcile/freebayes \
		    src_reconcile/find_repeats.pl \
		    src_reconcile/insert_repeats.pl \
		    src_reconcile/fix_consensus_from_vcf.pl \
		    src_reconcile/extract_single_best_match_coords_file.pl \
		    src_reconcile/chromosome_scaffolder.sh \
		    src_reconcile/break_contigs.pl src_reconcile/merge_sequences_simple.pl \
		    src_mega_reads/get_super_read_sizes.pl src_mega_reads/translate_reduced_reads.pl \
		    src_ref/recompute_astat_superreads_CA8.sh src_mega_reads/add_pb_seq.pl \
		    src_mega_reads/deduplicate_unitigs.sh src_mega_reads/deduplicate_contigs.sh src_mega_reads/parallel_delta-filter.sh \
		    src_reconcile/polish_with_illumina_assembly.sh src_reconcile/extract_merges.pl \
		    src_reconcile/merge_contigs.pl src_reconcile/create_merged_sequences.pl \
		    src_reconcile/make_rejoin_links.pl \
		    src_mega_reads/find_contained_reads.pl src_reconcile/reconcile_consensus.pl \
		    src_mega_reads/fragScaff.sh src_mega_reads/fragScaff.pl \
		    src_mega_reads/blast_self_alignment_filter.pl src_mega_reads/fasta_make_Nbase_bed.pl \
		    src_mega_reads/create_merged_mega-reads.pl src_mega_reads/extract_merges_mega-reads.pl \
		    src_mega_reads/merge_mega-reads.pl src_mega_reads/split_reads_to_join.pl src_mega_reads/filter_delta_file_for_qrys.pl \
		    src_mega_reads/process_scaffold_gaps.pl src_mega_reads/splitScaffoldsAtNs.pl src_reconcile/splitScaffoldsAtNs.sh \
		    src_trio/classify_reads.pl src_trio/resolve_trio.sh \
		    src_mega_reads/extract_unique_mega-reads.pl src_mega_reads/mega_reads_assemble_polish.sh \
		    src_mega_reads/correct_with_k_unitigs_fast.pl src_mega_reads/correct_with_k_unitigs.sh


YAGGO_SOURCES += src_jf_aligner/jf_aligner_cmdline.hpp			\
	         src_jf_aligner/create_mega_reads_cmdline.hpp		\
	         src_jf_aligner/longest_path_overlap_graph2_cmdline.hpp	\
	         src_jf_aligner/merge_coords_cmdline.hpp		\
	         src_mega_reads/findGapsInCoverageOfPacbios_cmdline.hpp \
		 src_mega_reads/get_super_read_sizes_cmdline.hpp


# # Print the value of a variable
# print-%:
# 	@echo -n $($*)

#########
# Tests #
#########
include tests/Makefile.am

