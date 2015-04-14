ACLOCAL_AMFLAGS = -I m4

AM_LDFLAGS = -lpthread
AM_CPPFLAGS = -Wall -I$(top_srcdir) -I$(top_srcdir)/include
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
	$(V_YAGGO)$(YAGGO) --license $(srcdir)/header-license -o $@ $<

YAGGO_SOURCES = # Append all file to be built by yaggo

##################
# Build binaries #
##################
bin_PROGRAMS += bin/jf_aligner bin/create_mega_reads bin/longest_path bin/merge_coords

AL_SRC = src_jf_aligner
bin_jf_aligner_SOURCES = src_jf_aligner/jf_aligner.cc src_jf_aligner/pb_aligner.cc	\
                      src_jf_aligner/superread_parser.cc			\
                      src_jf_aligner/super_read_name.cc			\
                      src_jf_aligner/coarse_aligner.cc			\
                      src_jf_aligner/fine_aligner.cc src_jf_aligner/misc.cc
bin_create_mega_reads_SOURCES = src_jf_aligner/create_mega_reads.cc	\
                             src_jf_aligner/pb_aligner.cc		\
                             src_jf_aligner/superread_parser.cc	\
                             src_jf_aligner/union_find.cc		\
                             src_jf_aligner/overlap_graph.cc		\
                             src_jf_aligner/super_read_name.cc	\
                             src_jf_aligner/coarse_aligner.cc		\
                             src_jf_aligner/fine_aligner.cc		\
                             src_jf_aligner/misc.cc
bin_longest_path_SOURCES = src_jf_aligner/longest_path_overlap_graph2.cc		\
                        src_jf_aligner/pb_aligner.cc				\
                        src_jf_aligner/superread_parser.cc			\
                        src_jf_aligner/union_find.cc				\
                        src_jf_aligner/overlap_graph.cc			\
                        src_jf_aligner/super_read_name.cc			\
                        src_jf_aligner/coords_parsing.cc src_jf_aligner/misc.cc
bin_merge_coords_SOURCES = src_jf_aligner/merge_coords.cc

YAGGO_SOURCES += src_jf_aligner/jf_aligner_cmdline.hpp			\
	         src_jf_aligner/create_mega_reads_cmdline.hpp		\
	         src_jf_aligner/longest_path_overlap_graph2_cmdline.hpp	\
	         src_jf_aligner/merge_coords_cmdline.hpp

# ###############
# # Build tests #
# ###############
# bin_generate_sequence_SOURCES = jellyfish/src_jf_aligner/generate_sequence.cc		\
#                                 jellyfish/src_jf_aligner/mersenne.cpp			\
#                                 jellyfish/src_jf_aligner/backtrace.cc			\
#                                 jellyfish/src_jf_aligner/dbg.cc
# YAGGO_SOURCES += jellyfish/src_jf_aligner/generate_sequence_cmdline.hpp

# #########
# # Tests #
# #########
# TEST_EXTENSIONS = .sh
# SH_LOG_COMPILER = /bin/sh
# AM_SH_LOG_FLAGS =

# TESTS = tests/generate_sequence.sh tests/parallel_hashing.sh	\
#         tests/merge.sh tests/bloom_filter.sh tests/big.sh	\
#         tests/subset_hashing.sh tests/multi_file.sh		\
#         tests/bloom_counter.sh tests/large_key.sh

# EXTRA_DIST += $(TESTS)
# clean-local: clean-local-check
# .PHONY: clean-local-check
# clean-local-check:
# 	-cd tests; rm -f *

# tests/parallel_hashing.log: tests/generate_sequence.log
# tests/subset_hashing.log: tests/generate_sequence.log
# tests/bloom_filter.log: tests/generate_sequence.log
# tests/bloom_counter.log: tests/generate_sequence.log
# tests/multi_file.log: tests/generate_sequence.log
# tests/merge.log: tests/generate_sequence.log
# tests/min_qual.log: tests/generate_fastq_sequence.log
# tests/large_key.log: tests/generate_sequence.log
# tests/quality_filter.log: tests/generate_sequence.log

# # SWIG tests
# TESTS += tests/swig_python.sh tests/swig_ruby.sh tests/swig_perl.sh
# tests/swig_python.log: tests/generate_sequence.log
# tests/swig_ruby.log: tests/generate_sequence.log
# tests/swig_perl.log: tests/generate_sequence.log
# EXTRA_DIST += swig/python/test_mer_file.py swig/python/test_hash_counter.py
# EXTRA_DIST += swig/ruby/test_mer_file.rb swig/ruby/test_hash_counter.rb
# EXTRA_DIST += swig/perl5/t/test_mer_file.t swig/perl5/t/test_hash_counter.t


# ##############
# # Unit tests #
# ##############
# TESTS += unit_tests/unit_tests.sh
# check_PROGRAMS += bin/test_all

# bin_test_all_SOURCES = unit_tests/src_jf_aligner/test_main.cc				\
# 	               unit_tests/src_jf_aligner/test_misc.cc				\
# 	               unit_tests/src_jf_aligner/test_offsets_key_value.cc		\
# 	               unit_tests/src_jf_aligner/test_simple_circular_buffer.cc	\
# 	               unit_tests/src_jf_aligner/test_rectangular_binary_matrix.cc	\
# 	               unit_tests/src_jf_aligner/test_mer_dna.cc			\
# 	               unit_tests/src_jf_aligner/test_large_hash_array.cc		\
# 	               unit_tests/src_jf_aligner/test_mer_overlap_sequence_parser.cc	\
# 	               unit_tests/src_jf_aligner/test_file_header.cc			\
# 	               unit_tests/src_jf_aligner/test_mer_iterator.cc			\
# 	               unit_tests/src_jf_aligner/test_hash_counter.cc			\
# 	               unit_tests/src_jf_aligner/test_mer_heap.cc			\
# 	               unit_tests/src_jf_aligner/test_stream_iterator.cc		\
# 	               unit_tests/src_jf_aligner/test_token_ring.cc			\
# 	               unit_tests/src_jf_aligner/test_text_dumper.cc			\
# 	               unit_tests/src_jf_aligner/test_dumpers.cc			\
# 	               unit_tests/src_jf_aligner/test_mapped_file.cc			\
# 	               unit_tests/src_jf_aligner/test_int128.cc			\
# 	               unit_tests/src_jf_aligner/test_mer_dna_bloom_counter.cc		\
# 	               unit_tests/src_jf_aligner/test_whole_sequence_parser.cc		\
# 	               unit_tests/src_jf_aligner/test_allocators_mmap.cc		\
# 	               unit_tests/src_jf_aligner/test_cooperative_pool2.cc		\
# 	               unit_tests/src_jf_aligner/test_generator_manager.cc		\
# 	               unit_tests/src_jf_aligner/test_atomic_bits_array.cc		\
# 	               unit_tests/src_jf_aligner/test_stdio_filebuf.cc
# bin_test_all_SOURCES += jellyfish/src_jf_aligner/backtrace.cc

# bin_test_all_CPPFLAGS = -Dprotected=public -Dprivate=public -DJSON_IS_AMALGAMATION=1
# bin_test_all_CXXFLAGS = $(AM_CXXFLAGS) -I$(top_srcdir)/unit_tests/gtest/include -I$(top_srcdir)/unit_tests -I$(top_srcdir)/include
# bin_test_all_LDADD = libgtest.la $(LDADD)
# YAGGO_SOURCES += unit_tests/src_jf_aligner/test_main_cmdline.hpp
# noinst_HEADERS += unit_tests/src_jf_aligner/test_main.hpp

# include gtest.mk
# -include $(top_srcdir)/development.mk

# # Print the value of a variable
# print-%:
# 	@echo -n $($*)
