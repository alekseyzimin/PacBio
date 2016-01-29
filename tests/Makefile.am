TESTS = $(UNITTESTS)
UNITTESTS = %D%/test_super_read_name %D%/test_union_find		\
%D%/test_pb_aligner %D%/test_superread_parser				\
%D%/test_align_coordinates %D%/test_kmers_info				\
%D%/test_coords_parsing %D%/test_tiling %D%/test_compute_coords_info
# check_PROGRAMS = %D%/speed_super_read_name $(UNITTESTS)
check_PROGRAMS = $(UNITTESTS)
check_LIBRARIES = %D%/libgtest.a %D%/libgtest_main.a

AM_CPPFLAGS += -I$(srcdir)/%D%
TEST_LDADD = -L%D% -lgtest_main
%C%_libgtest_a_SOURCES = %D%/gtest/src/gtest-all.cc %D%/misc.cc %D%/gtest/gtest.h %D%/misc.hpp
%C%_libgtest_main_a_SOURCES = %D%/gtest/src/gtest-all.cc %D%/gtest/src/gtest_main.cc %D%/misc.cc

%C%_speed_super_read_name_SOURCES = %D%/speed_super_read_name.cc src_jf_aligner/super_read_name.cc

%C%_test_super_read_name_SOURCES = %D%/test_super_read_name.cc src_jf_aligner/super_read_name.cc
%C%_test_super_read_name_LDADD = $(TEST_LDADD)

%C%_test_union_find_SOURCES = %D%/test_union_find.cc src_jf_aligner/union_find.cc
%C%_test_union_find_LDADD = $(TEST_LDADD)

%C%_test_pb_aligner_SOURCES = %D%/test_pb_aligner.cc			\
                              src_jf_aligner/pb_aligner.cc		\
                              src_jf_aligner/coarse_aligner.cc		\
                              src_jf_aligner/superread_parser.cc	\
                              src_jf_aligner/super_read_name.cc		\
                              lib/global_timer.cc
%C%_test_pb_aligner_LDADD = $(TEST_LDADD)

%C%_test_superread_parser_SOURCES = %D%/test_superread_parser.cc	\
                                    src_jf_aligner/superread_parser.cc	\
                                    src_jf_aligner/super_read_name.cc	\
                                    lib/global_timer.cc
%C%_test_superread_parser_LDADD = $(TEST_LDADD)

%C%_test_align_coordinates_SOURCES = %D%/test_align_coordinates.cc	\
                                     src_jf_aligner/pb_aligner.cc	\
                                     src_jf_aligner/superread_parser.cc	\
                                     src_jf_aligner/super_read_name.cc
%C%_test_align_coordinates_LDADD = $(TEST_LDADD)

%C%_test_kmers_info_SOURCES = %D%/test_kmers_info.cc		\
                              src_jf_aligner/pb_aligner.cc	\
                              src_jf_aligner/super_read_name.cc
%C%_test_kmers_info_LDADD = $(TEST_LDADD)

%C%_test_coords_parsing_SOURCES = %D%/test_coords_parsing.cc		\
                                  src_jf_aligner/coords_parsing.cc	\
                                  src_jf_aligner/super_read_name.cc
%C%_test_coords_parsing_LDADD = $(TEST_LDADD)

%C%_test_tiling_SOURCES = %D%/test_tiling.cc			\
                          src_jf_aligner/overlap_graph.cc	\
                          src_jf_aligner/super_read_name.cc	\
                          src_jf_aligner/union_find.cc
%C%_test_tiling_LDADD = $(TEST_LDADD)

%C%_test_compute_coords_info_SOURCES = %D%/test_compute_coords_info.cc		\
                                       src_jf_aligner/pb_aligner.cc		\
                                       src_jf_aligner/super_read_name.cc
%C%_test_compute_coords_info_LDADD = $(TEST_LDADD)
