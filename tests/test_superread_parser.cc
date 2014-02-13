#include <unistd.h>
#include <cstdlib>
#include <string>

#include <gtest/gtest.h>
#include <src_jf_aligner/superread_parser.hpp>

namespace {
struct remove_file {
  const char* path;
  remove_file(const char* p) : path(p) { }
  ~remove_file() { unlink(path); }
};

void create_sequence(const char* path) {
  static const char bases[4] = { 'A', 'C', 'G','T' };
  std::ofstream f(path);
  if(!f.good())
    throw std::runtime_error("Failed to open sequence file");
  f << ">superread\n";
  for(int i = 0; i < 100; ++i)
    f << bases[random() % 4];
  f << "\n";
}

TEST(SuperReadParser, Parse) {
  remove_file file("/tmp/superread.fa");
  create_sequence(file.path);

  mer_dna::k(17);

  mer_pos_hash_type hash(1024);
  name_lists names;
  superread_parse(1, hash, names, file.path);

  // Check every k-mer in the sequence
  std::ifstream f(file.path);
  if(!f.good())
    throw std::runtime_error("Failed to open sequence file");
  std::string line;
  std::getline(f, line); // Skip header
  std::getline(f, line);

  for(size_t i = 0; i < line.size() - mer_dna::k() + 1; ++i) {
    mer_dna m(line.substr(i, mer_dna::k()));
    mer_dna rm(m.get_reverse_complement());
    bool is_canonical = m < rm;

    const mer_pos_hash_type::mapped_type* list = hash.find_pos(is_canonical ? m : rm);
    ASSERT_TRUE(list != 0);
    EXPECT_EQ(1, std::distance(list->cbegin(), list->cend()));
    const mer_pos_hash_type::position_type& pos = *list->cbegin();
    EXPECT_EQ("superread", *pos.frag);
    EXPECT_EQ((int)i * (is_canonical ? 1 : -1), pos.offset);
  }
}
}
