#include <cstdlib>

#include <gtest/gtest.h>
#include <src_jf_aligner/mer_pos_hash.hpp>

namespace {
using jellyfish::mer_dna;
typedef mer_pos_hash<mer_dna> mer_pos_hash_type;
typedef mer_pos_hash_type::position_type position_type;
typedef mer_pos_hash_type::mapped_type list_type;

TEST(MerPosHash, Insert) {
  mer_dna::k(17);
  static const int nb_mers = 1000;
  std::vector<mer_dna> mers(nb_mers);
  for(int i = 0; i < nb_mers; ++i)
    mers[i].randomize();

  static const int nb_frags = 37;
  std::vector<std::string> frags(nb_frags);
  for(int i = 0; i < nb_frags; ++i)
    frags[i] += "frag_" + std::to_string(i);

  mer_pos_hash_type hash(2 * nb_mers);

  static const int nb_inserts = 10000;
  # pragma omp parallel for
  for(int i = 0; i < nb_inserts; ++i)
    hash.push_front(mers[i % nb_mers], frags[i % nb_frags], i * (2 * (i % 2) - 1));

  for(int i = 0; i < nb_mers; ++i) {
    const list_type& list = hash[mers[i]];
    const list_type* found_list = hash.find_pos(mers[i]);
    EXPECT_EQ(&list, found_list);
    for(auto it = list.cbegin(); it != list.cend(); ++it) {
      EXPECT_TRUE(abs(it->offset) >= 0 && it->offset < nb_inserts);
      EXPECT_EQ(i, abs(it->offset) % nb_mers);
      EXPECT_EQ(i % 2,  it->offset > 0);
      EXPECT_EQ(frags[abs(it->offset) % nb_frags], *it->frag);
    }
  }
}
}
