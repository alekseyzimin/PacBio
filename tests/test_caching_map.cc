#include <set>

#include <gtest/gtest.h>
#include <src_jf_aligner/caching_map.hpp>

namespace {
struct del_count {
  static int dels;
  int i;
  ~del_count() { ++dels; }
};
int del_count::dels = 0;

TEST(CachingMap, Insert) {
  caching_map<int, del_count> map;
  std::set<del_count*>        value_set;

  {
    auto res = map.insert(std::make_pair(0, del_count()));
    EXPECT_TRUE(res.second);
    EXPECT_EQ(0, (*res.first).first);
    (*res.first).second.i = 0;
    auto val_res = value_set.insert(&(*res.first).second);
    EXPECT_TRUE(val_res.second);
    std::cerr << (void*)&(*res.first).second << std::endl;
  }

  {
    std::pair<int, del_count> elt(1, del_count());
    auto res = map.insert(elt);
    EXPECT_TRUE(res.second);
    EXPECT_EQ(1, (*res.first).first);
    (*res.first).second.i = 1;
    auto val_res = value_set.insert(&(*res.first).second);
    EXPECT_TRUE(val_res.second);
    std::cerr << (void*)&(*res.first).second << std::endl;
  }

  {
    auto res = map.insert(std::make_pair(0, del_count()));
    EXPECT_FALSE(res.second);
    EXPECT_EQ(0, (*res.first).first);
    EXPECT_EQ(0, (*res.first).second.i);
  }

  del_count::dels = 0;
  map.clear();
  EXPECT_EQ(0, del_count::dels);

  // { // test reuse of val element
  //   auto res = map.insert(std::make_pair(2, del_count()));
  //   EXPECT_TRUE(res.second);
  //   EXPECT_EQ(2, (*res.first).first);
  //   (*res.first).second.i = 2;
  //   auto val_res = value_set.insert(&(*res.first).second);
  //   EXPECT_FALSE(val_res.second);
  //   std::cerr << (void*)&(*res.first).second << std::endl;
  // }

}
} // empty namespace
