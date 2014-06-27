#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <random>
#include <gtest/gtest.h>
#include <src_jf_aligner/union_find.hpp>

namespace {
TEST(UnionFind, Set) {
  union_find::set s;
  EXPECT_EQ(0, s.rank);
  EXPECT_EQ(&s, s.parent);
  EXPECT_EQ(&s, s.root());

  union_find::union_sets(s, s);
  EXPECT_EQ(0, s.rank);
  EXPECT_EQ(&s, s.parent);

  union_find::set s2;
  EXPECT_NE(s2.root(), s.root());
  union_find::union_sets(s, s2);
  EXPECT_EQ(s2.root(), s.root());
  EXPECT_TRUE((s.rank == 0 && s2.rank == 1 && s.root() == &s2) ||
              (s.rank == 1 && s2.rank == 0 && s2.root() == &s));
} // UnionFind.Set



struct elt {
  union_find::set component_;
  int i_;
  //  explicit elt(int i) : i_(i) { }
};

template<typename Iterator>
std::map<union_find::set*, std::set<int> > find_components(Iterator start, Iterator end) {
  std::map<union_find::set*, std::set<int> > components;
  EXPECT_EQ((size_t)0, components.size());
  for(auto it = start; it != end; ++it)
    components[it->component_.root()].insert(it->i_);
  return components;
}


TEST(UnionFind, Test) {
  static const int size = 100;
  static const int comp = 15;
  std::vector<elt> elements(size);
  std::vector<int> order;
  auto rng = std::default_random_engine();
  for(int i = 0; i < size; ++i) {
    elements[i].i_ = i;
    order.push_back(i);
  }
  shuffle(order.begin(), order.end(), rng);

  {
    auto components = find_components(elements.begin(), elements.end());
    EXPECT_EQ(size, components.size());
  }

  for(int i = 0; i < size; ++i) {
    const int j = order[i];
    if(j % comp != (comp - 1) && j < size -1)
      elements[j].component_ |= elements[j+1].component_;
  }

  auto components = find_components(elements.begin(), elements.end());
  EXPECT_EQ(size / comp + (size % comp > 0), components.size());

  for(auto& c : components) {
    auto& elements = c.second;
    auto it = elements.cbegin();
    EXPECT_EQ(0, *it % comp);
    for(auto pit = it++; it != elements.cend(); pit = it++)
      EXPECT_EQ(*pit + 1, *it);
  }
}
} // empty namespace
