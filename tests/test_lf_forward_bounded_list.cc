#include <gtest/gtest.h>
#include <src_jf_aligner/lf_forward_bounded_list.hpp>

namespace {
TEST(LfForwardBoundedList, LockFree) {
  typedef lf_forward_bounded_list<int>                                       list_type;
  typedef std::iterator_traits<list_type::iterator>::difference_type distance_type;

  static const int size = 10000;
  list_type        list;
  std::cout << sizeof(list) << ' ' << sizeof(lf_forward_list_base<int, true>) << std::endl;

  EXPECT_TRUE(list.cbegin() == list.cend());
  EXPECT_FALSE(list.cbegin() != list.cend());
  EXPECT_EQ((distance_type)0, std::distance(list.cbegin(), list.cend()));

  bool it_equal = true;
#pragma omp parallel for reduction(&&:it_equal)
  for(int i = 0; i < size; ++i) {
    auto it = list.push_front(i);
    it_equal = it_equal && (*it == i);
  }
  EXPECT_TRUE(it_equal);

  EXPECT_FALSE(list.cbegin() == list.cend());
  EXPECT_TRUE(list.cbegin() != list.cend());
  EXPECT_EQ((distance_type)size, std::distance(list.cbegin(), list.cend()));

  std::set<int> unique;
  for(auto it = list.cbegin(); it != list.cend(); ++it) {
    auto res = unique.insert(*it);
    EXPECT_TRUE(res.second); // Means it was not present in the set
    EXPECT_TRUE(*it >= 0 && *it < size);
  }
  EXPECT_EQ((std::set<int>::size_type)size, unique.size());

  list.clear();
  EXPECT_EQ((distance_type)0, std::distance(list.cbegin(), list.cend()));
}


TEST(LfForwardBoundedList, SingleThread) {
  typedef lf_forward_bounded_list<int, false>                                list_type;
  typedef std::iterator_traits<list_type::iterator>::difference_type distance_type;

  static const int size = 10000;
  list_type        list;
  std::cout << sizeof(list) << std::endl;

  EXPECT_TRUE(list.cbegin() == list.cend());
  EXPECT_FALSE(list.cbegin() != list.cend());
  EXPECT_EQ((distance_type)0, std::distance(list.cbegin(), list.cend()));

  bool it_equal = true;
  for(int i = 0; i < size; ++i) {
    auto it = list.push_front(i);
    it_equal = it_equal && (*it == i);
  }
  EXPECT_TRUE(it_equal);

  EXPECT_FALSE(list.cbegin() == list.cend());
  EXPECT_TRUE(list.cbegin() != list.cend());
  EXPECT_EQ((distance_type)size, std::distance(list.cbegin(), list.cend()));

  std::set<int> unique;
  for(auto it = list.cbegin(); it != list.cend(); ++it) {
    auto res = unique.insert(*it);
    EXPECT_TRUE(res.second); // Means it was not present in the set
    EXPECT_TRUE(*it >= 0 && *it < size);
  }
  EXPECT_EQ((std::set<int>::size_type)size, unique.size());

  list.clear();
  EXPECT_EQ((distance_type)0, std::distance(list.cbegin(), list.cend()));
}

}
