#include <gtest/gtest.h>

#include <src_jf_aligner/bounded_counter.hpp>


namespace {
TEST(BoundedCounter, Serial) {
  bounded_counter<false, uint16_t> counter;
  static uint16_t max = 100;

  for(int i = 0; i < 1000; ++i)
    EXPECT_EQ((uint16_t)i < max, counter.inc(max));
  EXPECT_EQ(max, (uint16_t)counter);
} // BoundedCounter.Serial

TEST(BoundedCounter, Parallel) {
  static uint16_t max = 10000;

  {
    bounded_counter<true, uint16_t> counter;

# pragma omp parallel for
    for(int i = 0; i < (uint16_t)max; ++i)
      EXPECT_TRUE(counter.inc(max));

# pragma omp parallel for
    for(int i = 0; i < 10000; ++i)
      EXPECT_FALSE(counter.inc(max));

    EXPECT_EQ(max, (uint16_t)counter);
  }

  {
    bounded_counter<true, uint16_t> counter;
# pragma omp parallel for
    for(int i = 0; i < 10000; ++i)
      counter.inc(max);

    EXPECT_EQ(max, (uint16_t)counter);
  }

  {
    bounded_counter<true, uint16_t> counter;
# pragma omp parallel for
    for(int i = 0; i < 10000; ++i)
      counter.inc(0);

    EXPECT_EQ(max, (uint16_t)10000);
  }

} // BoundedCounter.Parallel


} // empty namespace

