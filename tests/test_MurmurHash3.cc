#include <gtest/gtest.h>
#include <MurmurHash3.h>

namespace {
template<typename T>
void add_elt(MurmurHash3A& h, char* data, int& len) {
  T x = random();
  h.add(x);
  memcpy(data + len, &x, sizeof(T));
  len += sizeof(T);
}

TEST(MurmurHash3A, Incremental) {
  static const int nb_elt     = 100;
  static const int iterations = 100;

  for(int j = 0; j < iterations; ++j) {
    union {
      uint64_t       b64[nb_elt];
      char           b8[nb_elt * sizeof(uint64_t)];
    } buffer;
    int              len    = 0;
    const int        seed   = random();
    MurmurHash3A     hash(seed);

    // Add many random elements of random size
    for(int i = 0; i < nb_elt; ++i) {
      switch(random() % 4) {
      case 0: add_elt<uint8_t>(hash, buffer.b8, len);  break;
      case 1: add_elt<uint16_t>(hash, buffer.b8, len); break;
      case 2: add_elt<uint32_t>(hash, buffer.b8, len); break;
      case 3: add_elt<uint64_t>(hash, buffer.b8, len); break;
      }
    }

    uint64_t r1, r2;
    hash.finalize(r1, r2);
    uint64_t out[2];
    MurmurHash3_x64_128(buffer.b64, len, seed, out);
    EXPECT_EQ(out[0], r1);
    EXPECT_EQ(out[1], r2);
  }
}
} // empty namespace
