//-----------------------------------------------------------------------------
// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.

#ifndef _MURMURHASH3_H_
#define _MURMURHASH3_H_

#include <string.h>

//-----------------------------------------------------------------------------
// Platform-specific functions and macros

// Microsoft Visual Studio

#if defined(_MSC_VER)

typedef unsigned char uint8_t;
typedef unsigned long uint32_t;
typedef unsigned __int64 uint64_t;

#define FORCE_INLINE	__forceinline

#include <stdlib.h>

#define ROTL32(x,y)	_rotl(x,y)
#define ROTL64(x,y)	_rotl64(x,y)

#define BIG_CONSTANT(x) (x)


// Other compilers

#else	// defined(_MSC_VER)

#include <stdint.h>


#define	FORCE_INLINE inline __attribute__((always_inline))

template<typename T>
T rotl(T x, int8_t r) {
  return (x << r) | (x >> (sizeof(T) * 8 - r));
}

#define	ROTL32(x,y)	rotl<uint32_t>(x,y)
#define ROTL64(x,y)	rotl<uint64_t>(x,y)

#define BIG_CONSTANT(x) (x##LLU)

#endif // !defined(_MSC_VER)

//-----------------------------------------------------------------------------

void MurmurHash3_x86_32  ( const void * key, int len, uint32_t seed, void * out );

void MurmurHash3_x86_128 ( const void * key, int len, uint32_t seed, void * out );

void MurmurHash3_x64_128 ( const void * key, int len, uint32_t seed, void * out );

// Generic Mumur hash implementation. Key must support the subscripte
// operator: uint64_t operator[](unsigned int i), where i indexes the
// ith 64bit word. The len parameter is still the length in bytes (not
// in uint64_t).
template<typename T>
void MurmurHash3_T_128(const T& key, const int len, const uint32_t seed, void * out );
//-----------------------------------------------------------------------------

template<typename T>
void MurmurHash3_T_128(const T& key, const int len,
                       const uint32_t seed, void * out )
{
  const int nblocks = len / 16;
  uint64_t  h1      = seed;
  uint64_t  h2      = seed;
  uint64_t  c1      = BIG_CONSTANT(0x87c37b91114253d5);
  uint64_t  c2      = BIG_CONSTANT(0x4cf5ad432745937f);

  //----------
  // body

  for(int i = 0; i < nblocks; i++)
  {
    uint64_t k1 = key[i*2+0];
    uint64_t k2 = key[i*2+1];

    k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;

    h1 = ROTL64(h1,27); h1 += h2; h1 = h1*5+0x52dce729;

    k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

    h2 = ROTL64(h2,31); h2 += h1; h2 = h2*5+0x38495ab5;
  }

  //----------
  // tail

  const int remainder = len & 15;
  if(remainder) {
    uint64_t tail_words[2];
    tail_words[0] = key[nblocks*2];
    if(remainder > 8)
      tail_words[1] = key[nblocks*2+1];
    const uint8_t* tail = (const uint8_t*)tail_words;
    
    uint64_t k1 = 0;
    uint64_t k2 = 0;
    
    switch(remainder) {
    case 15: k2 ^= uint64_t(tail[14]) << 48;
    case 14: k2 ^= uint64_t(tail[13]) << 40;
    case 13: k2 ^= uint64_t(tail[12]) << 32;
    case 12: k2 ^= uint64_t(tail[11]) << 24;
    case 11: k2 ^= uint64_t(tail[10]) << 16;
    case 10: k2 ^= uint64_t(tail[ 9]) << 8;
    case  9: k2 ^= uint64_t(tail[ 8]) << 0;
             k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;
             
    case  8: k1 ^= uint64_t(tail[ 7]) << 56;
    case  7: k1 ^= uint64_t(tail[ 6]) << 48;
    case  6: k1 ^= uint64_t(tail[ 5]) << 40;
    case  5: k1 ^= uint64_t(tail[ 4]) << 32;
    case  4: k1 ^= uint64_t(tail[ 3]) << 24;
    case  3: k1 ^= uint64_t(tail[ 2]) << 16;
    case  2: k1 ^= uint64_t(tail[ 1]) << 8;
    case  1: k1 ^= uint64_t(tail[ 0]) << 0;
             k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
    };
  }

  //----------
  // finalization

  h1 ^= len; h2 ^= len;

  h1 += h2;
  h2 += h1;

  // Manually inlined
  // h1 = fmix(h1);
  // h2 = fmix(h2);
  h1 ^= h1 >> 33;
  h1 *= BIG_CONSTANT(0xff51afd7ed558ccd);
  h1 ^= h1 >> 33;
  h1 *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  h1 ^= h1 >> 33;

  h2 ^= h2 >> 33;
  h2 *= BIG_CONSTANT(0xff51afd7ed558ccd);
  h2 ^= h2 >> 33;
  h2 *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  h2 ^= h2 >> 33;

  h1 += h2;
  h2 += h1;

  ((uint64_t*)out)[0] = h1;
  ((uint64_t*)out)[1] = h2;
}

// Incremental version of MurmurHash3. Only support adding 32bits or
// 64 bits word at a time.
class MurmurHash3A {
  static const uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
  static const uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);
  uint64_t              h1, h2;
  union {
    uint8_t             b8[3 * sizeof(uint64_t)];
    uint64_t            b64[3];
  } buffer;
  unsigned int          filled;
  unsigned int          len;

public:
  explicit MurmurHash3A(const uint32_t seed = 590632053) :
    h1(seed), h2(seed), filled(0), len(0)
  { }

  // Add a new word of type T
  template<typename T>
  void add(const T k) {
    memcpy(buffer.b8 + filled, &k, sizeof(T));
    filled += sizeof(T);
    len    += sizeof(T);
    if(filled >= 2 * sizeof(uint64_t)) {
      add__(buffer.b64[0], buffer.b64[1]);
      filled -= 2 * sizeof(uint64_t);
      memcpy(buffer.b8, buffer.b8 + 2 * sizeof(uint64_t), filled);
    }
  }

  // Finalize the value of the hash. The hash object should not be
  // used after this call.
  void finalize(uint64_t& out0, uint64_t& out1) {
    const int remainder = len & 0xf;
    uint64_t k1 = 0;
    uint64_t k2 = 0;

    switch(remainder) {
    case 15: k2 ^= uint64_t(buffer.b8[14]) << 48;
    case 14: k2 ^= uint64_t(buffer.b8[13]) << 40;
    case 13: k2 ^= uint64_t(buffer.b8[12]) << 32;
    case 12: k2 ^= uint64_t(buffer.b8[11]) << 24;
    case 11: k2 ^= uint64_t(buffer.b8[10]) << 16;
    case 10: k2 ^= uint64_t(buffer.b8[ 9]) << 8;
    case  9: k2 ^= uint64_t(buffer.b8[ 8]) << 0;
      h2 ^= ROTL64(k2 * c2, 33) * c1;
      //        k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

    case  8: k1 ^= uint64_t(buffer.b8[ 7]) << 56;
    case  7: k1 ^= uint64_t(buffer.b8[ 6]) << 48;
    case  6: k1 ^= uint64_t(buffer.b8[ 5]) << 40;
    case  5: k1 ^= uint64_t(buffer.b8[ 4]) << 32;
    case  4: k1 ^= uint64_t(buffer.b8[ 3]) << 24;
    case  3: k1 ^= uint64_t(buffer.b8[ 2]) << 16;
    case  2: k1 ^= uint64_t(buffer.b8[ 1]) << 8;
    case  1: k1 ^= uint64_t(buffer.b8[ 0]) << 0;
      h1 ^= ROTL64(k1 * c1, 31) * c2;
      //        k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
    }
    h1 ^= len; h2 ^= len;
    h1 += h2;  h2 += h1;
    h1 = (h1 ^ (h1 >> 33)) * BIG_CONSTANT(0xff51afd7ed558ccd);
    h1 = (h1 ^ (h1 >> 33)) * BIG_CONSTANT(0xc4ceb9fe1a85ec53);
    h1 ^= h1 >> 33;
    h2 = (h2 ^ (h2 >> 33)) * BIG_CONSTANT(0xff51afd7ed558ccd);
    h2 = (h2 ^ (h2 >> 33)) * BIG_CONSTANT(0xc4ceb9fe1a85ec53);
    h2 ^= h2 >> 33;

    out0 = h1 = h1 + h2;
    out1 = h2 = h2 + h1;
  }

  void add__(const uint64_t k1, const uint64_t k2) {
    h1 ^= ROTL64(k1 * c1, 31) * c2;
    h1  = (ROTL64(h1, 27) + h2) * 5 + 0x52dce729;
    h2 ^= ROTL64(k2 * c2, 33) * c1;
    h2  = (ROTL64(h2, 31) + h1) * 5 + 0x38495ab5;
  }

};


#endif // _MURMURHASH3_H_
