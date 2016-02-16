/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef __COMPACT_DNA_H__
#define __COMPACT_DNA_H__

#include "compact_iterator.hpp"
#include <string>
#include <memory>


struct compact_dna {
  static constexpr const char* bases = "ACGT";
  typedef uint64_t W;
  const size_t size;
  std::unique_ptr<W[]> mem;

  compact_dna(size_t s) : size(s), mem(new W[s / 32 + (s % 32 != 0)]) { }

  typedef compact_iterator<uint8_t, W>       iterator;
  typedef const_compact_iterator<uint8_t, W> const_iterator;

  static iterator iterator_at(W* ptr) { return iterator(ptr, 2, 0); }
  static const_iterator iterator_at(const W* ptr) { return const_iterator(ptr, 2, 0); }
  static const_iterator const_iterator_at(const W* ptr) { return const_iterator(ptr, 2, 0); }

  class const_base_iterator : public std::iterator<std::random_access_iterator_tag, char> {
    const_iterator m_it;
  public:
    const_base_iterator() = default;
    const_base_iterator(const_iterator it) : m_it(it) { }
    const_base_iterator(const const_base_iterator& rhs) : m_it(rhs.m_it) { }

    const_base_iterator& operator++() { ++m_it; return *this; }
    const_base_iterator operator++(int) { const_base_iterator res(*this); operator++(); return res; }
    const_base_iterator& operator--() { --m_it; return *this; }
    const_base_iterator operator--(int) { const_base_iterator res(*this); operator--(); return res; }
    bool operator==(const const_base_iterator& rhs) const { return m_it == rhs.m_it; }
    bool operator!=(const const_base_iterator& rhs) const { return m_it != rhs.m_it; }
    bool operator<(const const_base_iterator& rhs) const { return m_it < rhs.m_it; }
    bool operator<=(const const_base_iterator& rhs) const { return m_it <= rhs.m_it; }
    bool operator>(const const_base_iterator& rhs) const { return m_it > rhs.m_it; }
    bool operator>=(const const_base_iterator& rhs) const { return m_it >= rhs.m_it; }
    const_base_iterator& operator+=(ptrdiff_t n) { m_it += n; return *this; }
    const_base_iterator& operator-=(ptrdiff_t n) { m_it -= n; return *this; }
    const_base_iterator operator+(ptrdiff_t n) const { const_base_iterator res(*this); return res += n; }
    const_base_iterator operator-(ptrdiff_t n) const { const_base_iterator res(*this); return res -= n; }
    ptrdiff_t operator-(const const_base_iterator& rhs) const { return m_it - rhs.m_it; }
    char operator*() const { return bases[*m_it]; }
    char operator[](size_t i) const { return bases[m_it[i]]; }
    friend void swap(const_base_iterator& a, const_base_iterator& b);
  };

  static const_base_iterator const_base_iterator_at(const W* ptr) { return const_base_iterator(const_iterator_at(ptr)); }

  const_iterator begin() const { return const_iterator(mem.get(), 2, 0); }
  iterator begin() { return iterator(mem.get(), 2, 0); }
  const_iterator end() const { return begin() + size; }
  iterator end() { return begin() + size; }
  const_iterator cbegin() const { return begin(); }
  const_iterator cend() const { return end(); }

  // Copy DNA-ASCII from <s> to <to>. Returns the number of bases
  // copied. If the value is less than the size of s, an error
  // occurred at that character.
  static size_t copy_from_str(iterator to, const std::string& s) {
    return copy_from_str(to, s.c_str(), s.size());
  }
  static size_t copy_from_str(iterator to, const char* start, size_t n) {
    return copy_from_str(to, start, start + n);
  }
  // static size_t copy_from_str(iterator to, const char* start, const char* end) {
  //   unsigned char c;
  //   auto p = start;
  //   for( ; p < end; ++p, ++to) {
  //     switch(*p) {
  //     case 'a': case 'A': c = 0; break;
  //     case 'c': case 'C': c = 1; break;
  //     case 'g': case 'G': c = 2; break;
  //     case 't': case 'T': c = 3; break;
  //     default: goto done;
  //     }
  //     *to = c;
  //   }
  // done:
  //   return p - start;
  // }

  static void copy_from_str_slow(iterator& to, const char* start, const char* const end) {
    char c = 0;
    for( ; start != end; ++to, ++start) {
      switch(*start) {
      case 'a': case 'A': c = 0; break;
      case 'c': case 'C': c = 1; break;
      case 'g': case 'G': c = 2; break;
      case 't': case 'T': c = 3; break;
      }
      *to = c;
    }
  }

  inline static W char_to_code8(W x) {
    x = ((x >> 1) ^ (x >> 2)) & 0x303030303030303UL;
    x |= (x >> 6);
    x |= (x >> 12);
    return (x & 0xff) | ((x >> 24) & 0xff00);
  }

  static size_t copy_from_str(iterator to, const char* const start, const char* const end) {
    // p is first W aligned pointer after start
    const W* p = reinterpret_cast<W*>((reinterpret_cast<std::uintptr_t>(start) + alignof(W) - 1) & -alignof(W));
    copy_from_str_slow(to, start, std::min((const char*)p, end));
    const W* pend = p + (end - (char*)p) / sizeof(W);
    W x;
    while(p + 4 <= pend) {
      x  = char_to_code8(*p);       ++p;
      x |= char_to_code8(*p) << 16; ++p;
      x |= char_to_code8(*p) << 32; ++p;
      x |= char_to_code8(*p) << 48; ++p;
      to.set_bits(x, 8 * sizeof(uint64_t));
      to += 32;
    }

    x = 0;
    unsigned int l = pend - p;
    switch(l) {
    case 3: x  = char_to_code8(*p);                   ++p;
    case 2: x |= char_to_code8(*p) << (16 * (l - 2)); ++p;
    case 1: x |= char_to_code8(*p) << (16 * (l - 1)); ++p;
      to.set_bits(x, 16 * l);
      to += l * 8;
    }

    copy_from_str_slow(to, (char*)p, end);
    return end - start;
  }

  template<typename Iterator>
  static void copy_to_str(Iterator to, const_iterator start, const_iterator end) {
    for( ; start != end; ++start, ++to)
      *to = bases[*start];
  }

  inline static W swap_word(W w) {
    return ((w & 0x5555555555555555UL) << 1) | ((w & 0xaaaaaaaaaaaaaaaaUL) >> 1);
  }
  inline static bool compare_swap_words(W w1, W w2) {
    w1 = swap_word(w1);
    w2 = swap_word(w2);
    W bmask = w1 ^ w2;
    bmask &= -bmask;
    return (w1 & bmask) == 0;
  }

  // static bool lexicographical_compare(const const_iterator first1, const const_iterator last1,
  //                                     const const_iterator first2, const const_iterator last2) {
  //   static const unsigned int Wbits = 8 * sizeof(W);
  //   static const unsigned int bits = 2;
  //   const W*           f1   = first1.get_ptr();
  //   const W* const     l1   = last1.get_ptr() + (last1.get_offset() + bits >= Wbits);
  //   unsigned int       o1   = first1.get_offset();
  //   const W*           f2   = first2.get_ptr();
  //   const W* const     l2   = last2.get_ptr() + (last2.get_offset() + bits >= Wbits);
  //   unsigned int       o2   = first2.get_offset();

  //   while(f1 < l1 && f2 < l2) {
  //     const W mask = ~(W)0 >> std::max(o1, o2);
  //     const W w1 = (*f1 >> o1) & mask;
  //     const W w2 = (*f2 >> o2) & mask;
  //     if(w1 != w2) return compare_swap_words(w1, w2);
  //     o1 += bits;
  //     if(o1 >= Wbits) {
  //       ++f1;
  //       o1 -= Wbits;
  //     }
  //     o2 += bits;
  //     if(o2 >= Wbits) {
  //       ++f2;
  //       o2 -= Wbits;
  //     }
  //   }

  //   const unsigned int len1 = (f1 == l1 ? last1.get_offset() : Wbits) - o1;
  //   const unsigned int len2 = (f2 == l2 ? last2.get_offset() : Wbits) - o2;
  //   if(len1 == 0 || len2 == 0) return len1 == 0 && len2 != 0;
  //   const W            mask = ~(W)0 >> (Wbits - std::min(len1, len2));
  //   const W            w1   = (*f1 >> o1) & mask;
  //   const W            w2   = (*f2 >> o2) & mask;

  //   return w1 == w2 ? len1 < len2 : compare_swap_words(w1, w2);
  // }
};

inline compact_dna::const_base_iterator operator+(ptrdiff_t n, const compact_dna::const_base_iterator& rhs) { return rhs + n; }
inline void swap(compact_dna::const_base_iterator& a, compact_dna::const_base_iterator& b) { swap(a.m_it, b.m_it); }

#endif /* __COMPACT_DNA_H__ */
