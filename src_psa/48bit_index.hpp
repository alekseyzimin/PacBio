/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef __48BIT_INDEX_H__
#define __48BIT_INDEX_H__

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include "48bit_iterator.hpp"

template<typename IDX>
struct fortyeight_index {
  size_t    m_size;
  uint32_t* m_base32;
  uint16_t* m_base16;

  fortyeight_index(size_t s)
    : m_size(s)
    , m_base32(new uint32_t[(s * 3 + 1) / 2 + 3])
    , m_base16((uint16_t*)(m_base32 + s))
  { }

  ~fortyeight_index() {
    delete [] m_base32;
  }

  typedef fortyeight_iterator<IDX>       iterator;
  typedef const_fortyeight_iterator<IDX> const_iterator;

  const_iterator begin() const { return const_iterator(m_base32, m_base16); }
  iterator begin() { return iterator(m_base32, m_base16); }
  const_iterator end() const { return begin() + m_size; }
  iterator end() { return begin() + m_size; }
  const_iterator cbegin() const { return begin(); }
  const_iterator cend() const { return end(); }

  auto operator[](size_t i) -> decltype(this->begin()[i]) { return begin()[i]; }
  IDX operator[](size_t i) const { return cbegin()[i]; }
};

#endif /* __48BIT_INDEX_H__ */
