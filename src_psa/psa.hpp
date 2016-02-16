/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef __PSA_H__
#define __PSA_H__

#include <memory>
#include <thread>

#include <boost/asio/coroutine.hpp>
#include <boost/asio/yield.hpp>
#include "compact_index.hpp"
#include "48bit_index.hpp"
#include "mer_sa.hpp"
#include "prefetch_iterator_traits.hpp"

template<typename CHARPTR, typename SAIDX = uint64_t>
class PSA {
  typedef mer_sa_imp::SA<CHARPTR, typename fortyeight_index<SAIDX>::iterator> SA;
  typedef typename fortyeight_index<SAIDX>::iterator                          SAIDPTR;
  typedef typename fortyeight_index<SAIDX>::const_iterator                    CSAIDPTR;
  typedef typename compactsufsort_imp::prefetch_iterator_traits<uint64_t*>    counts_prefetch;
  typedef typename compactsufsort_imp::prefetch_iterator_traits<CSAIDPTR>     CSAID_prefetch;

  template<typename CHARPTR2>
  struct search_element {
    search_element(const PSA& psa) : m_psa(psa) { }
    bool is_complete() const { return m_coro.is_complete(); }
    // Start a new search in PSA for pattern P. Only use if Psize <= PSA.min_size
    void start_small(const CHARPTR2 P, SAIDX Psize);
    // Call advance while .is_complete() returns false. Then results
    // are in index and nb
    void advance_small();

    void start_large(const CHARPTR2 P, SAIDX Psize);
    void advance_large();

    SAIDX nb, index;

  protected:
    boost::asio::coroutine m_coro;
    const PSA&             m_psa;
    uint64_t               m_mer, m_left_mer, m_start_mer, m_end_mer;
    unsigned int           m_min_size, m_left;
    CSAIDPTR               m_first, m_it, m_last;
    SAIDX                  m_step;

    // Return mer and length at this position it
    std::pair<uint64_t, SAIDX> get_mer(CSAIDPTR it) const {
      const SAIDX pos = *it + m_psa.min_size();
      const SAIDX cs  = pos < m_psa.size() ? std::min((SAIDX)m_left, m_psa.size() - pos) : 0;
      return std::make_pair(mer_sa_imp::str_to_mer(m_psa.m_T + pos, cs), cs);
    }

    // Compare mer returned by get_mer with left_mer. Return -1, 0, or
    // 1 if mer is less, equal or greater than left_mer.
    int compare_mer(const std::pair<uint64_t, SAIDX>& mer) const {
      if(m_left == mer.second)
        return mer.first < m_left_mer ? -1 : (mer.first > m_left_mer ? 1 : 0);
      // Only possibility is mer.second < left_mer. Then mer and
      // left_mer cannot be equal.
      const uint64_t lmer = m_left_mer >> (2 * (m_left - mer.second));
      return mer.first <= lmer ? -1 : 1;
    }
  };

  template<typename CHARPTR2>
  friend struct search_element;

  template<typename Iterator>
  struct search_iterator : boost::asio::coroutine {
    // Number of searches to carry out in "parallel"
    static const int                                            nb_searches = 8;
    typedef typename std::iterator_traits<Iterator>::value_type pattern_type;

    search_iterator(PSA& psa, Iterator first, Iterator last, SAIDX Psize)
      : m_it(first)
      , m_last(last)
      , m_Psize(Psize)
      , m_searches{ psa, psa, psa, psa, psa, psa, psa, psa }
      , m_live(1)
      , m_is_small(Psize <= psa.min_size())
    { }

    bool next() {
      reenter(this) for( ; m_live > 0; ) {
        // Initialize the searches
        for(m_live = 0; m_it != m_last && m_live < nb_searches; ++m_live, ++m_it)
          m_is_small ? m_searches[m_live].start_small(*m_it, m_Psize) : m_searches[m_live].start_large(*m_it, m_Psize);
        bool done;
        do { // Advances all searches in parallel
          done = true;
          for(m_cur = 0; m_cur < m_live; ++m_cur) {
            if(!m_searches[m_cur].is_complete()) {
              done = false;
              m_is_small ? m_searches[m_cur].advance_small() : m_searches[m_cur].advance_large();
            }
          }
        } while(!done);
        for(m_cur = 0; m_cur < m_live; ++m_cur) { // Return searches 1 by 1
          nb    = m_searches[m_cur].nb;
          index = m_searches[m_cur].index;
          yield return true;
        }
      }
      return false;
    }

    SAIDX nb;
    SAIDX index;

    Iterator                     m_it;
    const Iterator               m_last;
    const SAIDX                  m_Psize;
    search_element<pattern_type> m_searches[nb_searches];
    int                          m_live; // Number of searches that are alive
    int                          m_cur; // Current search result
    bool                         m_is_small;
  };
  template<typename Iterator>
  friend struct search_iterator;

public:
  PSA(CHARPTR T, SAIDX n, unsigned int min_size = 10, int max_size = 15, unsigned threads = std::thread::hardware_concurrency())
    : m_T(T)
    , m_n(n)
    , m_min_size(min_size)
    , m_max_size(max_size)
    , m_nb_counts(((size_t)1 << (2 * min_size)) + 1)
    , m_mer_counts(new uint64_t[m_nb_counts])
    , m_psa(n - min_size + 1)
  {
    mer_sa::create_mt(m_T, m_psa.begin(), m_n, threads, m_mer_counts.get(), m_min_size, m_max_size);
  }

  bool check() const {
    return mer_sa::check(m_T, m_n, m_psa.begin(), m_n - m_min_size + 1, m_mer_counts.get(), m_min_size, m_max_size);
  }

  bool check_suffixes() const {
    return mer_sa::check_suffixes(m_T, m_n, m_psa.begin(), m_n - m_min_size + 1, m_mer_counts.get(), m_min_size, m_max_size);
  }

  template<typename CHARPTR2>
  std::pair<SAIDX, SAIDX> search(const CHARPTR2 P, SAIDX Psize) const {
    return mer_sa::search(m_T, m_n, m_psa.begin(), m_n, m_mer_counts.get(), m_min_size, m_max_size, P, Psize);
  }

  template<typename Iterator>
  search_iterator<Iterator> multi_search(Iterator first, Iterator last, SAIDX Psize) const {
    return search_iterator<Iterator>(*this, first, last, Psize);
  }

  SAIDX operator[](size_t n) const { return m_psa[n]; }

  unsigned int min_size() const { return m_min_size; }
  unsigned int max_size() const { return m_max_size; }
  SAIDX size() const { return m_n; }
  void* sa_ptr() { return m_psa.get(); }

private:
  const CHARPTR              m_T;
  const SAIDX                m_n;
  const unsigned int         m_min_size, m_max_size;
  const size_t               m_nb_counts;
  const std::unique_ptr<uint64_t[]>      m_mer_counts;
  fortyeight_index<SAIDX>    m_psa;
};

// Implementation of search_element methods
template<typename CHARPTR, typename SAIDX>
template<typename CHARPTR2>
void PSA<CHARPTR, SAIDX>::search_element<CHARPTR2>::start_small(const CHARPTR2 P, SAIDX Psize) {
  assert(Psize <= m_psa.min_size());
  m_coro = boost::asio::coroutine();
  m_left = 2 * (m_psa.min_size() - Psize);
  m_mer  = mer_sa_imp::str_to_mer(P, Psize);
}

template<typename CHARPTR, typename SAIDX>
template<typename CHARPTR2>
void PSA<CHARPTR, SAIDX>::search_element<CHARPTR2>::advance_small() {
  reenter(m_coro) {
    counts_prefetch::read(&m_psa.m_mer_counts[m_mer << m_left]);
    yield;
    index = m_psa.m_mer_counts[m_mer << m_left];
    counts_prefetch::read(&m_psa.m_mer_counts[~(~m_mer << m_left) + 1]);
    yield;
    nb    = m_psa.m_mer_counts[~(~m_mer << m_left) + 1] - index;
  }
}

template<typename CHARPTR, typename SAIDX>
template<typename CHARPTR2>
void PSA<CHARPTR, SAIDX>::search_element<CHARPTR2>::start_large(const CHARPTR2 P, SAIDX Psize) {
  assert(Psize > m_psa.min_size());
  m_coro     = boost::asio::coroutine();
  m_left     = Psize - m_psa.min_size();
  m_mer      = mer_sa_imp::str_to_mer(P, m_psa.min_size());
  m_left_mer = mer_sa_imp::str_to_mer(P + m_psa.min_size(), m_left);
  index      = 0;
}

template<typename CHARPTR, typename SAIDX>
template<typename CHARPTR2>
void PSA<CHARPTR, SAIDX>::search_element<CHARPTR2>::advance_large() {
  reenter(m_coro) {
    counts_prefetch::read(&m_psa.m_mer_counts[m_mer]);
    yield;
    nb = m_psa.m_mer_counts[m_mer + 1] - m_psa.m_mer_counts[m_mer];
    if(nb == 0) { yield break; }
    m_first     = m_psa.m_psa.cbegin() + m_psa.m_mer_counts[m_mer];
    CSAID_prefetch::read(m_first);
    yield;
    m_start_mer = get_mer(m_first).first;
    CSAID_prefetch::read(m_first + nb - 1);
    yield;
    m_end_mer   = get_mer(m_first + nb - 1).first;
    while(nb > 0) { // Find mer in psa
      m_step = std::min(nb - 1,
                        (SAIDX)std::lrint(nb * (double)(2 * (m_left_mer - m_start_mer) + 1) /
                                          (double)(2 * (m_end_mer - m_start_mer + 1))));
      m_it   = m_first + m_step;
      CSAID_prefetch::read(m_it);
      yield;
      const auto cmer = get_mer(m_it);
      const int cmp = compare_mer(cmer);
      if(cmp < 0) {
        m_first      = ++m_it;
        nb          -= m_step + 1;
        m_start_mer  = cmer.first + 1;
      } else if(cmp > 0) {
        nb        = m_step;
        m_end_mer = cmer.first;
      } else {
        break;
      }
    }
    if(nb == 0) { yield break; }
    index = nb - m_step - 1;
    m_last = m_first + m_step + 1;
    m_mer = m_left_mer;
    for(nb = m_step; nb > 0; ) { // Find lower range
      m_step = std::min(nb - 1,
                        (SAIDX)std::rint(nb * (double)(m_left_mer - m_start_mer) / (double)(m_mer - m_start_mer + 1)));
      m_it   = m_first + m_step;
      CSAID_prefetch::read(m_it);
      yield;
      const auto cmer = get_mer(m_it);
      if(compare_mer(cmer) < 0) {
        m_first      = ++m_it;
        nb          -= m_step + 1;
        m_start_mer  = cmer.first;
      } else {
        nb    = m_step;
        m_mer = cmer.first;
      }
    }

    m_mer = m_left_mer;
    for(nb = index; nb > 0; ) { // Find upper range
      m_step = std::min(nb - 1,
                        (SAIDX)std::lrint(nb * (double)(m_left_mer + 1 - m_mer) / (double)(m_end_mer - m_mer + 1)));
      m_step = std::min(nb - 1, m_step);
      m_it   = m_last + m_step;
      CSAID_prefetch::read(m_it);
      yield;
      const auto cmer = get_mer(m_it);
      if(compare_mer(cmer) <= 0) {
        m_last  = ++m_it;
        nb     -= m_step + 1;
        m_mer   = cmer.first;
      } else {
        nb        = m_step;
        m_end_mer = cmer.first;
      }
    }
    nb    = m_last - m_first;
    index = m_first - m_psa.m_psa.cbegin();
  }
}

#endif /* __PSA_H__ */
