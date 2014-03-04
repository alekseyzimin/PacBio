#ifndef __LIS2_H__
#define __LIS2_H__

#include <vector>
#include <iterator>
#include <type_traits>
#include <algorithm>

#include <src_jf_aligner/lf_forward_list.hpp>


namespace lis_align {
/**
 * A forward list defined for trivial data types. Per thread bulk
 * memory management.
 */
template <typename T>
class forward_list : public lf_forward_list_base<T> {
  typedef lf_forward_list_base<T> super;
  typedef typename super::node    node;

public:
  typedef typename super::iterator       iterator;
  typedef typename super::const_iterator const_iterator;
  typedef T                              value_type;

  forward_list() : super(0) {
    static_assert(std::is_trivial<T>::value, "Forward list defined only for trivial types");
  }
  ~forward_list() { }

  iterator insert_after(const_iterator position, const value_type& val) {
    struct data_page {
      node*  data;
      size_t used;
    };
    static const size_t       node_per_page = 1024 * 1024;
    static __thread data_page thread_data   = { 0, 0 };

    if(!thread_data.data || thread_data.used >= node_per_page) {
      thread_data.used = 0;
      thread_data.data = new node[node_per_page];
      data_pages.push_front(thread_data.data);
    }
    node* cnode = &thread_data.data[thread_data.used++];
    std::copy(&val, &val + 1, &cnode->val_);
    super::insert_after_(position, cnode);
    return iterator(cnode);
  }
  //  iterator insert_after(const_iterator position, value_type&& val) 

private:
  lf_forward_list<node*> data_pages;
};


/**
 * Compute an alignment on an array X where each element is a pair of
 * offsets. It returns the longest alignment (in term of number of
 * elements) the second offset are all increasing and where the spans
 * in first offsets is bounded by an affine relation to the spans of
 * the second offsets (and conversely).
 *
 * The algorithm is quadratic in the length of X.
 */

/**
 * Compute the L and P intermediary arrays, defined as follows:
 *
 * L[i] contains: the length of the longest increasing subsequence ending at
 * element X[i], and the span in first and second offsets of the subsequence.
 *
 * P[i] is the index of the previous element to X[i] in the longest
 * increasing subsequence ending at X[i].
 *
 * a and b define the affine relation. So both (span1 < b * span2 + a)
 * and (span2 < b * span1 + a) are satisfied.
 *
 * The InputIterator is assumed to be a random access iterator with
 * element having first and second members of type T. (e.g. the
 * iterator of the vector type std::vector<std::pair<T, T> >).
 */
template<typename T>
struct element {
  size_t       elt;             // Index of element in X
  unsigned int len;             // Length of LIS
  T            span1;
  T            span2;
};
template<typename T>
std::ostream& operator<<(std::ostream& os, const element<T>& x) {
  return os << "<" << x.len << ", " << x.span1 << ", " << x.span2 << ">";
}
template<typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& x) {
  return os << "<" << x.first << ", " << x.second << ">";
}

template<typename InputIterator, typename T>
std::pair<unsigned int, unsigned int> compute_L_P(const InputIterator X, const InputIterator Xend,
                                                  forward_list<element<T> >& L, std::vector<unsigned int>& P,
                                                  T a, T b) {
  unsigned int longest = 0, longest_ind = 0;
  const size_t N = std::distance(X, Xend);

  for(unsigned int i = 0 ; i < N; ++i) {
    element<T>    e_longest = { i, 1, 0, 0 };
    unsigned int& j_longest = P[i];
    j_longest       = N;

    // Go in decreasing order of subsequence length. Stop at first
    // possible extension, as further subsequence length are too short
    // to improve on this extension.
    auto prev = L.cbefore_begin();
    auto it   = L.cbegin();
    for( ; it != L.cend() && it->len >= e_longest.len; ++it) {
      const unsigned int j = it->elt;
      if(X[i].second > X[j].second && e_longest.len < it->len + 1) {
        T new_span1 = it->span1 + (X[i].first - X[j].first);
        T new_span2 = it->span2 + (X[i].second - X[j].second);
        if(new_span1 <= a + b * new_span2 &&
           new_span2 <= a + b * new_span1) {
          e_longest.len   = it->len + 1;
          j_longest       = j;
          e_longest.span1 = new_span1;
          e_longest.span2 = new_span2;
          break;
        }
      }
      if(it->len < prev->len)
        prev = it;
    }
    L.insert_after(prev, e_longest);
    if(longest < e_longest.len) {
      longest     = e_longest.len;
      longest_ind = i;
    }
  }
  return std::make_pair(longest, longest_ind);
}

template<typename OutputIterator, typename P_type>
void indices_reversed(P_type& P, const unsigned int len, unsigned int start, OutputIterator out) {
  for(unsigned int i = 0; i < len; ++i, ++out, start = P[start])
    *out = start;
}

template<typename InputIterator, typename T>
unsigned int indices(const InputIterator X, const InputIterator Xend, std::vector<unsigned int>& res,
                     T a, T b) {
  const size_t N = std::distance(X, Xend);
  //  std::vector<element<T> > L(N);
  forward_list<element<T> > L;
  std::vector<unsigned int> P(N);
  const std::pair<unsigned int, unsigned int> lis = compute_L_P(X, Xend, L, P, a, b);

  if(res.size() < lis.first)
    res.resize(lis.first);
  indices_reversed(P, lis.first, lis.second, res.rbegin());
  return lis.first;
}

template<typename InputIterator, typename T>
std::vector<unsigned int> indices(const InputIterator X, const InputIterator Xend,
                                  T a, T b) {
  std::vector<unsigned int> res;

  indices(X, Xend, res, a, b);
  return res;
}
} // namespace lis2


#endif /* __LIS2_H__ */
