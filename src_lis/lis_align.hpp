#ifndef __LIS2_H__
#define __LIS2_H__

#include <vector>
#include <iterator>

namespace lis_align {
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
                                                  std::vector<element<T> >& L, std::vector<unsigned int>& P,
                                                  T a, T b) {
  unsigned int longest = 0, longest_ind = 0;
  const size_t N = std::distance(X, Xend);

  for(size_t i = 0 ; i < N; ++i) {
    element<T>&   e_longest = L[i];
    unsigned int& j_longest = P[i];

    e_longest.len   = 1;
    e_longest.span1 = 0;
    e_longest.span2 = 0;
    j_longest       = N;
    for(size_t j = 0; j < i; ++j) {
      if(X[i].second > X[j].second && e_longest.len < L[j].len + 1) {
        T new_span1 = L[j].span1 + (X[i].first - X[j].first);
        T new_span2 = L[j].span2 + (X[i].second - X[j].second);
        if(new_span1 <= a + b * new_span2 &&
           new_span2 <= a + b * new_span1) {
          e_longest.len   = L[j].len + 1;
          j_longest       = j;
          e_longest.span1 = new_span1;
          e_longest.span2 = new_span2;
        }
      }
    }
    if(longest < e_longest.len) {
      longest = e_longest.len;
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
  std::vector<element<T> > L(N);
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
