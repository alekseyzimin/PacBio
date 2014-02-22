#ifndef __LIS_HPP__
#define __LIS_HPP__

#include <iterator>
#include <algorithm>

namespace lis {
/* Compute the Longest Increasing Subsequence on the data in is the
   array X, given by the start and end random access iterators, and
   contains N elements. It outputs the length of the longest
   increasing subsequence and fills up the array M and P. M and P must
   point to array of the of length N. They are defined as
   follows. Each item i is processed in order (0<=i<N) and after
   processing item i:

   M[j] — stores the position k of the smallest value X[k] such that
   there is an increasing subsequence of length j+1 ending at X[k] on
   the range k ≤ i.

   P[k] — stores the position of the predecessor of X[k] in the
   longest increasing subsequence ending at X[k].

   M_type should be an output random iterator of size_t. P_type needs
   to support the subscript operator ([]) for output (for example an
   output random iterator or a map) for size_t. Both should be able to
   contain up to the number of elements in the X input range.
 */
template<typename InputIterator, typename M_type, typename P_type>
size_t compute_M_P(const InputIterator X, const InputIterator Xend, M_type M, P_type& P) {
  if(X == Xend)
    return 0;
  const size_t N = std::distance(X, Xend);

  struct comp_type {
    const InputIterator& X_;
    comp_type(const InputIterator& x) : X_(x) { }
    bool operator()(size_t x, size_t y) const {
      return this->X_[x] < this->X_[y];
    }
  };

  const comp_type comp(X);
  size_t          L = 1;
  M[0]              = 0;

  for(size_t i = 1; i < N; ++i) {
    auto   it     = std::lower_bound(M, M + L, i, comp);
    bool   extend = it == M + L;
    if(extend || X[i] < X[*it]) {
      if(extend)
        ++L;
      *it = i;
      if(it > M)
        P[i] = *(it - 1);
    }
  }
  return L;
}

/* Given the X array and the P, M array of indices computed by compute_M_P,
   output the longest increasing subsequence in out (an output forward
   iterator). The sequence is written largest to smallest element. out
   must contain enough space to contain the lis (no check).
 */
template<typename InputIterator, typename OutputIterator, typename M_type, typename P_type>
void sequence_reversed(const InputIterator X, const size_t L, M_type M, const P_type& P,
                       OutputIterator out) {
  size_t j = M[L - 1];
  for(size_t i = 0; i < L; ++i, ++out, j = P[j])
    *out = X[j];
}

/* Similar to sequence_reversed but output the indices within the X range
 */
template<typename OutputIterator, typename M_type, typename P_type>
void indices_reversed(const size_t L, M_type M, const P_type& P,
                       OutputIterator out) {
  size_t j = M[L - 1];
  for(size_t i = 0; i < L; ++i, ++out, j = P[j])
    *out = j;
}



/** Given a sequence given with iterator [X, Xend), return the length
    of the longest increasing subsequence.
*/
template<typename InputIterator>
size_t length(const InputIterator X, const InputIterator Xend) {
  size_t N = std::distance(X, Xend);
  std::vector<size_t> M(N), P(N);
  return compute_M_P(X, Xend, M.begin(), P);
}

/** Given a sequence given with iterator [X, Xend), return a vector of
    elements in the range X of the longest increasing subsequence.
 */
template<typename InputIterator, typename T = typename std::iterator_traits<InputIterator>::value_type>
std::vector<T> sequence(const InputIterator X, const InputIterator Xend) {
  size_t N = std::distance(X, Xend);
  std::vector<size_t> M(N), P(N);
  size_t L = compute_M_P(X, Xend, M.begin(), P);

  std::vector<T> result(L);
  sequence_reversed(X, L, M.cbegin(), P, result.rbegin());
  return result;
}

/** Given a sequence given with iterator [X, Xend), return a vector of
    indices in that range corresponding to the longest increasing
    subsequence.
 */
template<typename InputIterator>
std::vector<size_t> indices(const InputIterator X, const InputIterator Xend) {
  size_t N = std::distance(X, Xend);
  std::vector<size_t> M(N), P(N);
  size_t L = compute_M_P(X, Xend, M.begin(), P);

  std::vector<size_t> result(L);
  indices_reversed(L, M.cbegin(), P, result.rbegin());
  return result;
}



}
#endif /* __LIS_HPP__ */

