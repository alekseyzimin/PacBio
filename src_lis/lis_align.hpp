#ifndef __LIS2_H__
#define __LIS2_H__

#include <vector>
#include <iterator>
#include <type_traits>
#include <algorithm>
#include <forward_list>

#include <src_jf_aligner/lf_forward_list.hpp>


namespace lis_align {
template<typename T>
class sum_buffer {
  std::vector<T> v_;
  size_t         next_;
  bool           filled_;
  T              sum_;

public:
  sum_buffer(size_t size) : v_(size), next_(0), filled_(false), sum_(0) { }
  //  sum_buffer(const sum_buffer& rhs) : v_(rhs.v_), next_(rhs.next_), filled_(rhs.filled_), sum_(rhs.sum_) { }

  bool filled() const { return filled_; }
  bool will_be_filled() const { return filled_ || next_ == size() - 1; }
  T sum() const { return sum_; }
  size_t size() const { return v_.size(); }
  void push_back(const T& x) {
    if(v_.size()) {
      sum_       = test_sum(x);
      v_[next_]  = x;
      next_      = (next_ + 1) % v_.size();
      filled_    = filled_ || (next_ == 0);
    }
  }
  T test_sum(const T& x) const {
    T res = sum_ + x;
    if(filled_ || next_ > 0) res -= v_[next_];
    return res;
  }
};

template<typename T>
struct sum_pair : public std::pair<T, T> {
  sum_pair(const T x, const T y) : std::pair<T, T>(x, y) { }
  sum_pair(const T x) : std::pair<T, T>(x, x) { }
  sum_pair() : std::pair<T, T>(T(), T()) { }

  sum_pair& operator+=(const sum_pair& rhs) {
    this->first  += rhs.first;
    this->second += rhs.second;
    return *this;
  }
  sum_pair operator+(const sum_pair& rhs) const {
    sum_pair res(*this);
    return res += rhs;
  }

  sum_pair& operator-=(const sum_pair& rhs) {
    this->first  -= rhs.first;
    this->second -= rhs.second;
    return *this;
  }
  sum_pair operator-(const sum_pair& rhs) const {
    sum_pair res(*this);
    return res -= rhs;
  }
};

struct affine_capped {
  const double a, b, C;
  affine_capped(double a_, double b_, double C_) : a(a_), b(b_), C(C_) {
    std::cerr << "afine_capped a:" << a << " b:" << b << " C:" << C << '\n';
  }
  template<typename T>
  bool operator()(const sum_pair<T>& s) const {
    return (s.first <= b + a * s.second) && (s.second <= b + a * s.first) && s.first <= C && s.second <= C;
  }
};

struct linear {
  const double a;
  linear(double a_) : a(a_) { }
  template<typename T>
  bool operator()(const sum_pair<T>& s) const {
    return (s.first <= a * s.second) && (s.second <= a * s.first);
  }
};

/**
 * Compute an alignment on an array X where each element is a pair of
 * offsets. It returns the longest alignment (in term of number of
 * elements) the second offset are all increasing and where the spans
 * in first offsets is bounded by an affine relation to the spans offsets
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
  size_t                   elt;             // Index of element in X
  unsigned int             len;             // Length of LIS
  sum_buffer<sum_pair<T> > span_window; // Span in reference and query in a window
  sum_pair<T>              span_full;      // Span from first to last k-mer
};
template<typename T>
std::ostream& operator<<(std::ostream& os, const element<T>& x) {
  return os << "<" << x.len << ", " << x.span_window.first << ", " << x.span_window.second << ">";
}
template<typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& x) {
  return os << "<" << x.first << ", " << x.second << ">";
}

template<typename InputIterator, typename F1, typename F2,
         typename T=typename std::iterator_traits<InputIterator>::value_type::first_type>
std::pair<unsigned int, unsigned int> compute_L_P(const InputIterator X, const InputIterator Xend,
                                                  std::forward_list<element<T> >& L, std::vector<unsigned int>& P,
                                                  size_t window_size, F1& accept_mer, F2& accept_sequence) {
  unsigned int longest = 0, longest_ind = 0;
  const size_t N       = std::distance(X, Xend);
  std::cout << "---\n";
  for(unsigned int i = 0 ; i < N; ++i) {
    element<T>    e_longest = { i, 1, window_size };
    unsigned int& j_longest = P[i];
    j_longest       = N;

    // Go in decreasing order of subsequence length. Stop at first
    // possible extension, as further subsequence length are too short
    // to improve on this extension.
    auto prev = L.cbefore_begin();
    const auto before_begin = prev;
    auto it   = L.cbegin();
    for( ; it != L.cend() && it->len >= e_longest.len; ++it) {
      const unsigned int j = it->elt;
      if(X[i].second > X[j].second && e_longest.len < it->len + 1) {
        sum_pair<T>       add(X[i].first - X[j].first, X[i].second - X[j].second);
        const sum_pair<T> new_span = it->span_window.test_sum(add);
        std::cout << "test " << X[i] << ',' << X[j] << " len:" << it->len << " span:" << it->span_window.sum()
                  << " new_span:" << new_span
                  << ' ' << it->span_window.will_be_filled() << ' ' << accept_mer(new_span) << '\n';
        if(!it->span_window.will_be_filled() || accept_mer(new_span)) {
          e_longest.len         = it->len + 1;
          j_longest             = j;
          e_longest.span_window = it->span_window;
          e_longest.span_window.push_back(add);
          e_longest.span_full   = it->span_full + add;
          break;
        }
      }
      if(prev == before_begin || it->len < prev->len)
        prev = it;
    }
    L.insert_after(prev, e_longest);
    std::cout << "longest " << X[i] << " span:" << e_longest.span_full << " len:" << e_longest.len << " prev:" << j_longest << X[j_longest] << ' ' << accept_sequence(e_longest.span_full) << '\n';
    if(longest < e_longest.len && accept_sequence(e_longest.span_full)) {
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

template<typename InputIterator, typename F1, typename F2,
         typename T=typename std::iterator_traits<InputIterator>::value_type::first_type>
unsigned int indices(const InputIterator X, const InputIterator Xend,
                     std::forward_list<element<T> >& L, std::vector<unsigned int>& res,
                     size_t window_size, F1& accept_mer, F2& accept_sequence) {
  const size_t N = std::distance(X, Xend);
  std::vector<unsigned int> P(N);
  const std::pair<unsigned int, unsigned int> lis = compute_L_P(X, Xend, L, P, window_size, accept_mer, accept_sequence);

  if(res.size() < lis.first)
    res.resize(lis.first);
  indices_reversed(P, lis.first, lis.second, res.rbegin());
  return lis.first;
}

template<typename InputIterator, typename F1, typename F2,
         typename T=typename std::iterator_traits<InputIterator>::value_type::first_type>
unsigned int indices(const InputIterator X, const InputIterator Xend, std::vector<unsigned int>& res,
                     size_t window_size, F1& accept_mer, F2& accept_sequence) {
  std::forward_list<element<T> > L;
  return indices(X, Xend, L, res, window_size, accept_mer, accept_sequence);
}

template<typename InputIterator, typename F1, typename F2>
std::vector<unsigned int> indices(const InputIterator X, const InputIterator Xend,
                                  size_t window_size, F1& accept_mer, F2& accept_sequence) {
  std::vector<unsigned int> res;

  indices(X, Xend, res, window_size, accept_mer, accept_sequence);
  return res;
}
} // namespace lis2


#endif /* __LIS2_H__ */
