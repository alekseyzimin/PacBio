#ifndef __MER_SA_H__
#define __MER_SA_H__

#include "mer_sa_imp.hpp"

namespace mer_sa {
/**
 * Constructs the partial suffix array of a given string.
 * @param T[0..n-1] The input string.
 * @param SA[0..n-mer_size+1] The output array of suffixes.
 * @param n The length of the given string.
 * @param offsets[0..(2**mer_size)] Offsets of mers in suffix array
 * @param mer_size Length of mer used for pre-sorting
 * @param max_size Maximum pattern length for query (i.e. length of suffix sorting)
 * @return total number of mers in T if no error occurred, a negative number otherwise
 */
template<typename CHARPTR, typename SAIDPTR, typename SAIDX = typename compactsufsort::type_traits<SAIDPTR>::SAIDX>
SAIDX
create(CHARPTR T, SAIDPTR SA, SAIDX n, uint64_t* mer_counts, unsigned int mer_size, unsigned int max_size) {
  return mer_sa_imp::SA<CHARPTR, SAIDPTR>::create(T, SA, n, mer_counts, mer_size, max_size);
}

/**
 * Constructs the partial suffix array of a given string in a multi-threaded way.
 * @param T[0..n-1] The input string.
 * @param SA[0..n-mer_size+1] The output array of suffixes.
 * @param n The length of the given string.
 * @param t Number of threads
 * @param offsets[0..(2**mer_size)] Offsets of mers in suffix array
 * @param mer_size Length of mer used for pre-sorting
 * @param max_size Maximum pattern length for query (i.e. length of suffix sorting)
 * @return total number of mers in T if no error occurred, a negative number otherwise
 */
template<typename CHARPTR, typename SAIDPTR, typename SAIDX = typename compactsufsort::type_traits<SAIDPTR>::SAIDX>
SAIDX
create_mt(CHARPTR T, SAIDPTR SA, SAIDX n, unsigned int t, uint64_t* mer_counts, unsigned int mer_size, unsigned int max_size) {
  return mer_sa_imp::SA<CHARPTR, SAIDPTR>::create_mt(T, SA, n, t, mer_counts, mer_size, max_size);
}


/**
 * Search for the pattern P in the string T.
 * @param T[0..Tsize-1] The input string.
 * @param Tsize The length of the given string.
 * @param SA[0..Tsize-1] The input suffix array.
 * @param SAsize The length of the suffix array.
 * @param P[0..Psize-1] The input pattern string.
 * @param Psize The length of the given pattern string.
 * @return The count of matches and output index if no error occurred. -1 otherwise.
 */
template<typename CHARPTR, typename CSAIDPTR, typename CHARPTR2,
         typename SAIDX = typename compactsufsort::type_traits<CSAIDPTR>::SAIDX>
std::pair<SAIDX, SAIDX>
search(const CHARPTR T, SAIDX Tsize, const CSAIDPTR SA, SAIDX SAsize,
       const uint64_t* const counts, unsigned int mer_size, unsigned int max_size,
       const CHARPTR2 P, SAIDX Psize) {
  return mer_sa_imp::SA<CHARPTR, CSAIDPTR>::search(T, Tsize, SA, SAsize, counts, mer_size, max_size, P, Psize);
}

/**
 * Check the consistency of the suffix array.
 * @param T[0..Tsize-1] The input string.
 * @param Tsize The length of the given string.
 * @param SA[0..SAsize-1] The suffix array.
 * @param SAsize The length of the suffix array.
 * @param counts[0..Tsize-mer_size] mer_size-mer counts
 * @param mer_size The minimum length of a pattern to query
 * @param max_size The maximum length of a pattern to query
 * @return True if suffix array is correct.
 */
template<typename CHARPTR, typename CSAIDPTR,
         typename SAIDX = typename compactsufsort::type_traits<CSAIDPTR>::SAIDX>
bool
check(const CHARPTR T, SAIDX Tsize, const CSAIDPTR SA, SAIDX TAsize,
      const uint64_t* const counts, unsigned int mer_size, unsigned int max_size) {
  return mer_sa_imp::SA<CHARPTR, CSAIDPTR>::check(T, Tsize, SA, TAsize, counts, mer_size, max_size);
}

/**
 * Check that every suffixes in the text is in the suffix array.
 * @param T[0..Tsize-1] The input string.
 * @param Tsize The length of the given string.
 * @param SA[0..SAsize-1] The suffix array.
 * @param SAsize The length of the suffix array.
 * @param counts[0..Tsize-mer_size] mer_size-mer counts
 * @param mer_size The minimum length of a pattern to query
 * @param max_size The maximum length of a pattern to query
 * @return True if all string of length max_size from T are in SA at the correct
 *         position.
 */
template<typename CHARPTR, typename CSAIDPTR,
         typename SAIDX = typename compactsufsort::type_traits<CSAIDPTR>::SAIDX>
bool
check_suffixes(const CHARPTR T, SAIDX Tsize, const CSAIDPTR SA, SAIDX TAsize,
              const uint64_t* const counts, unsigned int mer_size, unsigned int max_size) {
  return mer_sa_imp::SA<CHARPTR, CSAIDPTR>::check(T, Tsize, SA, TAsize, counts, mer_size, max_size);
}
} // namespace mer_sa

#endif /* __MER_SA_H__ */
