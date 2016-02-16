/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef _JF_ALIGNER_HPP_
#define _JF_ALIGNER_HPP_

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <multiplexer.hpp>
#include <src_jf_aligner/frag_info.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <src_jf_aligner/output_file.hpp>

namespace err = jellyfish::err;
using jellyfish::mer_dna;
typedef mer_dna long_mer_type;
typedef jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1> short_mer_type;
typedef std::vector<const char*> file_vector;
typedef jellyfish::stream_manager<file_vector::const_iterator> stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> read_parser;

// Parse a DNA sequence and break the input into k-mers. The offset is
// 1-based
namespace parse_sequence_ns {
template<typename MT>
struct mer_info {
  typedef MT mer_type;
  mer_type     m, rm;
  unsigned int len;
  bool         valid;

  void reset() {
    len    = 0;
  }

  bool append(char base) {
    int code = mer_type::code(base);
    if(mer_type::not_dna(code)) {
      len = 0;
      return false;
    }
    ++len;
    m.shift_left(code);
    rm.shift_right(mer_type::complement(code));
    valid = (len >= mer_type::k());
    return valid;
  }

  bool is_canonical() const { return m < rm; }
};

template<typename TupleType, size_t I = std::tuple_size<TupleType>::value>
struct mers_function {
  // Reset all the mers
  static void reset(TupleType& t) {
    std::get<I-1>(t).reset();
    mers_function<TupleType, I-1>::reset(t);
  }

  // Append a new character to all the mers. Return true if at least
  // one of the mers is valid (i.e. one of the append method returned
  // true).
  static bool append(TupleType& t, const char c, const bool res = false) {
    const bool success = std::get<I-1>(t).append(c);
    return mers_function<TupleType, I-1>::append(t, c, res || success);
  }
};
template<typename TupleType>
struct mers_function<TupleType, 0> {
  static void reset(const TupleType& t) { }
  static bool append(const TupleType&t, const char c, const bool res = false) { return res; }
};

template<typename... MTs>
struct parser_base {
  typedef std::tuple<mer_info<MTs>...> mer_list_type;
  const bool                  compress;
  mer_list_type               mers;
  unsigned int                seq_offset;
  char                        prev_base;
  std::string::const_iterator base;
  std::string::const_iterator end;

  parser_base(const bool c = false) : compress(c) { }
  parser_base(const std::string& s, const bool c = false) : compress(c) { reset(s); }
  parser_base(const std::string* s, const bool c = false) : compress(c) { reset(s); }

  template<size_t I>
  typename std::tuple_element<I, mer_list_type>::type& mer() { return std::get<I>(mers); }
  template<size_t I>
  typename std::tuple_element<I, mer_list_type>::type const& mer() const { return std::get<I>(mers); }

  template<size_t I>
  int offset() const { return seq_offset - std::tuple_element<I, mer_list_type>::type::mer_type::k() + 1; }

  size_t size() const { return end - base; }


  void reset(const std::string& s) {
    seq_offset = 0;
    prev_base  = '\0';
    base       = s.begin();
    end        = s.end();
    mers_function<mer_list_type>::reset(mers);
  }
  void reset(const std::string* s) { reset(*s); }

  bool next() {
    while(base < end) {
      char new_base = *base;
      ++base; ++seq_offset;

      if(compress && prev_base == new_base) continue;
      prev_base = new_base;
      if(mers_function<mer_list_type>::append(mers, new_base)) return true;
    }
    return false;
  }
};
} // namespace parse_sequence_ns

typedef parse_sequence_ns::parser_base<mer_dna> parse_sequence;
typedef parse_sequence_ns::parser_base<short_mer_type> short_parse_sequence;
typedef parse_sequence_ns::parser_base<long_mer_type, short_mer_type> parse_sequence2;

#endif /* _JF_ALIGNER_HPP_ */
