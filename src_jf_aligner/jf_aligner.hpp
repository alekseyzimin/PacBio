#ifndef _JF_ALIGNER_HPP_
#define _JF_ALIGNER_HPP_

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <src_jf_aligner/frag_info.hpp>
#include <src_jf_aligner/mer_pos_hash.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>

using jellyfish::mer_dna;
typedef mer_pos_hash<frag_lists::frag_info, mer_dna> mer_pos_hash_type;
typedef std::vector<const char*> file_vector;
typedef jellyfish::stream_manager<file_vector::const_iterator> stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> read_parser;

// Parse a DNA sequence and break the input into k-mers. The offset is
// 1-based
struct parse_sequence {
  const bool                  compress;
  mer_dna                     m, rm;
  int                         offset;
  int                         prev;
  unsigned int                len;
  std::string::const_iterator base;
  std::string::const_iterator end;

  parse_sequence(const bool c = false) : compress(c) { }
  parse_sequence(const std::string& s, const bool c = false) : compress(c) { reset(s); }
  parse_sequence(const std::string* s, const bool c = false) : compress(c) { reset(s); }

  void reset(const std::string& s) {
    offset = -mer_dna::k() + 1;
    len    = 0;
    base   = s.begin();
    end    = s.end();
    prev   = -1;
  }
  void reset(const std::string* s) { reset(*s); }

  bool next() {
    while(base < end) {
      int code = mer_dna::code(*base);
      ++base; ++offset;
      int oprev = prev;
      prev      = code;
      if(mer_dna::not_dna(code)) {
        len = 0;
        continue;
      }
      if(compress && code == oprev) continue;
      ++len;
      m.shift_left(code);
      rm.shift_right(mer_dna::complement(code));
      if(len >= mer_dna::k())
        return true;
    }
    return false;
  }
};

#endif /* _JF_ALIGNER_HPP_ */
