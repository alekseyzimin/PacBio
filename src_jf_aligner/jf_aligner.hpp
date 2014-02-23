#ifndef _JF_ALIGNER_HPP_
#define _JF_ALIGNER_HPP_

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <src_jf_aligner/mer_pos_hash.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>

using jellyfish::mer_dna;
typedef mer_pos_hash<mer_dna> mer_pos_hash_type;
// typedef std::vector<std::vector<const char*> > name_lists;
typedef std::vector<const char*> file_vector;
typedef jellyfish::stream_manager<file_vector::const_iterator> stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> read_parser;

// Multiple vector to hold names of fragments. There is a vector for
// each thread and the string is copied when push_back is called.
class name_lists {
  std::vector<std::vector<const char*> > names_;
public:
  name_lists(size_t threads) : names_(threads) { }
  ~name_lists() {
    for(auto it = names_.cbegin(); it != names_.cend(); ++it)
      for(auto it2 = it->cbegin(); it2 != it->cend(); ++it2)
        free((void*)*it2);
  }

  size_t size() const { return names_.size(); }

  void ensure(size_t threads) {
    if(names_.size() < threads)
      names_.resize(threads);
  }

  const char* push_back(int thid, const char* s) {
    const char* res = strdup(s);
    names_[thid].push_back(res);
    return res;
  }

  const char* push_back(int thid, const std::string& s) {
    return push_back(thid, s.c_str());
  }

  const std::vector<const char*> operator[](int i) const {
    return names_[i];
  }
};

// Parse a DNA sequence and break the input into k-mers. The offset is
// 1-based
struct parse_sequence {
  mer_dna                     m, rm;
  int                         offset;
  unsigned int                len;
  std::string::const_iterator base;
  std::string::const_iterator end;

  parse_sequence() = default;
  parse_sequence(const std::string& s) { reset(s); }
  parse_sequence(const std::string* s) { reset(s); }

  void reset(const std::string& s) {
    offset = -mer_dna::k() + 1;
    len    = 0;
    base   = s.begin();
    end    = s.end();
  }
  void reset(const std::string* s) { reset(*s); }

  bool next() {
    while(base < end) {
      int code = mer_dna::code(*base);
      ++base; ++offset;
      if(mer_dna::not_dna(code)) {
        len = 0;
        continue;
      }
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
