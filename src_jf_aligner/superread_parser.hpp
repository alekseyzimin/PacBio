#ifndef _SUPERREAD_PARSER_HPP_
#define _SUPERREAD_PARSER_HPP_

#include <vector>

#include <jellyfish/thread_exec.hpp>
#include <src_jf_aligner/jf_aligner.hpp>

class superreads_read_mers : public jellyfish::thread_exec {
  mer_pos_hash_type*       ary_long_;
  short_mer_pos_hash_type* ary_short_;
  read_parser              parser_;
  frag_lists&              names_; // super reads names
  const bool               compress_;

public:
  superreads_read_mers(int nb_threads, mer_pos_hash_type* ary_long, short_mer_pos_hash_type* ary_short,
                       frag_lists& names, stream_manager& streams,
                       bool compress = false) :
    ary_long_(ary_long),
    ary_short_(ary_short),
    parser_(4 * nb_threads, 100, 1, streams),
    names_(names),
    compress_(compress)
  { }

  virtual void start(int thid);
  void parse_one_mer(int thid);
  void parse_two_mers(int thid);
};

void superread_parse(int threads, mer_pos_hash_type& hash, frag_lists& names,
                     file_vector::const_iterator begin, file_vector::const_iterator end,
                     bool compress = false);

void superread_parse(int threads, mer_pos_hash_type& hash, frag_lists& names, const char* file, bool compress = false);

#endif /* _SUPERREAD_PARSER_HPP_ */
