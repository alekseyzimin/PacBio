#ifndef _SUPERREAD_PARSER_HPP_
#define _SUPERREAD_PARSER_HPP_

#include <src_jf_aligner/mer_pos_hash.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>


using jellyfish::mer_dna;
typedef mer_pos_hash<mer_dna> mer_pos_hash_type;

typedef std::vector<const char*> file_vector;
typedef jellyfish::stream_manager<file_vector::const_iterator> stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> read_parser;
typedef std::vector<std::vector<std::string> > name_lists;

class superreads_read_mers : public jellyfish::thread_exec {
  mer_pos_hash_type& ary_;
  read_parser        parser_;
  name_lists& names_; // super reads names

public:
  superreads_read_mers(int nb_threads, mer_pos_hash_type& ary, name_lists& names, stream_manager& streams) :
    ary_(ary),
    parser_(4 * nb_threads, 100, 1, streams),
    names_(names)
  { }

  virtual void start(int thid) {
    mer_dna m, rm;

    while(true) {
      read_parser::job job(parser_);
      if(job.is_empty()) break;

      for(size_t i = 0; i < job->nb_filled; ++i) { // Process each read
        names_[thid].push_back(job->data[i].header);
        const std::string& name = names_[i].back();
        std::string& seq = job->data[i].seq;

        auto         base   = seq.begin();
        unsigned int len    = 0; // Length of low quality stretch
        int offset = -mer_dna::k() + 1;
        for( ; base != seq.end(); ++base, ++offset) {
          int code = mer_dna::code(*base);
          if(mer_dna::not_dna(code)) {
            len = 0;
            continue;
          }
          m.shift_left(code);
          rm.shift_right(mer_dna::complement(code));
          ++len;
          if(len >= mer_dna::k()) {
            const bool is_canonical = m < rm;
            ary_.push_front(is_canonical ? m : rm, name, is_canonical ? offset : -offset);
          }
        }
      }
    }
  }
};

void superread_parse(int threads, mer_pos_hash_type& hash, name_lists& names,
                     file_vector::const_iterator begin, file_vector::const_iterator end) {
  if(names.size() < (size_t)threads)
    names.resize(threads);
  stream_manager streams(begin, end);
  superreads_read_mers reader(threads, hash, names, streams);
  reader.exec_join(threads);
}

void superread_parse(int threads, mer_pos_hash_type& hash, name_lists& names, const char* file) {
  file_vector files;
  files.push_back(file);
  superread_parse(threads, hash, names, files.cbegin(), files.cend());
}

#endif /* _SUPERREAD_PARSER_HPP_ */
