#ifndef _SUPERREAD_PARSER_HPP_
#define _SUPERREAD_PARSER_HPP_

#include <src_jf_aligner/jf_aligner.hpp>

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
    parse_sequence parser;

    while(true) {
      read_parser::job job(parser_);
      if(job.is_empty()) break;

      for(size_t i = 0; i < job->nb_filled; ++i) { // Process each read
        const char* header = names_.push_back(thid, job->data[i].header);
        parser.reset(job->data[i].seq);

        while(parser.next()) { // Process each k-mer
            const bool is_canonical = parser.m < parser.rm;
            ary_.push_front(is_canonical ? parser.m : parser.rm,
                            header,
                            is_canonical ? parser.offset : -parser.offset);
        }
      }
    }
  }
};

void superread_parse(int threads, mer_pos_hash_type& hash, name_lists& names,
                     file_vector::const_iterator begin, file_vector::const_iterator end) {
  names.ensure(threads);
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
