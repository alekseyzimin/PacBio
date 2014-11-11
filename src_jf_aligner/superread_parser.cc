#include <tuple>
#include <src_jf_aligner/superread_parser.hpp>

void superreads_read_mers::start(int thid) {
  if(ary_short_ == 0) {
    parse_one_mer(thid);
  } else {
    parse_two_mers(thid);
  }
}

void superreads_read_mers::parse_one_mer(int thid) {
  parse_sequence            parser(compress_);
  mer_pos_hash_type::thread at_long(*ary_long_);

  while(true) {
    read_parser::job job(parser_);
    if(job.is_empty()) break;

    for(size_t i = 0; i < job->nb_filled; ++i) { // Process each read
      auto name_end = job->data[i].header.find_first_of(" \t\n\v\f\r");
      auto header   = names_.push_back(thid, job->data[i].seq.length(), job->data[i].header.substr(0, name_end));
      parser.reset(job->data[i].seq);

      while(parser.next()) { // Process each k-mer
        if(!parser.mer<0>().m.is_homopolymer()) {
          const bool is_canonical = parser.mer<0>().is_canonical();
          at_long.push_front(is_canonical ? parser.mer<0>().m : parser.mer<0>().rm,
                             header,
                             is_canonical ? parser.offset<0>() : -parser.offset<0>());
        }
      }
    }
  }
}

void superreads_read_mers::parse_two_mers(int thid) {
  parse_sequence2                 parser(compress_);
  mer_pos_hash_type::thread       at_long(*ary_long_);
  short_mer_pos_hash_type::thread at_short(*ary_short_);

  while(true) {
    read_parser::job job(parser_);
    if(job.is_empty()) break;

    for(size_t i = 0; i < job->nb_filled; ++i) { // Process each read
      auto name_end = job->data[i].header.find_first_of(" \t\n\v\f\r");
      auto header   = names_.push_back(thid, job->data[i].seq.length(), job->data[i].header.substr(0, name_end));
      parser.reset(job->data[i].seq);

      while(parser.next()) { // Process each k-mer
        if(parser.mer<0>().valid && !parser.mer<0>().m.is_homopolymer()) {
          const bool is_canonical = parser.mer<0>().is_canonical();
          at_long.push_front(is_canonical ? parser.mer<0>().m : parser.mer<0>().rm,
                             header,
                             is_canonical ? parser.offset<0>() : -parser.offset<0>());
        }
        if(parser.mer<0>().valid && !parser.mer<1>().m.is_homopolymer()) {
          const bool is_canonical = parser.mer<1>().is_canonical();
          at_short.push_front(is_canonical ? parser.mer<1>().m : parser.mer<1>().rm,
                              header,
                              is_canonical ? parser.offset<1>() : -parser.offset<1>());
        }
      }
    }
  }
}

void superread_parse(int threads, mer_pos_hash_type& hash, frag_lists& names,
                     file_vector::const_iterator begin, file_vector::const_iterator end,
                     bool compress) {
  names.ensure(threads);
  stream_manager streams(begin, end);
  superreads_read_mers reader(threads, &hash, 0, names, streams, compress);
  reader.exec_join(threads);
}

void superread_parse(int threads, mer_pos_hash_type& hash, frag_lists& names, const char* file, bool compress) {
  file_vector files;
  files.push_back(file);
  superread_parse(threads, hash, names, files.cbegin(), files.cend(), compress);
}
