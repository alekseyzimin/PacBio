#include <src_jf_aligner/superread_parser.hpp>

void superreads_read_mers::start(int thid) {
  parse_sequence parser(compress_);

  while(true) {
    read_parser::job job(parser_);
    if(job.is_empty()) break;

    for(size_t i = 0; i < job->nb_filled; ++i) { // Process each read
      auto name_end = job->data[i].header.find_first_of(" \t\n\v\f\r");
      auto header   = names_.push_back(thid, job->data[i].seq.length(), job->data[i].header.substr(0, name_end));
      parser.reset(job->data[i].seq);

      while(parser.next()) { // Process each k-mer
        if(parser.m.is_homopolymer()) continue;
        const bool is_canonical = parser.m < parser.rm;
        ary_.push_front(is_canonical ? parser.m : parser.rm,
                        header,
                        is_canonical ? parser.offset : -parser.offset);
      }
    }
  }
}

void superread_parse(int threads, mer_pos_hash_type& hash, frag_lists& names,
                     file_vector::const_iterator begin, file_vector::const_iterator end,
                     bool compress) {
  names.ensure(threads);
  stream_manager streams(begin, end);
  superreads_read_mers reader(threads, hash, names, streams, compress);
  reader.exec_join(threads);
}

void superread_parse(int threads, mer_pos_hash_type& hash, frag_lists& names, const char* file, bool compress) {
  file_vector files;
  files.push_back(file);
  superread_parse(threads, hash, names, files.cbegin(), files.cend(), compress);
}
