#ifndef _PB_ALIGNER_HPP_
#define _PB_ALIGNER_HPP_

#include <src_jf_aligner/jf_aligner.hpp>
#include <src_lis/lis.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jflib/multiplexed_io.hpp>

class align_pb : public jellyfish::thread_exec {
  const mer_pos_hash_type& ary_;
  read_parser              parser_;
  jflib::o_multiplexer&    out_multiplexer_;

  typedef const mer_pos_hash_type::mapped_type list_type;
  typedef const mer_pos_hash_type::position_type position_type;


public:
  // For each super reads, lists, in order the apparition of the
  // k-mers in the PacBio read, the offsets in the super read. The
  // LIS (longest increasing subsequence) on these offsets is then
  // compute, both forward and backward. The longest such
  // subsequence is stored.
  struct mer_lists {
    std::vector<int> offsets;
    std::vector<int> lis;
    bool             rev;
  };
  typedef std::map<const char*, mer_lists> frags_pos_type;

  align_pb(int nb_threads, const mer_pos_hash_type& ary, stream_manager& streams, jflib::o_multiplexer& out) :
    ary_(ary),
    parser_(4 * nb_threads, 100, 1, streams),
    out_multiplexer_(out)
  { }

  virtual void start(int thid) {
    mer_dna         tmp_m;
    parse_sequence  parser;
    jflib::omstream out(out_multiplexer_);
    frags_pos_type  frags_pos;

    while(true) {
      read_parser::job job(parser_);
      if(job.is_empty()) break;

      for(size_t i = 0; i < job->nb_filled; ++i) { // Process each read
        parser.reset(job->data[i].seq);
        frags_pos.clear();
        process_read(ary_, parser, frags_pos);
      }
    }
  }

  static void process_read(const mer_pos_hash_type& ary, parse_sequence& parser, frags_pos_type& frags_pos) {
    while(parser.next()) { // Process each k-mer
      const bool is_canonical = parser.m < parser.rm;
      list_type* list = ary.find_pos(is_canonical ? parser.m : parser.rm);
      if(!list) // mer not found in superreads
        continue;
      for(auto it = list->cbegin(); it != list->cend(); ++it) {
        frags_pos[it->frag].offsets.push_back(is_canonical ? it->offset : -it->offset);
      }
    }

    // Compute LIS forward and backward on every super reads.
    for(auto it = frags_pos.begin(); it != frags_pos.end(); ++it) {
      mer_lists& mer_list = it->second;
      std::vector<int> fwd_indices = lis::sequence(mer_list.offsets.cbegin(), mer_list.offsets.cend());
      std::vector<int> bwd_indices = lis::sequence(mer_list.offsets.crbegin(), mer_list.offsets.crend());
      mer_list.rev = bwd_indices.size() > fwd_indices.size();
      if(mer_list.rev)
        mer_list.lis = std::move(bwd_indices);
      else
        mer_list.lis = std::move(fwd_indices);
    }
  }
};


#endif /* _PB_ALIGNER_HPP_ */
