#ifndef _PB_ALIGNER_HPP_
#define _PB_ALIGNER_HPP_

#include <src_jf_aligner/jf_aligner.hpp>
#include <src_lis/lis_align.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jflib/multiplexed_io.hpp>

class align_pb : public jellyfish::thread_exec {
  const mer_pos_hash_type& ary_;
  read_parser              parser_;
  jflib::o_multiplexer&    out_multiplexer_;
  int                      stretch_constant_, stretch_factor_;

  typedef const mer_pos_hash_type::mapped_type list_type;
  typedef const mer_pos_hash_type::position_type position_type;


public:
  // For each super reads, lists, in order the apparition of the
  // k-mers in the PacBio read, the offsets in the super read. The LIS
  // (longest increasing subsequence) on these offsets The longest
  // such subsequence is stored.
  typedef std::pair<int, int> pb_sr_offsets; // first = pb_offset, second = sr_offset

  struct mer_lists {
    std::vector<pb_sr_offsets> offsets;
    std::vector<unsigned int>  lis;
    bool                       rev;
  };
  typedef std::map<const char*, mer_lists> frags_pos_type;

  align_pb(int nb_threads, const mer_pos_hash_type& ary, stream_manager& streams, jflib::o_multiplexer& out,
           int stretch_constant, int stretch_factor) :
    ary_(ary),
    parser_(4 * nb_threads, 100, 1, streams),
    out_multiplexer_(out),
    stretch_constant_(stretch_constant),
    stretch_factor_(stretch_factor)
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
        process_read(ary_, parser, frags_pos, stretch_constant_, stretch_factor_);
        out << ">" << job->data[i].header << "\n";
        for(auto it = frags_pos.cbegin(); it != frags_pos.cend(); ++it) {
          out << "+" << it->first << "\n";
          const align_pb::mer_lists& ml      = it->second;
          auto                       lisit   = ml.lis.cbegin();
          auto                       lisend  = ml.lis.cend();
          size_t                     current = 0;
          for(auto offit = ml.offsets.cbegin(); offit != ml.offsets.cend(); ++offit, ++current) {
            bool part_of_lis = (lisit < lisend && *lisit == current);
            out << " " <<(part_of_lis ? "[" : "")
                << offit->first << ":" << offit->second
                << (part_of_lis ? "]" : "");
            if(part_of_lis)
              ++lisit;
          }
          out << "\n";
        }
        out << jflib::endr;
      }
    }
  }

  static void process_read(const mer_pos_hash_type& ary, parse_sequence& parser, frags_pos_type& frags_pos,
                           int a, int b) {
    while(parser.next()) { // Process each k-mer
      const bool is_canonical = parser.m < parser.rm;
      list_type* list = ary.find_pos(is_canonical ? parser.m : parser.rm);
      if(!list) // mer not found in superreads
        continue;
      for(auto it = list->cbegin(); it != list->cend(); ++it) {
        frags_pos[it->frag].offsets.push_back(pb_sr_offsets(parser.offset + 1, is_canonical ? it->offset : -it->offset));
      }
    }

    // Compute LIS forward and backward on every super reads.
    for(auto it = frags_pos.begin(); it != frags_pos.end(); ++it) {
      mer_lists& mer_list = it->second;
      //      mer_list.lis = lis::indices(mer_list.offsets.cbegin(), mer_list.offsets.cend());
      mer_list.lis = lis_align::indices(mer_list.offsets.cbegin(), mer_list.offsets.cend(),
                                        a, b);
    }
  }
};


#endif /* _PB_ALIGNER_HPP_ */
