#ifndef _PB_ALIGNER_HPP_
#define _PB_ALIGNER_HPP_

#include <src_jf_aligner/jf_aligner.hpp>
#include <src_lis/lis_align.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jflib/multiplexed_io.hpp>

using jflib::o_multiplexer;
using jflib::omstream;

/**
 * A unique ptr to a multiplexed output stream with a convenient
 * constructor.
 */
class mstream : public std::unique_ptr<omstream> {
public:
  mstream(o_multiplexer* m) :
    std::unique_ptr<omstream>(m ? new omstream(*m) : 0)
  { }
};

// struct frags_pos_type {
//   std::map
// };

class align_pb : public jellyfish::thread_exec {
  const mer_pos_hash_type& ary_;
  read_parser              parser_;
  int                      stretch_constant_, stretch_factor_;
  int                      consecutive_, nmers_;
  const bool               compress_;

  o_multiplexer* details_multiplexer_;
  o_multiplexer* coords_multiplexer_;


  typedef const mer_pos_hash_type::mapped_type list_type;
  typedef const mer_pos_hash_type::position_type position_type;


public:
  // For each super reads, lists, in order the apparition of the
  // k-mers in the PacBio read, the offsets in the super read. The LIS
  // (longest increasing subsequence) on these offsets The longest
  // such subsequence is stored.
  typedef std::pair<int, int> pb_sr_offsets; // first = pb_offset, second = sr_offset

  struct mer_lists {
    std::vector<pb_sr_offsets>   offsets;
    std::vector<unsigned int>    lis;
    const frag_lists::frag_info* frag;
  };
  typedef std::map<const char*, mer_lists> frags_pos_type;

  align_pb(int nb_threads, const mer_pos_hash_type& ary, stream_manager& streams,
           int stretch_constant, int stretch_factor, int consecutive, int nmers,
           bool compress = false) :
    ary_(ary),
    parser_(4 * nb_threads, 100, 1, streams),
    stretch_constant_(stretch_constant),
    stretch_factor_(stretch_factor),
    consecutive_(consecutive),
    nmers_(nmers),
    compress_(compress),
    details_multiplexer_(0),
    coords_multiplexer_(0)
  { }

  align_pb& details_multiplexer(o_multiplexer* m) { details_multiplexer_ = m; return *this; }
  align_pb& coords_multiplexer(o_multiplexer* m) {
    coords_multiplexer_ = m;
    omstream o(*coords_multiplexer_); // Write header
    o << "Rstart Rend Qstart Qend Nmers Qcons Rcons Qcover Rcover Rlen Qlen Qname Rname\n";
    o << jflib::endr;
    return *this;
  }

  virtual void start(int thid) {
    mer_dna         tmp_m;
    parse_sequence  parser(compress_);
    frags_pos_type  frags_pos;

    mstream     details(details_multiplexer_);
    mstream     coords(coords_multiplexer_);
    std::string name;

    while(true) {
      read_parser::job job(parser_);
      if(job.is_empty()) break;

      for(size_t i = 0; i < job->nb_filled; ++i) { // Process each read
        auto name_end = job->data[i].header.find_first_of(" \t\n\v\f\r");
        name = job->data[i].header.substr(0, name_end);
        parser.reset(job->data[i].seq);
        frags_pos.clear();
        process_read(ary_, parser, frags_pos, stretch_constant_, stretch_factor_);
        if(details) print_details(*details, name, frags_pos);
        if(coords) print_coords(*coords, name, job->data[i].seq.size(), frags_pos);
      }
    }
  }

  void print_details(omstream& out, const std::string& pb_name, const frags_pos_type& frags_pos) {
    bool first = true;
    for(auto it = frags_pos.cbegin(); it != frags_pos.cend(); ++it, first = false) {
      if(first)
        out << pb_name;
      else
        out << ".";
      out << " " << it->first;
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

  void print_coords(omstream& out, const std::string& pb_name, size_t pb_size, const frags_pos_type& frags_pos) {
    for(auto it = frags_pos.cbegin(); it != frags_pos.cend(); ++it) {
      const align_pb::mer_lists& ml    = it->second;
      const auto nb_mers               = std::distance(ml.lis.begin(), ml.lis.end());
      if(nb_mers < nmers_) continue; // Enough matching mers

      // Compute consecutive mers and covered bases
      unsigned int pb_cover = mer_dna::k(), sr_cover = mer_dna::k(), pb_cons = 0, sr_cons = 0;
      {
        auto lisit = ml.lis.cbegin();
        pb_sr_offsets prev = ml.offsets[*lisit];
        pb_sr_offsets cur;
        for(++lisit; lisit != ml.lis.cend(); prev = cur, ++lisit) {
          cur                   = ml.offsets[*lisit];
          unsigned int pb_diff  = cur.first - prev.first;
          pb_cons              += pb_diff == 1;
          pb_cover             += std::min(mer_dna::k(), pb_diff);
          unsigned int sr_diff  = cur.second - prev.second;
          sr_cons              += sr_diff == 1;
          sr_cover             += std::min(mer_dna::k(), sr_diff);
        }
      }
      if(pb_cons < (unsigned int)consecutive_ && sr_cons < (unsigned int)consecutive_) continue;

      auto                       first = ml.offsets[ml.lis.front()];
      auto                       last  = ml.offsets[ml.lis.back()];
      auto                       s2    = first.second;
      auto                       e2    = last.second;
      if(s2 < 0)
        s2 -= mer_dna::k() - 1;
      else
        e2 += mer_dna::k() - 1;
      out << first.first << " " << (last.first + mer_dna::k() - 1) // RS RE
          << " " << s2 << " " << e2 // QS QE
          << " " << nb_mers // N
          << " " << pb_cons << " " << sr_cons
          << " " << pb_cover << " " << sr_cover
          << " " << pb_size << " " << ml.frag->len // RL QL
          << " " << pb_name << " " << ml.frag->name << "\n"; // Q R
    }
    out << jflib::endr;
  }

  static void process_read(const mer_pos_hash_type& ary, parse_sequence& parser, frags_pos_type& frags_pos,
                           int a, int b) {
    while(parser.next()) { // Process each k-mer
      const bool is_canonical = parser.m < parser.rm;
      auto it = ary.find_pos(is_canonical ? parser.m : parser.rm);
      for( ; it != ary.pos_end(); ++it) {
        mer_lists& ml = frags_pos[it->frag->name];
        ml.frag       = it->frag;
        ml.offsets.push_back(pb_sr_offsets(parser.offset + 1, is_canonical ? it->offset : -it->offset));
      }
    }
    if(frags_pos.empty()) return;

    // Compute LIS forward and backward on every super reads.
    for(auto it = frags_pos.begin(); it != frags_pos.end(); ++it) {
      mer_lists& mer_list = it->second;
      mer_list.lis = lis_align::indices(mer_list.offsets.cbegin(), mer_list.offsets.cend(),
                                        a, b);
    }
  }
};


#endif /* _PB_ALIGNER_HPP_ */
