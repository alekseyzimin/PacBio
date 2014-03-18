#ifndef _PB_ALIGNER_HPP_
#define _PB_ALIGNER_HPP_

#include <src_jf_aligner/jf_aligner.hpp>
#include <src_lis/lis_align.hpp>
#include <jellyfish/thread_exec.hpp>

/**
 * A unique ptr to a multiplexed output stream with a convenient
 * constructor.
 */
class mstream : public std::unique_ptr<Multiplexer::ostream> {
public:
  mstream(Multiplexer* m) :
    std::unique_ptr<Multiplexer::ostream>(m ? new Multiplexer::ostream(m) : 0)
  { }
};

class align_pb : public jellyfish::thread_exec {
  const mer_pos_hash_type& ary_;
  read_parser              parser_;
  double                   stretch_constant_, stretch_factor_;
  int                      consecutive_, nmers_;
  const bool               compress_;

  Multiplexer* details_multiplexer_;
  Multiplexer* coords_multiplexer_;


  typedef const mer_pos_hash_type::mapped_type list_type;
  typedef const mer_pos_hash_type::position_type position_type;


public:
  // For each super reads, lists, in order the apparition of the
  // k-mers in the PacBio read, the offsets in the super read. The LIS
  // (longest increasing subsequence) on these offsets The longest
  // such subsequence is stored.
  typedef std::pair<int, int> pb_sr_offsets; // first = pb_offset, second = sr_offset

  struct mer_lists {
    std::vector<pb_sr_offsets>   fwd_offsets;
    std::vector<pb_sr_offsets>   bwd_offsets;
    std::vector<unsigned int>    fwd_lis;
    std::vector<unsigned int>    bwd_lis;
    const frag_lists::frag_info* frag;
  };
  typedef std::map<const char*, mer_lists> frags_pos_type;

  align_pb(int nb_threads, const mer_pos_hash_type& ary, stream_manager& streams,
           double stretch_constant, double stretch_factor, int consecutive, int nmers,
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

  align_pb& details_multiplexer(Multiplexer* m) { details_multiplexer_ = m; return *this; }
  align_pb& coords_multiplexer(Multiplexer* m, bool header) {
    coords_multiplexer_ = m;
    if(header) {
      Multiplexer::ostream o(coords_multiplexer_); // Write header
      o << "Rstart Rend Qstart Qend Nmers Rcons Qcons Rcover Qcover Rlen Qlen Rname Qname\n";
      o.end_record();
    }
    return *this;
  }

  virtual void start(int thid) {
    mer_dna         tmp_m;
    parse_sequence  parser(compress_);
    frags_pos_type  frags_pos;
    lis_align::forward_list<lis_align::element<double> > L; // L buffer

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
        process_read(ary_, parser, frags_pos, L, stretch_constant_, stretch_factor_);
        if(details) print_details(*details, name, frags_pos);
        if(coords) print_coords(*coords, name, job->data[i].seq.size(), frags_pos);
      }
    }
  }

  void print_details(Multiplexer::ostream& out, const std::string& pb_name, const frags_pos_type& frags_pos) {
    for(auto it = frags_pos.cbegin(); it != frags_pos.cend(); ++it) {
      out << pb_name << " " << it->first;
      const align_pb::mer_lists& ml              = it->second;
      const bool                       fwd_align =
        std::distance(ml.fwd_lis.cbegin(), ml.fwd_lis.cend()) > std::distance(ml.bwd_lis.cbegin(), ml.bwd_lis.cend());
      auto                       lisit           = fwd_align ? ml.fwd_lis.cbegin() : ml.bwd_lis.cbegin();
      const auto                 lisend          = fwd_align ? ml.fwd_lis.cend() : ml.bwd_lis.cend();
      auto                       fwd_offit       = ml.fwd_offsets.cbegin();
      const auto                 fwd_offbegin    = fwd_offit;
      const auto                 fwd_offend      = ml.fwd_offsets.cend();
      auto                       bwd_offit       = ml.bwd_offsets.cbegin();
      const auto                 bwd_offbegin    = bwd_offit;
      const auto                 bwd_offend      = ml.bwd_offsets.cend();

      while(fwd_offit != fwd_offend || bwd_offit != bwd_offend) {
        std::pair<int, int> pos;
        bool                part_of_lis = false;
        if(fwd_offit != fwd_offend && (bwd_offit == bwd_offend || fwd_offit->first < bwd_offit->first)) {
          pos = *fwd_offit;
          part_of_lis = fwd_align && (lisit < lisend) && (*lisit == std::distance(fwd_offbegin, fwd_offit));
          ++fwd_offit;
        } else if(bwd_offit != bwd_offend && (fwd_offit == fwd_offend || bwd_offit->first < fwd_offit->first)) {
          pos = *bwd_offit;
          part_of_lis = !fwd_align && (lisit < lisend) && (*lisit == std::distance(bwd_offbegin, bwd_offit));
          ++bwd_offit;
        }
        out << " " <<(part_of_lis ? "[" : "")
            << pos.first << ":" << pos.second
            << (part_of_lis ? "]" : "");
        if(part_of_lis)
          ++lisit;
      }
      out << "\n";
    }
    out.end_record();
  }

  void print_coords(Multiplexer::ostream& out, const std::string& pb_name, size_t pb_size, const frags_pos_type& frags_pos) {
    struct coords_info {
      int          rs, re;
      int          qs, qe;
      int          nb_mers;
      unsigned int pb_cons, sr_cons;
      unsigned int pb_cover, sr_cover;
      size_t       ql;
      const char*  qname;
      coords_info(const char* name, size_t l, int n) :
        nb_mers(n),
        pb_cons(0), sr_cons(0), pb_cover(mer_dna::k()), sr_cover(mer_dna::k()),
        ql(l), qname(name)
      { }
      bool operator<(const coords_info& rhs) const { return rs < rhs.rs || (rs == rhs.rs && re < rhs.re); }
    };
    std::vector<coords_info> coords;

    for(auto it = frags_pos.cbegin(); it != frags_pos.cend(); ++it) {
      const align_pb::mer_lists& ml          = it->second;
      const auto                 fwd_nb_mers = std::distance(ml.fwd_lis.begin(), ml.fwd_lis.end());
      const auto                 bwd_nb_mers = std::distance(ml.bwd_lis.begin(), ml.bwd_lis.end());
      const bool                 fwd_align   = fwd_nb_mers >= bwd_nb_mers;
      const auto                 nb_mers     = fwd_align ? fwd_nb_mers : bwd_nb_mers;
      if(nb_mers < nmers_) continue; // Enough matching mers

      // Compute consecutive mers and covered bases
      coords_info info(ml.frag->name, ml.frag->len, nb_mers);
      const std::vector<pb_sr_offsets>& offsets = fwd_align ? ml.fwd_offsets : ml.bwd_offsets;
      const std::vector<unsigned int>&  lis     = fwd_align ? ml.fwd_lis : ml.bwd_lis;
      {
        auto                              lisit   = lis.cbegin();
        pb_sr_offsets                     prev    = offsets[*lisit];
        pb_sr_offsets                     cur;
        for(++lisit; lisit != lis.cend(); prev = cur, ++lisit) {
          cur                         = offsets[*lisit];
          const unsigned int pb_diff  = cur.first - prev.first;
          info.pb_cons               += pb_diff == 1;
          info.pb_cover              += std::min(mer_dna::k(), pb_diff);
          const unsigned int sr_diff  = cur.second - prev.second;
          info.sr_cons               += sr_diff == 1;
          info.sr_cover              += std::min(mer_dna::k(), sr_diff);
        }
      }
      if(info.pb_cons < (unsigned int)consecutive_ && info.sr_cons < (unsigned int)consecutive_) continue;

      auto                       first = offsets[lis.front()];
      auto                       last  = offsets[lis.back()];
      info.rs = first.first;
      info.re = last.first + mer_dna::k() - 1;

      info.qs = first.second;
      info.qe = last.second;
      if(info.qs < 0) {
        info.qs = -info.qs + mer_dna::k() - 1;
        info.qe = -info.qe;
      } else {
        info.qe += mer_dna::k() - 1;
      }

      coords.push_back(info);
    }

    std::sort(coords.begin(), coords.end());
    for(auto it = coords.cbegin(); it != coords.cend(); ++it) {
      out << it->rs << " " << it->re << " " << it->qs << " " << it->qe << " "
          << it->nb_mers << " "
          << it->pb_cons << " " << it->sr_cons << " "
          << it->pb_cover << " " << it->sr_cover << " "
          << pb_size << " " << it->ql << " "
          << pb_name << " " << it->qname << "\n";
    }
    out.end_record();
  }

  static void process_read(const mer_pos_hash_type& ary, parse_sequence& parser,
                           frags_pos_type& frags_pos, lis_align::forward_list<lis_align::element<double> >& L,
                           double a, double b) {
    while(parser.next()) { // Process each k-mer
      if(parser.m.is_homopolymer()) continue;
      const bool is_canonical = parser.m < parser.rm;
      auto it = ary.find_pos(is_canonical ? parser.m : parser.rm);
      for( ; it != ary.pos_end(); ++it) {
        mer_lists& ml = frags_pos[it->frag->name];
        ml.frag       = it->frag;
        const int offset = is_canonical ? it->offset : -it->offset;
        if(offset > 0)
          ml.fwd_offsets.push_back(pb_sr_offsets(parser.offset + 1, offset));
        else
          ml.bwd_offsets.push_back(pb_sr_offsets(parser.offset + 1, offset));
      }
    }
    if(frags_pos.empty()) return;

    // Compute LIS forward and backward on every super reads.
    for(auto it = frags_pos.begin(); it != frags_pos.end(); ++it) {
      mer_lists& mer_list = it->second;
      mer_list.fwd_lis.clear();
      mer_list.bwd_lis.clear();
      L.clear();
      lis_align::indices(mer_list.fwd_offsets.cbegin(), mer_list.fwd_offsets.cend(),
                         L, mer_list.fwd_lis, a, b);
      L.clear();
      lis_align::indices(mer_list.bwd_offsets.cbegin(), mer_list.bwd_offsets.cend(),
                         L, mer_list.bwd_lis, a, b);
    }
  }

  static void process_read(const mer_pos_hash_type& ary, parse_sequence& parser,
                           frags_pos_type& frags_pos, double a, double b) {
    lis_align::forward_list<lis_align::element<double> > L;
    process_read(ary, parser, frags_pos, L, a, b);
  }
};

void align_pb_reads(int threads, mer_pos_hash_type& hash, double stretch_const, double stretch_factor,
                    int consecutive, int nmers, bool compress,
                    const char* pb_path,
                    const char* coords_path = 0, const char* details_path = 0) {
  output_file details, coords;
  if(coords_path) coords.open(coords_path, threads);
  if(details_path) details.open(details_path, threads);
  file_vector files;
  files.push_back(pb_path);
  stream_manager streams(files.cbegin(), files.cend());
  align_pb aligner(threads, hash, streams, stretch_const, stretch_factor,
                   consecutive, nmers, compress);
  if(coords_path) aligner.coords_multiplexer(coords.multiplexer(), true);
  if(details_path) aligner.details_multiplexer(details.multiplexer());
  aligner.exec_join(threads);
}

#endif /* _PB_ALIGNER_HPP_ */
