#ifndef _PB_ALIGNER_HPP_
#define _PB_ALIGNER_HPP_

#include <unordered_map>

#include <src_jf_aligner/jf_aligner.hpp>
#include <src_jf_aligner/super_read_name.hpp>
#include <src_lis/lis_align.hpp>
#include <jellyfish/thread_exec.hpp>
#include <MurmurHash3.h>

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
  const bool               forward_;
  const bool               compress_;
  const bool               duplicated_;

  Multiplexer*             details_multiplexer_;
  Multiplexer*             coords_multiplexer_;
  const std::vector<int>*  unitigs_lengths_; // Lengths of unitigs
  unsigned int             k_len_; // k-mer length used for creating k-unitigs
  int                      max_mer_count_; // max mer count to be used for alignment
  bool                     compact_format_;


  typedef const mer_pos_hash_type::mapped_type list_type;
  typedef const mer_pos_hash_type::position_type position_type;


public:
  // For each super reads, lists, in order the apparition of the
  // k-mers in the PacBio read, the offsets in the super read. The LIS
  // (longest increasing subsequence) on these offsets The longest
  // such subsequence is stored.
  typedef std::pair<int, int> pb_sr_offsets; // first = pb_offset, second = sr_offset

  struct off_lis {
    std::vector<pb_sr_offsets>   offsets;
    std::vector<unsigned int>    lis;
  };
  struct mer_lists {
    off_lis fwd;
    off_lis bwd;
    const frag_lists::frag_info* frag;
  };
  typedef std::map<const char*, mer_lists> frags_pos_type;

  align_pb(int nb_threads, const mer_pos_hash_type& ary, stream_manager& streams,
           double stretch_constant, double stretch_factor, int consecutive, int nmers,
           bool forward = false, bool compress = false, bool duplicated = false) :
    ary_(ary),
    parser_(4 * nb_threads, 100, 1, streams),
    stretch_constant_(stretch_constant),
    stretch_factor_(stretch_factor),
    consecutive_(consecutive),
    nmers_(nmers),
    forward_(forward),
    compress_(compress),
    duplicated_(duplicated),
    details_multiplexer_(0),
    coords_multiplexer_(0),
    unitigs_lengths_(0), k_len_(0),
    max_mer_count_(0), compact_format_(false)
  { }

  align_pb& details_multiplexer(Multiplexer* m) { details_multiplexer_ = m; return *this; }
  align_pb& coords_multiplexer(Multiplexer* m, bool header) {
    coords_multiplexer_ = m;
    if(header) {
      Multiplexer::ostream o(coords_multiplexer_); // Write header
      o << "Rstart Rend Qstart Qend Nmers Rcons Qcons Rcover Qcover Rlen Qlen";
      if(!compact_format_)
        o << " Rname";
      o << " Qname\n";
      o.end_record();
    }
    return *this;
  }
  align_pb& unitigs_lengths(const std::vector<int>* m, unsigned int k_len) {
    if(!forward_) throw std::logic_error("Forward flag must be used if passing unitigs lengths");
    unitigs_lengths_ = m;
    k_len_           = k_len;
    return *this;
  }

  align_pb& max_mer_count(int m) { max_mer_count_ = m; return *this; }
  align_pb& compact_format(bool c) { compact_format_ = c; return *this; }

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
        process_read(ary_, parser, frags_pos, L, stretch_constant_, stretch_factor_, max_mer_count_);
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
        std::distance(ml.fwd.lis.cbegin(), ml.fwd.lis.cend()) > std::distance(ml.bwd.lis.cbegin(), ml.bwd.lis.cend());
      auto                       lisit           = fwd_align ? ml.fwd.lis.cbegin() : ml.bwd.lis.cbegin();
      const auto                 lisend          = fwd_align ? ml.fwd.lis.cend() : ml.bwd.lis.cend();
      auto                       fwd_offit       = ml.fwd.offsets.cbegin();
      const auto                 fwd_offbegin    = fwd_offit;
      const auto                 fwd_offend      = ml.fwd.offsets.cend();
      auto                       bwd_offit       = ml.bwd.offsets.cbegin();
      const auto                 bwd_offbegin    = bwd_offit;
      const auto                 bwd_offend      = ml.bwd.offsets.cend();

      while(fwd_offit != fwd_offend || bwd_offit != bwd_offend) {
        std::pair<int, int> pos;
        bool                part_of_lis = false;
        if(fwd_offit != fwd_offend && (bwd_offit == bwd_offend || fwd_offit->first <= bwd_offit->first)) {
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

  // Contains the summary information of an alignment of a pac-bio and a super read (named qname).
  struct coords_info {
    int              rs, re;
    int              qs, qe;
    int              nb_mers;
    unsigned int     pb_cons, sr_cons;
    unsigned int     pb_cover, sr_cover;
    size_t           ql;
    bool             rn;
    uint64_t         hash;
    std::string      qname;
    std::vector<int> kmers_info; // Number of k-mers in k-unitigs and common between unitigs
    coords_info(const std::string& name, size_t l, int n) :
      nb_mers(n),
      pb_cons(0), sr_cons(0), pb_cover(mer_dna::k()), sr_cover(mer_dna::k()),
      ql(l), rn(false), qname(name)
    { }

    bool operator<(const coords_info& rhs) const {
      return rs < rhs.rs ||
                  (rs == rhs.rs && (re < rhs.re ||
                                    (re == rhs.re && (hash < rhs.hash ||
                                                      (hash == rhs.hash && ql < rhs.ql)))));
    }
  };

  // Helper class that computes the number of k-mers in each k-unitigs
  // and that are shared between the k-unitigs. It handles errors
  // gracefully (the kmers_info array is then empty) and can be a NOOP
  // if the super read name is empty.
  template<typename AlignerT>
  struct compute_kmers_info {
    std::vector<int>&                info_;
    std::unique_ptr<super_read_name> sr_name_;
    unsigned int                     cunitig_;
    unsigned int                     cend_;
    const AlignerT&                  aligner_;

    compute_kmers_info(std::vector<int>& info, const std::string& name, const AlignerT& aligner) :
      info_(info), sr_name_(aligner.k_len_ && name.size() ? new super_read_name(name) : 0),
      cunitig_(0), cend_(0), aligner_(aligner)
    {
      if(sr_name_) {
        const auto unitig_id = sr_name_->unitig_id(0);
        if(unitig_id != super_read_name::invalid && unitig_id < aligner_.unitigs_lengths_->size()) {
          info_.resize(2 * sr_name_->size() - 1, 0);
          cend_ = unitigs_lengths(unitig_id);
        } else { // error
          sr_name_.reset();
        }
      }
    }

    int unitigs_lengths(long int id) const { return (*aligner_.unitigs_lengths_)[id]; }
    size_t nb_unitigs() const { return aligner_.unitigs_lengths_->size(); }
    unsigned int k_len() const { return aligner_.k_len_; }

    void add_mer(const int pos) {
      if(!sr_name_) return;

      unsigned int       cendi;
      const unsigned int sr_pos = abs(pos);
      while(sr_pos + mer_dna::k() > cend_ + 1) {
        const auto unitig_id = sr_name_->unitig_id(++cunitig_);
        if(unitig_id != super_read_name::invalid && unitig_id < nb_unitigs())
          cend_ += unitigs_lengths(unitig_id) - k_len();
        else
          goto error;
      }
      ++info_[2 * cunitig_];
      cendi = cend_;
      for(unsigned int i = cunitig_; (i < sr_name_->size() - 1) && (sr_pos + k_len() > cendi); ++i) {
        ++info_[2 * i + 1];
        ++info_[2 * i + 2];
        const auto unitig_id = sr_name_->unitig_id(i + 1);
        if(unitig_id != super_read_name::invalid && unitig_id < nb_unitigs())
          cendi += unitigs_lengths(unitig_id) - k_len();
        else
          goto error;
      }
      return;

    error:
      sr_name_.reset();
      info_.resize(0);
    }
  };


  // Compute the statistics of the matches in frags_pos (all the
  // matches to a given pac-bio read)
  std::vector<coords_info> compute_coordinates(const frags_pos_type& frags_pos) {
    std::vector<coords_info> coords;
    for(auto it = frags_pos.cbegin(); it != frags_pos.cend(); ++it) {
      const align_pb::mer_lists& ml          = it->second;
      const auto                 fwd_nb_mers = std::distance(ml.fwd.lis.begin(), ml.fwd.lis.end());
      const auto                 bwd_nb_mers = std::distance(ml.bwd.lis.begin(), ml.bwd.lis.end());
      const bool                 fwd_align   = fwd_nb_mers >= bwd_nb_mers;
      const auto                 nb_mers     = fwd_align ? fwd_nb_mers : bwd_nb_mers;
      if(nb_mers < nmers_) continue; // Enough matching mers?

      // Compute consecutive mers and covered bases. Compute hash of
      // positions in pb read of aligned mers. If k_len_ > 0, also
      // compute the number of k-mers aligned in each k-unitigs of the
      // super-read. If an error occurs (unknown k-unitigs, k-unitigs
      // too short, etc.), an empty vector is returned.
      MurmurHash3A                      hasher;
      coords_info                       info(forward_ && !fwd_align ? reverse_super_read_name(ml.frag->name) : ml.frag->name,
                                             ml.frag->len, nb_mers);
      const std::vector<pb_sr_offsets>& offsets   = fwd_align ? ml.fwd.offsets : ml.bwd.offsets;
      const std::vector<unsigned int>&  lis       = fwd_align ? ml.fwd.lis : ml.bwd.lis;
      compute_kmers_info<align_pb>      kmers_info(info.kmers_info, info.qname, *this);
      {
        auto                              lisit   = lis.cbegin();
        pb_sr_offsets                     prev    = offsets[*lisit];
        pb_sr_offsets                     cur;
        kmers_info.add_mer(fwd_align ? prev.second : info.ql + prev.second - mer_dna::k() + 2);
        for(++lisit; lisit != lis.cend(); prev = cur, ++lisit) {
          cur                         = offsets[*lisit];
          const unsigned int pb_diff  = cur.first - prev.first;
          info.pb_cons               += pb_diff == 1;
          info.pb_cover              += std::min(mer_dna::k(), pb_diff);
          const unsigned int sr_diff  = cur.second - prev.second;
          info.sr_cons               += sr_diff == 1;
          info.sr_cover              += std::min(mer_dna::k(), sr_diff);
          hasher.add((uint32_t)cur.first);
          kmers_info.add_mer(fwd_align ? cur.second : info.ql + cur.second - mer_dna::k() + 2);
        }
      }
      if(info.pb_cons < (unsigned int)consecutive_ && info.sr_cons < (unsigned int)consecutive_) continue;

      uint64_t r1, r2;
      hasher.finalize(r1, r2);
      info.hash  = r1 ^ r2;
      auto first = offsets[lis.front()];
      auto last  = offsets[lis.back()];
      info.rs    = first.first;
      info.re    = last.first + mer_dna::k() - 1;

      info.qs = first.second;
      info.qe = last.second;
      if(info.qs < 0) {
        info.qs = -info.qs + mer_dna::k() - 1;
        info.qe = -info.qe;
        if(forward_) {
          info.qs      = info.ql - info.qs + 1;
          info.qe      = info.ql - info.qe + 1;
          info.rn      = true;
        }
      } else {
        info.qe += mer_dna::k() - 1;
      }

      coords.push_back(info);
    }

    return coords;
  }

  void print_coords(Multiplexer::ostream& out, const std::string& pb_name, size_t pb_size, const frags_pos_type& frags_pos) {
    std::vector<coords_info> coords = compute_coordinates(frags_pos);
    auto nb_lines = std::distance(coords.cbegin(), coords.cend());
    if(nb_lines == 0) return;

    std::sort(coords.begin(), coords.end());
    if(compact_format_)
      out << ">" << nb_lines << " " << pb_name << "\n";
    for(auto it = coords.cbegin(), pit = it; it != coords.cend(); pit = it, ++it) {
      if(duplicated_ && it != pit && it->rs == pit->rs && it->re == pit->re && it->hash == pit->hash)
        continue;
      out << it->rs << " " << it->re << " " << it->qs << " " << it->qe << " "
          << it->nb_mers << " "
          << it->pb_cons << " " << it->sr_cons << " "
          << it->pb_cover << " " << it->sr_cover << " "
          << pb_size << " " << it->ql;
      if(!compact_format_)
        out << " " << pb_name;
      out << " " << it->qname;
      for(auto mit : it->kmers_info)
        out << " " << mit;
      out << "\n";
    }
    out.end_record();
  }

  // For each k-unitigs in the super read qname, output its length, the number of k-mers
  // void print_mers_in_unitigs(Multiplexer::ostream& out, const mer_lists* ml, const std::string qname) {
  //   // super_read_name unitigs(qname);
  // }

  static void process_read(const mer_pos_hash_type& ary, parse_sequence& parser,
                           frags_pos_type& frags_pos, lis_align::forward_list<lis_align::element<double> >& L,
                           double a, double b, const int max_mer_count = 0) {
    while(parser.next()) { // Process each k-mer
      const bool is_canonical = parser.m < parser.rm;
      auto it = ary.find_pos(is_canonical ? parser.m : parser.rm);
      if(max_mer_count && std::distance(it, ary.pos_end()) > max_mer_count) continue;
      for( ; it != ary.pos_end(); ++it) {
        mer_lists& ml = frags_pos[it->frag->name];
        ml.frag       = it->frag;
        const int offset = is_canonical ? it->offset : -it->offset;
        if(offset > 0)
          ml.fwd.offsets.push_back(pb_sr_offsets(parser.offset + 1, offset));
        else
          ml.bwd.offsets.push_back(pb_sr_offsets(parser.offset + 1, offset));
      }
    }
    if(frags_pos.empty()) return;

    // Compute LIS forward and backward on every super reads.
    for(auto it = frags_pos.begin(); it != frags_pos.end(); ++it) {
      mer_lists& mer_list = it->second;
      mer_list.fwd.lis.clear();
      mer_list.bwd.lis.clear();
      L.clear();
      lis_align::indices(mer_list.fwd.offsets.cbegin(), mer_list.fwd.offsets.cend(),
                         L, mer_list.fwd.lis, a, b);
      L.clear();
      lis_align::indices(mer_list.bwd.offsets.cbegin(), mer_list.bwd.offsets.cend(),
                         L, mer_list.bwd.lis, a, b);
    }
  }

  static void process_read(const mer_pos_hash_type& ary, parse_sequence& parser,
                           frags_pos_type& frags_pos, double a, double b) {
    lis_align::forward_list<lis_align::element<double> > L;
    process_read(ary, parser, frags_pos, L, a, b);
  }

  // Reverse the name of a super read. For example 1R_2F_3F becomes
  // 3R_2R_1F. If the name is not valid, name is returned unchanged.
  static std::string reverse_super_read_name(const std::string& name) {
    std::string res;
    const size_t nsize = name.size();
    if(nsize == 0)
      return res;

    size_t ppos = nsize;
    size_t pos  = name.find_last_of('_');
    while(true) {
      pos = (pos == std::string::npos) ? 0 : pos + 1;
      if(pos > ppos - 2) goto invalid_name;
      res += name.substr(pos, ppos - pos - 1);
      switch(name[ppos - 1]) {
      case 'R': res += 'F'; break;
      case 'F': res += 'R'; break;
      default: goto invalid_name;
      }
      if(pos == 0)
        break;
      res += '_';
      ppos = pos - 1;
      pos  = name.find_last_of('_', ppos - 1);
    }
    return res;

  invalid_name:
    res = name;
    return res;
  }

};

void align_pb_reads(int threads, mer_pos_hash_type& hash, double stretch_const, double stretch_factor,
                    int consecutive, int nmers, bool compress, bool forward, bool duplicated,
                    const char* pb_path,
                    const char* coords_path = 0, const char* details_path = 0,
                    const int* lengths = 0, const int nb_unitigs = 0, const unsigned int k_len = 0) {
  output_file details, coords;
  if(coords_path) coords.open(coords_path, threads);
  if(details_path) details.open(details_path, threads);
  file_vector files;
  files.push_back(pb_path);
  stream_manager streams(files.cbegin(), files.cend());
  align_pb aligner(threads, hash, streams, stretch_const, stretch_factor,
                   consecutive, nmers, forward, compress, duplicated);
  if(coords_path) aligner.coords_multiplexer(coords.multiplexer(), true);
  if(details_path) aligner.details_multiplexer(details.multiplexer());
  std::vector<int> unitigs_lengths;
  if(nb_unitigs && lengths && k_len) {
    unitigs_lengths.insert(unitigs_lengths.begin(), lengths, lengths + nb_unitigs);
    aligner.unitigs_lengths(&unitigs_lengths, k_len);
  }
  aligner.exec_join(threads);
}

#endif /* _PB_ALIGNER_HPP_ */
