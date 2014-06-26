#ifndef _PB_ALIGNER_HPP_
#define _PB_ALIGNER_HPP_

#include <cmath>
#include <unordered_map>
#include <limits>

#include <src_jf_aligner/jf_aligner.hpp>
#include <src_jf_aligner/super_read_name.hpp>
#include <src_jf_aligner/least_square_2d.hpp>
#include <src_lis/lis_align.hpp>
#include <jellyfish/thread_exec.hpp>

class align_pb {
  const mer_pos_hash_type& ary_;
  double                   stretch_constant_, stretch_factor_; // Maximum stretch in LIS
  const bool               forward_;
  int                      max_mer_count_; // max mer count to be used for alignment
  double                   matching_mers_factor_;
  double                   matching_bases_factor_;

  const std::vector<int>*  unitigs_lengths_; // Lengths of unitigs
  unsigned int             k_len_; // k-mer length used for creating k-unitigs


  typedef const mer_pos_hash_type::mapped_type list_type;
  typedef const mer_pos_hash_type::position_type position_type;


public:
  // For each super reads, lists, in order the apparition of the
  // k-mers in the PacBio read, the offsets in the super read. The LIS
  // (longest increasing subsequence) on these offsets The longest
  // such subsequence is stored. The offsets are 1-based
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
  typedef lis_align::forward_list<lis_align::element<double>> lis_buffer_type;

  align_pb(const mer_pos_hash_type& ary,
           double stretch_constant, double stretch_factor,
           bool forward = false, int max_mer_count = 0,
           double matching_mers = 0.0, double matching_bases = 0.0) :
    ary_(ary),
    stretch_constant_(stretch_constant),
    stretch_factor_(stretch_factor),
    forward_(forward),
    max_mer_count_(max_mer_count),
    matching_mers_factor_(matching_mers),
    matching_bases_factor_(matching_bases),
    unitigs_lengths_(0), k_len_(0)
  { }

  align_pb& unitigs_lengths(const std::vector<int>* m, unsigned int k_len) {
    if(!forward_) throw std::logic_error("Forward flag must be used if passing unitigs lengths");
    unitigs_lengths_ = m;
    k_len_           = k_len;
    return *this;
  }

  align_pb& max_mer_count(int m) { max_mer_count_ = m; return *this; }

  // Contains the summary information of an alignment of a pac-bio and a super read (named qname).
  struct coords_info {
    int              rs, re;
    int              qs, qe;
    int              nb_mers;
    unsigned int     pb_cons, sr_cons;
    unsigned int     pb_cover, sr_cover;
    size_t           rl, ql;
    bool             rn;
    std::string      qname;
    std::vector<int> kmers_info; // Number of k-mers in k-unitigs and common between unitigs
    std::vector<int> bases_info; // Number of bases in k-unitigs and common between unitigs
    double           stretch, offset, avg_err; // Least square stretch, offset and average error
    coords_info() = default;
    coords_info(const std::string& name, size_t rl, size_t ql, int n) :
      nb_mers(n),
      pb_cons(0), sr_cons(0), pb_cover(mer_dna::k()), sr_cover(mer_dna::k()),
      rl(rl), ql(ql), rn(false), qname(name),
      stretch(0), offset(0), avg_err(0)
    { }
    coords_info(int rs_, int re_, int qs_, int qe_, int nb_mers_,
                unsigned int pb_cons_, unsigned sr_cons_,
                unsigned int pb_cover_, unsigned int sr_cover_,
                size_t rl_, size_t ql_, bool rn_,
                double stretch_, double offset_, double avg_err_,
                std::string qname_) :
      rs(rs_), re(re_), qs(qs_), qe(qe_), nb_mers(nb_mers_),
      pb_cons(pb_cons_), sr_cons(sr_cons_),
      pb_cover(pb_cover_), sr_cover(sr_cover_),
      rl(rl_), ql(ql_), rn(rn_), qname(qname_),
      stretch(stretch_), offset(offset_), avg_err(avg_err_)
    { }

    bool operator<(const coords_info& rhs) const {
      return rs < rhs.rs || (rs == rhs.rs &&
                             (re < rhs.re || (re == rhs.re && ql < rhs.ql)));
    }

    // Transform raw coordinates obtained from matching to final
    // coordinates, encompassing the first and last bases matched (and
    // not the first base of the kmer mathing), properly reverse the
    // match if needed and forward is set to true.
    void canonicalize(bool forward) {
      if(qs < 0) {
        if(forward) {
          qs      = ql + qs - mer_dna::k() + 2;
          qe      = ql + qe + 1;
          rn      = true;
          offset -= stretch * (ql + 1) - mer_dna::k();
        } else {
          qs       = -qs + mer_dna::k() - 1;
          qe       = -qe;
          stretch  = -stretch;
          offset  += mer_dna::k() - 1;
        }
      } else {
        qe += mer_dna::k() - 1;
      }
    }

    double imp_s() const { return std::max(1.0, std::min((double)rl, stretch + offset)); }
    double imp_e() const { return std::max(1.0, std::min((double)rl, stretch * ql + offset)); }
    int imp_len() const { return abs(lrint(imp_e() - imp_s())) + 1; }
    bool min_bases(double factor) { return factor * (imp_len() - 2 * (int)mer_dna::k()) <= pb_cover; }

    bool min_mers(double factor) { return factor * (imp_len() - mer_dna::k() + 1) <= nb_mers; }
  };

  typedef std::vector<coords_info> coords_info_type;

  // Helper class that computes the number of k-mers in each k-unitigs
  // and that are shared between the k-unitigs. It handles errors
  // gracefully (the kmers_info array is then empty) and can be a NOOP
  // if the super read name is empty.
  template<typename AlignerT>
  struct compute_kmers_info {
    std::vector<int>&                mers_;
    std::vector<int>&                bases_;
    std::unique_ptr<super_read_name> sr_name_;
    unsigned int                     cunitig_;
    int                              cend_;
    int                              prev_pos_;
    const AlignerT&                  aligner_;

    compute_kmers_info(std::vector<int>& mers, std::vector<int>& bases, const std::string& name, const AlignerT& aligner) :
      mers_(mers), bases_(bases), sr_name_(aligner.k_len_ && name.size() ? new super_read_name(name) : 0),
      cunitig_(0), cend_(0), prev_pos_(-mer_dna::k()), aligner_(aligner)
    {
      if(sr_name_) {
        const auto unitig_id = sr_name_->unitig_id(0);
        if(unitig_id != super_read_name::invalid && unitig_id < nb_unitigs()) {
          mers_.resize(2 * sr_name_->size() - 1, 0);
          bases_.resize(2 * sr_name_->size() - 1, 0);
          cend_ = unitigs_lengths(unitig_id);
        } else // error
          sr_name_.reset();
      }
    }


    int unitigs_lengths(long int id) const { return (*aligner_.unitigs_lengths_)[id]; }
    size_t nb_unitigs() const { return aligner_.unitigs_lengths_->size(); }
    unsigned int k_len() const { return aligner_.k_len_; }

    void add_mer(const int pos) {
      if(!sr_name_) return;

      int       cendi;
      const int sr_pos = abs(pos);
      const int new_bases   = std::min((int)mer_dna::k(), sr_pos - prev_pos_);
      while(sr_pos + (int)mer_dna::k() > cend_ + 1) {
        if(cend_ >= sr_pos) {
          const int nb_bases = cend_ - std::max(sr_pos, prev_pos_ + (int)mer_dna::k()) + 1;
          bases_[2 * cunitig_]     += nb_bases;
          bases_[2 * cunitig_ + 1] += nb_bases;
        }
        const auto unitig_id = sr_name_->unitig_id(++cunitig_);
        if(unitig_id == super_read_name::invalid || unitig_id >= nb_unitigs())
          goto error;
        cend_ += unitigs_lengths(unitig_id) - k_len() + 1;
      }
      ++mers_[2 * cunitig_];
      bases_[2 * cunitig_] += new_bases;
      cendi                 = cend_;
      for(unsigned int i = cunitig_; (i < sr_name_->size() - 1) && (sr_pos + mer_dna::k() > cendi - k_len() + 1); ++i) {
        const int  full_mer   = sr_pos + (int)k_len() > cendi + 1;
        mers_[2 * i + 1]     += full_mer;
        mers_[2 * i + 2]     += full_mer;
        const int  nb_bases   = std::min(new_bases, sr_pos + (int)mer_dna::k() - cendi + (int)k_len() - 2);
        bases_[2 * i + 1]    += nb_bases;
        bases_[2 * i + 2]    += nb_bases;
        const auto unitig_id  = sr_name_->unitig_id(i + 1);
        if(unitig_id != super_read_name::invalid && unitig_id < nb_unitigs())
          cendi += unitigs_lengths(unitig_id) - k_len() + 1;
        else
          goto error;
      }
      prev_pos_ = sr_pos;
      return;

    error:
      sr_name_.reset();
      mers_.clear();
      bases_.clear();
    }
  };


  // Compute the statistics of the matches in frags_pos (all the
  // matches to a given pac-bio read)
  coords_info compute_coords_info(const mer_lists& ml, const size_t pb_size) const;
  void compute_coords(const frags_pos_type& frags_pos, const size_t pb_size, coords_info_type& coords) const;
  coords_info_type compute_coords(const frags_pos_type& frags_pos, const size_t pb_size) const {
    coords_info_type coords;
    compute_coords(frags_pos, pb_size, coords);
    return coords;
  }

  void align_sequence(parse_sequence& parser, const size_t pb_size,
                      coords_info_type& coords, frags_pos_type& frags, lis_buffer_type& L) const;
  std::pair<coords_info_type, frags_pos_type> align_sequence(parse_sequence& parser, const size_t pb_size) const {
    std::pair<coords_info_type, frags_pos_type> res;
    lis_buffer_type                             L;
    align_sequence(parser, pb_size, res.first, res.second, L);
    return res;
  }
  std::pair<coords_info_type, frags_pos_type> align_sequence(const std::string& seq) const {
    parse_sequence parser(seq);
    return align_sequence(parser, seq.size());
  }

  static void print_coords_header(Multiplexer* m, bool compact);
  static void print_coords(Multiplexer::ostream& out, const std::string& pb_name, size_t const pb_size,
                           const bool compact, const coords_info_type& coords);
  static void print_details(Multiplexer::ostream& out, const std::string& pb_name, const frags_pos_type& frags_pos);

  static void fetch_super_reads(const mer_pos_hash_type& ary, parse_sequence& parser,
                                frags_pos_type& frags_pos, const int max_mer_count = 0);

  static void do_LIS(frags_pos_type& frags_pos, lis_buffer_type& L, double a, double b);

  static void do_LIS(frags_pos_type& frags_pos, double a, double b) {
    lis_buffer_type L;
    do_LIS(frags_pos, L, a, b);
  }

  // Reverse the name of a super read. For example 1R_2F_3F becomes
  // 3R_2R_1F. If the name is not valid, name is returned unchanged.
  static std::string reverse_super_read_name(const std::string& name);

  class thread {
    const align_pb&  align_data_;
    frags_pos_type   frags_pos_;
    coords_info_type coords_;
    lis_buffer_type  L_;

  public:
    thread(const align_pb& a) : align_data_(a) { }

    void align_sequence(parse_sequence& parser, const size_t pb_size);
    void align_sequence(const std::string& seq) {
      parse_sequence parser(seq);
      align_sequence(parser, seq.size());
    }
    const coords_info_type& coords() const { return coords_; }
    const frags_pos_type& frags_pos() const { return frags_pos_; }
  };
  friend class thread;

  // static void align_pb_reads(mer_pos_hash_type& hash, double stretch_const, double stretch_factor,
  //                            bool forward, const char* pb_path,
  //                            const int* lengths = 0, const int nb_unitigs = 0, const unsigned int k_len = 0,
  //                            double matching_mers = 0.0, double matching_bases = 0.0);
};


#endif /* _PB_ALIGNER_HPP_ */
