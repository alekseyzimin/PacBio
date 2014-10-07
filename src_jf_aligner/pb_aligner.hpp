#ifndef _PB_ALIGNER_HPP_
#define _PB_ALIGNER_HPP_

#include <cmath>
#include <unordered_map>
#include <limits>
#include <map>

#include <src_jf_aligner/jf_aligner.hpp>
#include <src_jf_aligner/super_read_name.hpp>
#include <src_jf_aligner/least_square_2d.hpp>
#include <src_lis/lis_align.hpp>
#include <jellyfish/thread_exec.hpp>
#include <debug.hpp>

namespace align_pb {

typedef std::forward_list<lis_align::element<double>> lis_buffer_type;

// For each super reads, lists, in order the apparition of the
// k-mers in the PacBio read, the offsets in the super read. The LIS
// (longest increasing subsequence) on these offsets The longest
// such subsequence is stored. The offsets are 1-based
typedef std::pair<int, int> pb_sr_offsets; // first = pb_offset, second = sr_offset

struct off_lis {
  std::vector<pb_sr_offsets>   offsets;
  std::vector<unsigned int>    lis;

  template<typename F1, typename F2>
  void do_LIS(F1& accept_mer, F2& accept_sequence, size_t window_size, lis_buffer_type& L) {
    L.clear();
    lis.clear();
    lis_align::indices(offsets.cbegin(), offsets.cend(), L, lis, window_size, accept_mer, accept_sequence);
  }
  void discard_LIS() {
    auto lis_it = lis.cbegin();
    if(lis.cbegin() == lis.cend()) return;
    auto w_oit = offsets.begin() + *lis_it;
    ++lis_it;
    for(auto r_oit = w_oit + 1; r_oit != offsets.end(); ++r_oit) {
      if(lis_it != lis.cend() && (r_oit - offsets.begin()) == *lis_it) {
        ++lis_it;
      } else {
        *w_oit = *r_oit;
        ++w_oit;
      }
    }
    offsets.resize(offsets.size() - lis.size());
  }

  template<typename F1, typename F2>
  void discard_update_LIS(F1& accept_mer, F2& accept_sequence, size_t window_size, lis_buffer_type& L) {
    discard_LIS();
    do_LIS(accept_mer, accept_sequence, window_size, L);
  }
};

struct mer_lists {
  off_lis fwd;
  off_lis bwd;
  const frag_lists::frag_info* frag;
  template<typename F1, typename F2>
  void do_LIS(F1& accept_mer, F2& accept_sequence, size_t window_size, lis_buffer_type& L) {
    fwd.do_LIS(accept_mer, accept_sequence, window_size, L);
    bwd.do_LIS(accept_mer, accept_sequence, window_size, L);
  }
  template<typename F1, typename F2>
  void discard_update_LIS(F1& accept_mer, F2& accept_sequence, size_t window_size, lis_buffer_type& L) {
    if(fwd.lis.size() > bwd.lis.size())
      fwd.discard_update_LIS(accept_mer, accept_sequence, window_size, L);
    else
      bwd.discard_update_LIS(accept_mer, accept_sequence, window_size, L);
  }
};

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
  super_read_name  unitigs;
  std::vector<int> kmers_info; // Number of k-mers in k-unitigs and common between unitigs
  std::vector<int> bases_info; // Number of bases in k-unitigs and common between unitigs
  double           stretch, offset, avg_err; // Least square stretch, offset and average error
  unsigned int     align_k_;

  coords_info() = default;
  coords_info(const std::string& name, unsigned int align_k, size_t rl, size_t ql, int n) :
    nb_mers(n),
    pb_cons(0), sr_cons(0), pb_cover(align_k), sr_cover(align_k),
    rl(rl), ql(ql), rn(false), qname(name),
    unitigs(name),
    stretch(0), offset(0), avg_err(0),
    align_k_(align_k)
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
    unitigs(qname_),
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
        qs      = ql + qs - align_k_ + 2;
        qe      = ql + qe + 1;
        rn      = true;
        offset -= stretch * (ql + 1) - align_k_;
      } else {
        qs       = -qs + align_k_ - 1;
        qe       = -qe;
        stretch  = -stretch;
        offset  += align_k_ - 1;
      }
    } else {
      qe += align_k_ - 1;
    }
  }

  double imp_s() const { return std::max(1.0, std::min((double)rl, stretch + offset)); }
  double imp_e() const { return std::max(1.0, std::min((double)rl, stretch * ql + offset)); }
  int imp_len() const { return abs(lrint(imp_e() - imp_s())) + 1; }
  bool min_bases(double factor) { return factor * (imp_len() - 2 * (int)align_k_) <= pb_cover; }

  bool min_mers(double factor) { return factor * (imp_len() - align_k_ + 1) <= nb_mers; }
};

typedef std::vector<coords_info>         coords_info_type;
typedef std::map<const char*, mer_lists> frags_pos_type;

// Find all super reads that have k-mers in common with a PacBio
// read. ary contains the k-mers from the super reads, parser is
// initialized to parse the PacBio read. The results are stored in
// frags_pos.
void fetch_super_reads(const mer_pos_hash_type& ary, parse_sequence& parser,
                       frags_pos_type& frags_pos, const int max_mer_count = 0);

// Compute the LIS (Longest Increasing Subsequence) in a list of mers.
template<typename F1, typename F2>
void do_LIS(mer_lists& ml, lis_buffer_type& L, F1& accept_mer, F2& accept_sequence, size_t window_size);

template<typename F1, typename F2>
void do_all_LIS(frags_pos_type& frags_pos, lis_buffer_type& L, F1& accept_mer, F2& accept_sequence, size_t window_size) {
  // Compute LIS forward and backward on every super reads.
  for(auto& it : frags_pos)
    it.second.do_LIS(accept_mer, accept_sequence, window_size, L);
}

template<typename F1, typename F2>
static void do_all_LIS(frags_pos_type& frags_pos, F1& accept_mer, F2& accept_sequence, size_t window_size) {
  lis_buffer_type L;
  do_all_LIS(frags_pos, L, accept_mer, accept_sequence, window_size);
}

// Helper class that computes the number of k-mers in each k-unitigs
// and that are shared between the k-unitigs. It handles errors
// gracefully (the kmers_info array is then empty) and can be a NOOP
// if the super read name is empty.
struct compute_kmers_info {
  std::vector<int>&             mers_;
  std::vector<int>&             bases_;
  const super_read_name&        sr_name_;
  unsigned int                  cunitig_;
  int                           cend_;
  int                           prev_pos_;
  const unsigned int            align_k_;
  const unsigned int            unitigs_k_;
  const std::vector<int>* const unitigs_lengths_;

  compute_kmers_info(std::vector<int>& mers, std::vector<int>& bases, const super_read_name& sr_name,
                     unsigned int unitigs_k, unsigned int align_k, const std::vector<int>* ul);

  size_t nb_unitigs() const { return unitigs_lengths_->size(); }
  int unitig_length(int id) const { return (*unitigs_lengths_)[id]; }

  void add_mer(const int pos);
};


class coarse_aligner {
  const mer_pos_hash_type& ary_;
  const unsigned int       align_k_; // k-mer length used for alignment
  lis_align::affine_capped accept_mer_;
  lis_align::linear        accept_sequence_;
  size_t                   window_size_; // Window to compute stretch
  const bool               forward_;
  const bool               max_match_;
  int                      max_mer_count_; // max mer count to be used for alignment
  double                   matching_mers_factor_;
  double                   matching_bases_factor_;

  const std::vector<int>*  unitigs_lengths_; // Lengths of unitigs
  unsigned int             unitigs_k_; // k-mer length used for creating k-unitigs


  typedef const mer_pos_hash_type::mapped_type list_type;
  typedef const mer_pos_hash_type::position_type position_type;


public:
  coarse_aligner(const mer_pos_hash_type& ary, const unsigned int align_k,
           double stretch_factor, double stretch_constant, double stretch_cap, size_t window_size,
           bool forward = false, bool max_match = false, int max_mer_count = 0,
           double matching_mers = 0.0, double matching_bases = 0.0) :
    ary_(ary),
    align_k_(align_k),
    accept_mer_(stretch_factor, stretch_constant, stretch_cap),
    accept_sequence_(stretch_factor),
    window_size_(window_size),
    forward_(forward),
    max_match_(max_match),
    max_mer_count_(max_mer_count),
    matching_mers_factor_(matching_mers),
    matching_bases_factor_(matching_bases),
    unitigs_lengths_(0), unitigs_k_(0)
  { }

  coarse_aligner& unitigs_lengths(const std::vector<int>* m, unsigned int unitigs_k) {
    if(!forward_) throw std::logic_error("Forward flag must be used if passing unitigs lengths");
    unitigs_lengths_ = m;
    unitigs_k_       = unitigs_k;
    return *this;
  }

  coarse_aligner& max_mer_count(int m) { max_mer_count_ = m; return *this; }

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

  void align_sequence_max(parse_sequence& parser, const size_t pb_size,
                          coords_info_type& coords, frags_pos_type& frags, lis_buffer_type& L) const;
  std::pair<coords_info_type, frags_pos_type> align_sequence_max(parse_sequence& parser, const size_t pb_size) const {
    std::pair<coords_info_type, frags_pos_type> res;
    lis_buffer_type                             L;
    align_sequence_max(parser, pb_size, res.first, res.second, L);
    return res;
  }
  std::pair<coords_info_type, frags_pos_type> align_sequence_max(const std::string& seq) const {
    parse_sequence parser(seq);
    return align_sequence_max(parser, seq.size());
  }

  // Reverse the name of a super read. For example 1R_2F_3F becomes
  // 3R_2R_1F. If the name is not valid, name is returned unchanged.
  static std::string reverse_super_read_name(const std::string& name);

  class thread {
    const coarse_aligner& aligner_;
    frags_pos_type        frags_pos_;
    coords_info_type      coords_;
    lis_buffer_type       L_;

  public:
    thread(const coarse_aligner& a) : aligner_(a) { }

    void align_sequence(parse_sequence& parser, const size_t pb_size);
    void align_sequence(const std::string& seq) {
      parse_sequence parser(seq);
      align_sequence(parser, seq.size());
    }
    void align_sequence_max(parse_sequence& parser, const size_t pb_size);
    void align_sequence_max(const std::string& seq) {
      parse_sequence parser(seq);
      align_sequence_max(parser, seq.size());
    }
    const coords_info_type& coords() const { return coords_; }
    const frags_pos_type& frags_pos() const { return frags_pos_; }
  };
  friend class thread;
};

// class fine_aligner {
//   const unsigned int align_k_;
// public:
//   fine_aligner(unsigned int align_k) : align_k_(align_k) { }
// };

} // namespace align_pb

#endif /* _PB_ALIGNER_HPP_ */
