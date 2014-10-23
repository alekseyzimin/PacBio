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
  void do_LIS(F1& accept_mer, F2& accept_sequence, size_t window_size, lis_buffer_type& L, std::vector<unsigned int>& P) {
    L.clear();
    lis.clear();
    lis_align::indices(offsets.cbegin(), offsets.cend(), L, P, lis, window_size, accept_mer, accept_sequence);
  }


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
  off_lis                      fwd;
  off_lis                      bwd;
  const frag_lists::frag_info* frag;
  mer_lists() : frag(0) { }
  mer_lists(const frag_lists::frag_info* frag_) : frag(frag_) { }
  template<typename F1, typename F2>
  void do_LIS(F1& accept_mer, F2& accept_sequence, size_t window_size, lis_buffer_type& L) {
    fwd.do_LIS(accept_mer, accept_sequence, window_size, L);
    bwd.do_LIS(accept_mer, accept_sequence, window_size, L);
  }
  template<typename F1, typename F2>
  void do_LIS(F1& accept_mer, F2& accept_sequence, size_t window_size, lis_buffer_type& L, std::vector<unsigned int>& P) {
    fwd.do_LIS(accept_mer, accept_sequence, window_size, L, P);
    bwd.do_LIS(accept_mer, accept_sequence, window_size, L, P);
  }
  template<typename F1, typename F2>
  void discard_update_LIS(F1& accept_mer, F2& accept_sequence, size_t window_size, lis_buffer_type& L) {
    if(fwd.lis.size() > bwd.lis.size())
      fwd.discard_update_LIS(accept_mer, accept_sequence, window_size, L);
    else
      bwd.discard_update_LIS(accept_mer, accept_sequence, window_size, L);
  }
  template<typename F1, typename F2>
  void discard_update_LIS(F1& accept_mer, F2& accept_sequence, size_t window_size, lis_buffer_type& L, std::vector<unsigned int>& P) {
    if(fwd.lis.size() > bwd.lis.size())
      fwd.discard_update_LIS(accept_mer, accept_sequence, window_size, L, P);
    else
      bwd.discard_update_LIS(accept_mer, accept_sequence, window_size, L, P);
  }
};

// Contains the summary information of an alignment of a pac-bio and a super read (named qname).
struct coords_info {
  int                          rs, re;
  int                          qs, qe;
  int                          nb_mers;
  unsigned int                 pb_cons, sr_cons;
  unsigned int                 pb_cover, sr_cover;
  size_t                       rl, ql;
  bool                         rn;
  const frag_lists::frag_info* qfrag;
  super_read_name              unitigs;
  std::vector<int>             kmers_info; // Number of k-mers in k-unitigs and common between unitigs
  std::vector<int>             bases_info; // Number of bases in k-unitigs and common between unitigs
  double                       stretch, offset, avg_err; // Least square stretch, offset and average error
  unsigned int                 align_k_;

  coords_info() = default;
  coords_info(const frag_lists::frag_info* frag, unsigned int align_k, size_t rl, size_t ql, int n) :
    nb_mers(n),
    pb_cons(0), sr_cons(0), pb_cover(align_k), sr_cover(align_k),
    rl(rl), ql(ql), rn(false), qfrag(frag),
    unitigs(frag->name),
    stretch(0), offset(0), avg_err(0),
    align_k_(align_k)
  { }
  coords_info(int rs_, int re_, int qs_, int qe_, int nb_mers_,
              unsigned int pb_cons_, unsigned sr_cons_,
              unsigned int pb_cover_, unsigned int sr_cover_,
              size_t rl_, size_t ql_, bool rn_,
              double stretch_, double offset_, double avg_err_,
              const frag_lists::frag_info* frag) :
    rs(rs_), re(re_), qs(qs_), qe(qe_), nb_mers(nb_mers_),
    pb_cons(pb_cons_), sr_cons(sr_cons_),
    pb_cover(pb_cover_), sr_cover(sr_cover_),
    rl(rl_), ql(ql_), rn(rn_), qfrag(frag),
    unitigs(frag->name),
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

// Compute the coords info from the mer list (after LIS).
coords_info compute_coords_info(const mer_lists& ml, const size_t pb_size,
                                const unsigned int align_k,
                                const unsigned int unitigs_k,
                                const std::vector<int>* const unitigs_lengths,
                                const bool forward);

// Compute the LIS (Longest Increasing Subsequence) in a list of mers.
template<typename F1, typename F2>
void do_LIS(mer_lists& ml, lis_buffer_type& L, F1& accept_mer, F2& accept_sequence, size_t window_size);

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
} // namespace align_pb

#endif /* _PB_ALIGNER_HPP_ */
