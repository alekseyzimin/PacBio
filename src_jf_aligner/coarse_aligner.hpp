/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef __COARSE_ALIGNER_H__
#define __COARSE_ALIGNER_H__

#include <unordered_map>
#include <queue>
#include <algorithm>
#include <src_jf_aligner/pb_aligner.hpp>
#include <src_jf_aligner/superread_parser.hpp>

namespace align_pb {
typedef std::unordered_map<const char*, mer_lists> frags_pos_type;

template<typename F1, typename F2>
void do_all_LIS(frags_pos_type& frags_pos, lis_buffer_type& L, std::vector<unsigned int>& P,
                const F1& accept_mer, const F2& accept_sequence, size_t window_size) {
  // Compute LIS forward and backward on every super reads.
  for(auto& it : frags_pos)
    it.second.do_LIS(accept_mer, accept_sequence, window_size, L, P);
}

template<typename F1, typename F2>
static void do_all_LIS(frags_pos_type& frags_pos, const F1& accept_mer, const F2& accept_sequence, size_t window_size) {
  lis_buffer_type L;
  std::vector<unsigned int> P;
  do_all_LIS(frags_pos, L, P, accept_mer, accept_sequence, window_size);
}

// Find all super reads that have k-mers in common with a PacBio
// read. psa contains the (partial) suffix array from the super reads,
// parser is initialized to parse the PacBio read. The results are
// stored in frags_pos.
void fetch_super_reads(const sequence_psa& psa, parse_sequence& parser,
                       frags_pos_type& frags_pos, const int max_mer_count = 0, const float max_percent = 100.0f);

class coarse_aligner {
  const sequence_psa&      psa_;
  const unsigned int       align_k_; // k-mer length used for alignment
  lis_align::affine_capped accept_mer_;
  lis_align::linear        accept_sequence_;
  size_t                   window_size_; // Window to compute stretch
  const bool               forward_;
  const bool               max_match_;
  int                      max_mer_count_; // max mer count to be used for alignment
  double                   matching_mers_factor_;
  double                   matching_bases_factor_;
  double                   max_percent_; // Compute mer count threshold to keep max_percent_% mers

  const std::vector<int>*  unitigs_lengths_; // Lengths of unitigs
  unsigned int             unitigs_k_; // k-mer length used for creating k-unitigs

public:
  coarse_aligner(const sequence_psa& psa, const unsigned int align_k,
                 double stretch_factor, double stretch_constant, double stretch_cap, size_t window_size,
                 bool forward = false, bool max_match = false, int max_mer_count = 0,
                 double matching_mers = 0.0, double matching_bases = 0.0)
    : psa_(psa)
    , align_k_(align_k)
    , accept_mer_(stretch_factor, stretch_constant, stretch_cap)
    , accept_sequence_(stretch_factor)
    , window_size_(window_size)
    , forward_(forward)
    , max_match_(max_match)
    , max_mer_count_(max_mer_count)
    , matching_mers_factor_(matching_mers)
    , matching_bases_factor_(matching_bases)
    , max_percent_(100.0f)
    , unitigs_lengths_(0)
    , unitigs_k_(0)
  { }

  coarse_aligner max_percent(double v) { max_percent_ = v; return *this; }

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
                      coords_info_type& coords, frags_pos_type& frags, lis_buffer_type& L, std::vector<unsigned int>& P) const;
  std::pair<coords_info_type, frags_pos_type> align_sequence(parse_sequence& parser, const size_t pb_size) const {
    std::pair<coords_info_type, frags_pos_type> res;
    lis_buffer_type                             L;
    std::vector<unsigned int>                   P;
    align_sequence(parser, pb_size, res.first, res.second, L, P);
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
    const coarse_aligner&     aligner_;
    frags_pos_type            frags_pos_;
    coords_info_type          coords_;
    lis_buffer_type           L_;
    std::vector<unsigned int> P_;

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


// Compute the n-th lowest element
template<typename T>
struct nth_min {
  std::priority_queue<T> queue;
  const size_t n;
  nth_min(size_t s) : n(s) { }
  void add(const T& x) {
    if(queue.size() < n || x < queue.top()) {
      if(queue.size() == n)
        queue.pop();
      queue.push(x);
    }
  }
  T highest() const { return queue.top(); }
};

} // namespace align_pb

#endif /* __COARSE_ALIGNER_H__ */
