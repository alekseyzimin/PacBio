#ifndef __FINE_ALIGNER_H__
#define __FINE_ALIGNER_H__

#include <src_jf_aligner/jf_aligner.hpp>
#include <src_jf_aligner/pb_aligner.hpp>
#include <src_jf_aligner/superread_parser.hpp>

namespace align_pb {
struct sr_local_ml {
  double     begin; // Implied first start base of k-mer
  double     end;   // Implied last start base of k-mer
  mer_lists  ml;    // List of mers
  sr_local_ml(const frag_lists::frag_info* frag, double b, double e) : begin(b), end(e), ml(frag) { }
};
typedef std::vector<sr_local_ml> local_mer_lists_type;
typedef std::map<const char*, local_mer_lists_type> frags_local_pos_type;

void fetch_local_super_reads(const sequence_psa& ary, short_parse_sequence& parser,
                             frags_local_pos_type& frags_pos);

class fine_aligner {
  const sequence_psa&           psa_;
  const unsigned int            align_k_;
  const std::vector<int>* const unitigs_lengths_;
  const unsigned int            unitigs_k_;

public:
  fine_aligner(const sequence_psa& psa, unsigned int align_k,
               const std::vector<int>* unitigs_lengths = 0, unsigned int unitigs_k = 0)
    : psa_(psa)
    , align_k_(align_k)
    , unitigs_lengths_(unitigs_lengths)
    , unitigs_k_(unitigs_k)
  { }

  class thread {
    const fine_aligner&       aligner_;
    frags_local_pos_type      frags_pos_;
    lis_buffer_type           L_;
    std::vector<unsigned int> P_;
    coords_info_type          coords_;

  public:
    thread(const fine_aligner& a) : aligner_(a) { }

    template<typename Iterator>
    void prime_frags_pos(Iterator it, const Iterator end) {
      frags_pos_.clear();
      for( ; it != end; ++it) {
        auto&  local_mls = frags_pos_[it->qfrag->fwd.name.c_str()];
        double begin     = std::max((double)0, it->stretch + it->offset - it->avg_err);
        double end       = std::min((double)it->rl, it->stretch * it->ql + it->offset + it->avg_err - aligner_.align_k_);
        local_mls.push_back(sr_local_ml(it->qfrag, begin, end));
      }
    }
    void align_sequence(short_parse_sequence& parser, const size_t pb_size, const coords_info_type& coarse_coords);
    const coords_info_type& coords() { return coords_; }
  };
  friend class thread;
};
} // namespace align_pb

#endif /* __FINE_ALIGNER_H__ */
