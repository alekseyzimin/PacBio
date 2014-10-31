#include <src_jf_aligner/fine_aligner.hpp>

namespace align_pb {
void fetch_local_super_reads(const short_mer_pos_hash_type& ary, short_parse_sequence& parser,
                             frags_local_pos_type& frags_pos) {
  while(parser.next()) { // Process each k-mer
    const bool is_canonical = parser.mer<0>().is_canonical();
    auto list = ary.equal_range(is_canonical ? parser.mer<0>().m : parser.mer<0>().rm);

    for(auto it = list.first ; it != list.second; ++it) { // For each instance of the k-mer in a super read
      auto local_mls = frags_pos.find(it->frag->fwd.name.c_str());
      if(local_mls == frags_pos.end()) continue;
      const auto mls_end = local_mls->second.end();
      for(auto local_ml = local_mls->second.begin(); local_ml != mls_end; ++local_ml) {
        if(parser.offset<0>() >= local_ml->begin && parser.offset<0>() <= local_ml->end) {
          auto&     ml     = local_ml->ml;
          //          ml.frag          = it->frag;
          const int offset = is_canonical ? it->offset : -it->offset;
          if(offset > 0)
            ml.fwd.offsets.push_back(pb_sr_offsets(parser.offset<0>(), offset));
          else
            ml.bwd.offsets.push_back(pb_sr_offsets(parser.offset<0>(), offset));
        }
      }
    }
  }

}

void fine_aligner::thread::align_sequence(short_parse_sequence& parser, const size_t pb_size, const coords_info_type& coarse_coords) {
  prime_frags_pos(coarse_coords.cbegin(), coarse_coords.cend());
  coords_.clear();
  fetch_local_super_reads(aligner_.ary_, parser, frags_pos_);

  lis_align::accept_all all_mers;
  for(auto& it : frags_pos_) {
    for(auto& local_ml : it.second) {
      local_ml.ml.do_LIS(all_mers, all_mers, 1, L_, P_);
      coords_.push_back(compute_coords_info(local_ml.ml, pb_size, aligner_.align_k_, aligner_.unitigs_k_,
                                            aligner_.unitigs_lengths_, true));
    }
  }
}
} // namespace align_pb
