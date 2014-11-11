#include <src_jf_aligner/coarse_aligner.hpp>

namespace align_pb {
void coarse_aligner::compute_coords(const frags_pos_type& frags_pos, const size_t pb_size,
                              coords_info_type& coords) const {
  for(const auto& it : frags_pos) {
    coords_info info = compute_coords_info(it.second, pb_size);
    if(fabs(info.stretch) == 0.0) continue; // Very compressed. Collapsed repeats.
    if(matching_mers_factor_ && !info.min_mers(matching_mers_factor_)) continue;
    if(matching_bases_factor_ > 0.0 && !info.min_bases(matching_bases_factor_)) continue;
    coords.push_back(std::move(info));
  }
}

coords_info coarse_aligner::compute_coords_info(const mer_lists& ml, const size_t pb_size) const {
  return ::align_pb::compute_coords_info(ml, pb_size, align_k_, unitigs_k_,
                                         unitigs_lengths_, forward_);
}

void coarse_aligner::align_sequence(parse_sequence& parser, const size_t pb_size,
                                    coords_info_type& coords, frags_pos_type& frags_pos,
                                    lis_buffer_type& L, std::vector<unsigned int>& P) const {
  fetch_super_reads(ary_, parser, frags_pos, max_mer_count_);
  do_all_LIS(frags_pos, L, P, accept_mer_, accept_sequence_, window_size_);
  compute_coords(frags_pos, pb_size, coords);
}

void coarse_aligner::align_sequence_max (parse_sequence& parser, const size_t pb_size,
                                   coords_info_type& coords, frags_pos_type& frags_pos,
                                   lis_buffer_type& L) const {
  fetch_super_reads(ary_, parser, frags_pos, max_mer_count_);
  for(auto& it : frags_pos) {
    auto& ml = it.second;
    ml.do_LIS(accept_mer_, accept_sequence_, window_size_, L);
    while(true) {
      coords_info info = compute_coords_info(ml, pb_size);
      if(info.nb_mers == 0) break;
      if(fabs(info.stretch) == 0.0) break;
      if(matching_mers_factor_ && !info.min_mers(matching_mers_factor_)) break;
      if(matching_bases_factor_ > 0.0 && !info.min_bases(matching_bases_factor_)) break;
      coords.push_back(std::move(info));
      if(!max_match_) break;
      ml.discard_update_LIS(accept_mer_, accept_sequence_, window_size_, L);
    }
  }
}

void coarse_aligner::thread::align_sequence(parse_sequence& parser, const size_t pb_size) {
  frags_pos_.clear();
  coords_.clear();
  aligner_.align_sequence(parser, pb_size, coords_, frags_pos_, L_, P_);
}

void coarse_aligner::thread::align_sequence_max(parse_sequence& parser, const size_t pb_size) {
  frags_pos_.clear();
  coords_.clear();
  aligner_.align_sequence_max(parser, pb_size, coords_, frags_pos_, L_);
}

void fetch_super_reads(const mer_pos_hash_type& ary, parse_sequence& parser,
                       frags_pos_type& frags_pos, const int max_mer_count) {
  while(parser.next()) { // Process each k-mer
    const bool is_canonical = parser.mer<0>().is_canonical();
    auto list = ary.find_pos_size(is_canonical ? parser.mer<0>().m : parser.mer<0>().rm);
    if(max_mer_count && list.second >= max_mer_count) continue;
    const auto end = ary.pos_end();
    for(auto it = list.first ; it != end; ++it) { // For each instance of the k-mer in a super read
      mer_lists& ml = frags_pos[it->frag->fwd.name.c_str()];
      ml.frag       = it->frag;
      const int offset = is_canonical ? it->offset : -it->offset;
      if(offset > 0)
        ml.fwd.offsets.push_back(pb_sr_offsets(parser.offset<0>(), offset));
      else
        ml.bwd.offsets.push_back(pb_sr_offsets(parser.offset<0>(), offset));
    }
  }
}

} // namespace align_pb
