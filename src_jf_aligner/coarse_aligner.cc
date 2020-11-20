/******************************************
Copyright University of Maryland 2015
******************************************/
#include <src_jf_aligner/coarse_aligner.hpp>
#include <cmath>


bool is_ssr(const jellyfish::mer_dna& m,uint32_t l=2) {
  jellyfish::mer_dna nm(m);
  for(uint32_t i = 0; i < l; ++i) {
    nm.shift_right(nm.base(0));
    if(nm == m) return true;
  }
  return false;
}

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
  fetch_super_reads(psa_, parser, frags_pos, max_mer_count_, max_percent_);
  do_all_LIS(frags_pos, L, P, accept_mer_, accept_sequence_, window_size_);
  compute_coords(frags_pos, pb_size, coords);
}

void coarse_aligner::align_sequence_max (parse_sequence& parser, const size_t pb_size,
                                   coords_info_type& coords, frags_pos_type& frags_pos,
                                   lis_buffer_type& L) const {
  fetch_super_reads(psa_, parser, frags_pos, max_mer_count_, max_percent_);
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

struct mer_list_info {
  sequence_psa::pos_iterator list;
  size_t                     size;
  bool                       is_canonical;
  int                        offset;
};

void fetch_super_reads(const sequence_psa& psa, parse_sequence& parser,
                       frags_pos_type& frags_pos, const int max_mer_count, const float max_percent) {
  const auto end = psa.pos_end();

  std::vector<mer_list_info> lists_info;
  uint32_t counts[max_mer_count + 1];

  memset(counts, '\0', sizeof(counts));
  uint32_t flag=0;
  while(parser.next()) { // Process each k-mer

    //skip if k-mer is low complexity 2-simple sequence repeat (e.g. either AAAAAAAAAAAAAAA or ATATATATATATATATA)
    if(is_ssr(parser.mer<0>().m,2)) continue;
    //if(parser.mer<0>().m.is_homopolymer()) continue;
    //we take every other k-mer in the long read to reduce the number of calls to psa.find_pos_size
    flag=1-flag;
    if(flag==0) continue;

    const bool is_canonical = parser.mer<0>().is_canonical();
    auto list = is_canonical
      ? psa.find_pos_size(parser.mer<0>().m, parser.mer<0>().rm)
      : psa.find_pos_size(parser.mer<0>().rm, parser.mer<0>().m);
    if(list.second == 0 ||
       (max_mer_count && list.second >= (size_t)max_mer_count)) { // && std::distance(list.first, end) >= max_mer_count))
      continue;
    }
    ++counts[std::min(max_mer_count, (int)list.second)];
    lists_info.push_back({list.first, list.second, is_canonical, parser.offset<0>()});
  }

  uint32_t sum        = 0;
  uint32_t sum_thresh = std::round(lists_info.size() * 0.99);
  uint32_t threshold  = 1;

  for( ; threshold <= (uint32_t)max_mer_count; ++threshold) {
    sum += counts[threshold];
    if(sum > sum_thresh)
      break;
  }
  //std::cerr << threshold <<'\n';
  
  for(auto& info : lists_info) {
    if(info.size > threshold)
      continue;
    for(auto& it = info.list; it != end; ++it) { // For each instance of the k-mer in a super read
      mer_lists& ml = frags_pos[it->frag->fwd.name.c_str()];
      ml.frag       = it->frag;
      const int offset = info.is_canonical ? it->offset : -it->offset;
      if(offset > 0)
        ml.fwd.offsets.push_back(pb_sr_offsets(info.offset, offset));
      else
        ml.bwd.offsets.push_back(pb_sr_offsets(info.offset, offset));
    }
  }
}

} // namespace align_pb
