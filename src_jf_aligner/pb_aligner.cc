#include <ostream>
#include <string>
#include <vector>
#include <src_jf_aligner/pb_aligner.hpp>

void align_pb::compute_coords(const frags_pos_type& frags_pos, const size_t pb_size,
                              coords_info_type& coords) const {
  for(const auto& it : frags_pos) {
    coords_info info = compute_coords_info(it.second, pb_size);
    if(matching_mers_factor_ && !info.min_mers(matching_mers_factor_)) continue;
    if(matching_bases_factor_ > 0.0 && !info.min_bases(matching_bases_factor_)) continue;
    coords.push_back(std::move(info));
  }
}

align_pb::coords_info align_pb::compute_coords_info(const mer_lists& ml, const size_t pb_size) const {
    const auto                 fwd_nb_mers = std::distance(ml.fwd.lis.begin(), ml.fwd.lis.end());
    const auto                 bwd_nb_mers = std::distance(ml.bwd.lis.begin(), ml.bwd.lis.end());
    const bool                 fwd_align   = fwd_nb_mers >= bwd_nb_mers;
    const auto                 nb_mers     = fwd_align ? fwd_nb_mers : bwd_nb_mers;

    // Compute consecutive mers and covered bases. Compute hash of
    // positions in pb read of aligned mers. If k_len_ > 0, also
    // compute the number of k-mers aligned in each k-unitigs of the
    // super-read. If an error occurs (unknown k-unitigs, k-unitigs
    // too short, etc.), an empty vector is returned.
    coords_info                       info(ml.frag->name, pb_size, ml.frag->len, nb_mers);
    if(forward_ && !fwd_align)
      info.unitigs.reverse();
    if(nb_mers == 0) return info;
    const std::vector<pb_sr_offsets>& offsets = fwd_align ? ml.fwd.offsets : ml.bwd.offsets;
    const std::vector<unsigned int>&  lis     = fwd_align ? ml.fwd.lis : ml.bwd.lis;
    compute_kmers_info                kmers_info(info.kmers_info, info.bases_info, info.unitigs, k_len_, unitigs_lengths_);
    least_square_2d                   least_square;
    {
      auto          lisit = lis.cbegin();
      pb_sr_offsets prev  = offsets[*lisit];
      pb_sr_offsets cur;
      const int     pos   = fwd_align ? prev.second : info.ql + prev.second - mer_dna::k() + 2;
      kmers_info.add_mer(pos);
      least_square.add(prev.second, prev.first);
      for(++lisit; lisit != lis.cend(); prev = cur, ++lisit) {
        cur                         = offsets[*lisit];
        const unsigned int pb_diff  = cur.first - prev.first;
        info.pb_cons               += pb_diff == 1;
        info.pb_cover              += std::min(mer_dna::k(), pb_diff);
        const unsigned int sr_diff  = cur.second - prev.second;
        info.sr_cons               += sr_diff == 1;
        info.sr_cover              += std::min(mer_dna::k(), sr_diff);
        const int pos               = fwd_align ? cur.second : info.ql + cur.second - mer_dna::k() + 2;
        kmers_info.add_mer(pos);
        least_square.add(cur.second, cur.first);
      }
    }

    // Compute average error
    double e = 0;
    if(least_square.n == 1) {
      // In that case, compute offset to shift the 1 mer from the
      // super read coordinate to the PB coordinate
      info.stretch = 1.0;
      info.offset  = least_square.EY - least_square.EX;
      info.avg_err = 0;
    } else if(least_square.n > 1) {
      const double a = info.stretch = least_square.a();
      const double b = info.offset = least_square.b();
      for(auto v : lis) {
        auto& c = offsets[v];
        e += abs(a * c.second + b - c.first);
      }
      info.avg_err = e / least_square.n;
    }

    auto first = offsets[lis.front()];
    auto last  = offsets[lis.back()];
    info.rs    = first.first;
    info.re    = last.first + mer_dna::k() - 1;

    info.qs = first.second;
    info.qe = last.second;
    info.canonicalize(forward_);

    return info;
}

void align_pb::align_sequence(parse_sequence& parser, const size_t pb_size,
                              coords_info_type& coords, frags_pos_type& frags_pos,
                              lis_buffer_type& L) const {
  fetch_super_reads(ary_, parser, frags_pos, max_mer_count_);

  do {
    do_all_LIS(frags_pos, L, stretch_constant_, stretch_factor_);
    compute_coords(frags_pos, pb_size, coords);
  } while(false);
}

void align_pb::align_sequence_max (parse_sequence& parser, const size_t pb_size,
                                   coords_info_type& coords, frags_pos_type& frags_pos,
                                   lis_buffer_type& L) const {
  fetch_super_reads(ary_, parser, frags_pos, max_mer_count_);
  for(auto& it : frags_pos) {
    auto& ml = it.second;
    ml.do_LIS(stretch_constant_, stretch_factor_, L);
    while(true) {
      coords_info info = compute_coords_info(ml, pb_size);
      if(info.nb_mers == 0) break;
      if(matching_mers_factor_ && !info.min_mers(matching_mers_factor_)) break;
      if(matching_bases_factor_ > 0.0 && !info.min_bases(matching_bases_factor_)) break;
      coords.push_back(std::move(info));
      if(!max_match_) break;
      ml.discard_update_LIS(stretch_constant_, stretch_factor_, L);
    }
  }
}


void align_pb::fetch_super_reads(const mer_pos_hash_type& ary, parse_sequence& parser,
                                 frags_pos_type& frags_pos, const int max_mer_count) {
  while(parser.next()) { // Process each k-mer
    const bool is_canonical = parser.m < parser.rm;
    auto it = ary.find_pos(is_canonical ? parser.m : parser.rm);
    if(max_mer_count && std::distance(it, ary.pos_end()) > max_mer_count) continue;
    for( ; it != ary.pos_end(); ++it) { // For each instance of the k-mer in a super read
      mer_lists& ml = frags_pos[it->frag->name];
      ml.frag       = it->frag;
      const int offset = is_canonical ? it->offset : -it->offset;
      if(offset > 0)
        ml.fwd.offsets.push_back(pb_sr_offsets(parser.offset, offset));
      else
        ml.bwd.offsets.push_back(pb_sr_offsets(parser.offset, offset));
    }
  }
}

void align_pb::do_all_LIS(frags_pos_type& frags_pos, lis_buffer_type& L, double a, double b) {
  // Compute LIS forward and backward on every super reads.
  for(auto& it : frags_pos)
    it.second.do_LIS(a, b, L);
}

void align_pb::thread::align_sequence(parse_sequence& parser, const size_t pb_size) {
  frags_pos_.clear();
  coords_.clear();
  align_data_.align_sequence(parser, pb_size, coords_, frags_pos_, L_);
}

void align_pb::thread::align_sequence_max(parse_sequence& parser, const size_t pb_size) {
  frags_pos_.clear();
  coords_.clear();
  align_data_.align_sequence_max(parser, pb_size, coords_, frags_pos_, L_);
}

align_pb::compute_kmers_info::compute_kmers_info(std::vector<int>& mers, std::vector<int>& bases, const super_read_name& sr_name,
                                                 unsigned int k_len, const std::vector<int>* ul) :
  mers_(mers), bases_(bases), sr_name_(sr_name),
  cunitig_(0), cend_(0), prev_pos_(-mer_dna::k()),
  k_len_(k_len), unitigs_lengths_(ul)
{
  if(k_len_) {
    const auto unitig_id = sr_name_.unitig_id(0);
    if(unitig_id != super_read_name::invalid && unitig_id < nb_unitigs()) {
      mers_.resize(2 * sr_name_.size() - 1, 0);
      bases_.resize(2 * sr_name_.size() - 1, 0);
      cend_ = unitig_length(unitig_id);
    } else {
      mers_.clear();
      bases_.clear();
    }
  }
}

void align_pb::compute_kmers_info::add_mer(const int pos) {
  if(!k_len_) return;

  int       cendi;
  const int sr_pos = abs(pos);
  const int new_bases   = std::min((int)mer_dna::k(), sr_pos - prev_pos_);
  while(sr_pos + (int)mer_dna::k() > cend_ + 1) {
    if(cend_ >= sr_pos) {
      const int nb_bases = cend_ - std::max(sr_pos, prev_pos_ + (int)mer_dna::k()) + 1;
      bases_[2 * cunitig_]     += nb_bases;
      bases_[2 * cunitig_ + 1] += nb_bases;
    }
    const auto unitig_id = sr_name_.unitig_id(++cunitig_);
    if(unitig_id == super_read_name::invalid || unitig_id >= nb_unitigs())
      goto error;
    cend_ += unitig_length(unitig_id) - k_len_+ 1;
  }
  ++mers_[2 * cunitig_];
  bases_[2 * cunitig_] += new_bases;
  cendi                 = cend_;
  for(unsigned int i = cunitig_; (i < sr_name_.size() - 1) && (sr_pos + mer_dna::k() > cendi - k_len_ + 1); ++i) {
    const int  full_mer   = sr_pos + (int)k_len_ > cendi + 1;
    mers_[2 * i + 1]     += full_mer;
    mers_[2 * i + 2]     += full_mer;
    const int  nb_bases   = std::min(new_bases, sr_pos + (int)mer_dna::k() - cendi + (int)k_len_ - 2);
    bases_[2 * i + 1]    += nb_bases;
    bases_[2 * i + 2]    += nb_bases;
    const auto unitig_id  = sr_name_.unitig_id(i + 1);
    if(unitig_id != super_read_name::invalid && unitig_id < nb_unitigs())
      cendi += unitig_length(unitig_id) - k_len_ + 1;
    else
      goto error;
  }
  prev_pos_ = sr_pos;
  return;

 error:
  mers_.clear();
  bases_.clear();
}
