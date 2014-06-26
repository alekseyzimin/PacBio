#include <ostream>
#include <string>
#include <vector>
#include <src_jf_aligner/pb_aligner.hpp>


void align_pb::print_details(Multiplexer::ostream& out, const std::string& pb_name, const frags_pos_type& frags_pos) {
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

void align_pb::print_coords_header(Multiplexer* m, bool compact) {
  Multiplexer::ostream o(m); // Write header
  o << "Rstart Rend Qstart Qend Nmers Rcons Qcons Rcover Qcover Rlen Qlen Stretch Offset Err";
  if(!compact)
    o << " Rname";
  o << " Qname\n";
  o.end_record();
}


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
    coords_info                       info(forward_ && !fwd_align ? reverse_super_read_name(ml.frag->name) : ml.frag->name,
                                           pb_size, ml.frag->len, nb_mers);
    if(nb_mers == 0) return info;
    const std::vector<pb_sr_offsets>& offsets = fwd_align ? ml.fwd.offsets : ml.bwd.offsets;
    const std::vector<unsigned int>&  lis     = fwd_align ? ml.fwd.lis : ml.bwd.lis;
    compute_kmers_info<align_pb>      kmers_info(info.kmers_info, info.bases_info, info.qname, *this);
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
  std::sort(coords.begin(), coords.end());
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
  std::sort(coords.begin(), coords.end());
}

void align_pb::print_coords(Multiplexer::ostream& out, const std::string& pb_name, const size_t pb_size,
                            const bool compact, const coords_info_type& coords) {
  auto nb_lines = std::distance(coords.cbegin(), coords.cend());
  if(nb_lines == 0) return;

  if(compact)
    out << ">" << nb_lines << " " << pb_name << "\n";
  for(const auto& it : coords) {
    out << it.rs << " " << it.re << " " << it.qs << " " << it.qe << " "
        << it.nb_mers << " "
        << it.pb_cons << " " << it.sr_cons << " "
        << it.pb_cover << " " << it.sr_cover << " "
        << pb_size << " " << it.ql
        << " " << it.stretch << " " << it.offset << " " << it.avg_err;
    if(!compact)
      out << " " << pb_name;
    out << " " << it.qname;
    auto mit = it.kmers_info.cbegin();
    auto bit = it.bases_info.cbegin();
    for( ; mit != it.kmers_info.cend(); ++mit, ++bit)
      out << " " << *mit << ":" << *bit;
    out << "\n";
  }
  out.end_record();
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

// static void align_pb::align_pb_reads(mer_pos_hash_type& hash, double stretch_const, double stretch_factor,
//                                      bool forward, const char* pb_path,
//                                      const int* lengths, const int nb_unitigs, const unsigned int k_len,
//                                      double matching_mers, double matching_bases) {
//   file_vector files;
//   files.push_back(pb_path);
//   stream_manager streams(files.cbegin(), files.cend());
//   align_pb aligner(threads, hash, streams, stretch_const, stretch_factor,
//                    forward,
//                    matching_mers, matching_bases);
//   if(coords_path) aligner.coords_multiplexer(coords.multiplexer(), true);
//   if(details_path) aligner.details_multiplexer(details.multiplexer());
//   std::vector<int> unitigs_lengths;
//   if(nb_unitigs && lengths && k_len) {
//     unitigs_lengths.insert(unitigs_lengths.begin(), lengths, lengths + nb_unitigs);
//     aligner.unitigs_lengths(&unitigs_lengths, k_len);
//   }
//   aligner.exec_join(threads);
// }

std::string align_pb::reverse_super_read_name(const std::string& name) {
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
