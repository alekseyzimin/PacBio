#include <ostream>
#include <string>
#include <vector>
#include <src_jf_aligner/pb_aligner.hpp>


void align_pb::start(int thid) {
  mer_dna         tmp_m;
  parse_sequence  parser(compress_);
  frags_pos_type  frags_pos;
  lis_align::forward_list<lis_align::element<double> > L; // L buffer

  mstream     details(details_multiplexer_);
  mstream     coords(coords_multiplexer_);
  std::string name;

  while(true) {
    read_parser::job job(parser_);
    if(job.is_empty()) break;

    for(size_t i = 0; i < job->nb_filled; ++i) { // Process each read
      auto name_end = job->data[i].header.find_first_of(" \t\n\v\f\r");
      name = job->data[i].header.substr(0, name_end);
      parser.reset(job->data[i].seq);
      frags_pos.clear();
      process_read(ary_, parser, frags_pos, L, stretch_constant_, stretch_factor_, max_mer_count_);
      if(details) print_details(*details, name, frags_pos);
      if(coords) print_coords(*coords, name, job->data[i].seq.size(), frags_pos);
    }
  }
}

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

std::vector<align_pb::coords_info> align_pb::compute_coordinates(const frags_pos_type& frags_pos, const size_t pb_size) {
  std::vector<align_pb::coords_info> coords;
  for(auto it = frags_pos.cbegin(); it != frags_pos.cend(); ++it) {
    const align_pb::mer_lists& ml          = it->second;
    const auto                 fwd_nb_mers = std::distance(ml.fwd.lis.begin(), ml.fwd.lis.end());
    const auto                 bwd_nb_mers = std::distance(ml.bwd.lis.begin(), ml.bwd.lis.end());
    const bool                 fwd_align   = fwd_nb_mers >= bwd_nb_mers;
    const auto                 nb_mers     = fwd_align ? fwd_nb_mers : bwd_nb_mers;

    // Compute consecutive mers and covered bases. Compute hash of
    // positions in pb read of aligned mers. If k_len_ > 0, also
    // compute the number of k-mers aligned in each k-unitigs of the
    // super-read. If an error occurs (unknown k-unitigs, k-unitigs
    // too short, etc.), an empty vector is returned.
    MurmurHash3A                      hasher;
    coords_info                       info(forward_ && !fwd_align ? reverse_super_read_name(ml.frag->name) : ml.frag->name,
                                           ml.frag->len, nb_mers);
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
        hasher.add((uint32_t)cur.first);
        const int pos               = fwd_align ? cur.second : info.ql + cur.second - mer_dna::k() + 2;
        kmers_info.add_mer(pos);
        least_square.add(cur.second, cur.first);
      }
    }

    // Compute average error
    double e = 0;
    if(least_square.n > 0){
      const double a = least_square.a();
      const double b = least_square.b();
      for(auto v : lis) {
        auto& c = offsets[v];
        e += abs(a * c.second + b - c.first);
      }
      e = e / least_square.n;
    }

    uint64_t r1, r2;
    hasher.finalize(r1, r2);
    info.hash  = r1 ^ r2;
    auto first = offsets[lis.front()];
    auto last  = offsets[lis.back()];
    info.rs    = first.first;
    info.re    = last.first + mer_dna::k() - 1;

    info.qs = first.second;
    info.qe = last.second;
    info.stretch = least_square.a();
    info.offset  = least_square.b();
    info.avg_err = e;
    info.canonicalize(forward_);

    if(matching_mers_factor_ || matching_bases_factor_) {
      const double imp_s   = std::max(1.0, std::min((double)pb_size, info.stretch + info.offset));
      const double imp_e   = std::max(1.0, std::min((double)pb_size, info.stretch * info.ql + info.offset));
      const int    imp_len = abs(lrint(imp_e - imp_s)) + 1;
      if(matching_mers_factor_ && matching_mers_factor_ * (imp_len - mer_dna::k() + 1) > lis.size())
        continue;
      if(matching_bases_factor_ && matching_bases_factor_ * (imp_len - 2 * mer_dna::k()) > info.pb_cover)
        continue;
    }

    coords.push_back(info);
  }

  return coords;
}


void align_pb::print_coords(Multiplexer::ostream& out, const std::string& pb_name, size_t pb_size, const frags_pos_type& frags_pos) {
  std::vector<coords_info> coords = compute_coordinates(frags_pos, pb_size);
  auto nb_lines = std::distance(coords.cbegin(), coords.cend());
  if(nb_lines == 0) return;

  std::sort(coords.begin(), coords.end());
  if(compact_format_)
    out << ">" << nb_lines << " " << pb_name << "\n";
  for(auto it = coords.cbegin(), pit = it; it != coords.cend(); pit = it, ++it) {
    if(duplicated_ && it != pit && it->rs == pit->rs && it->re == pit->re && it->hash == pit->hash)
      continue;
    out << it->rs << " " << it->re << " " << it->qs << " " << it->qe << " "
        << it->nb_mers << " "
        << it->pb_cons << " " << it->sr_cons << " "
        << it->pb_cover << " " << it->sr_cover << " "
        << pb_size << " " << it->ql
        << " " << it->stretch << " " << it->offset << " " << it->avg_err;
    if(!compact_format_)
      out << " " << pb_name;
    out << " " << it->qname;
    if(k_len_ != 0 && it->kmers_info.empty()) {
      std::cerr << "Error while finding k-mers matching in k-unitigs. Most likely the lengths of k-unitigs are wrong or the super-read sequences are messed up.";
      exit(1);
    }
    auto mit = it->kmers_info.cbegin();
    auto bit = it->bases_info.cbegin();
    for( ; mit != it->kmers_info.cend(); ++mit, ++bit)
      out << " " << *mit << ":" << *bit;
    out << "\n";
  }
  out.end_record();
}

void align_pb::process_read(const mer_pos_hash_type& ary, parse_sequence& parser,
                            frags_pos_type& frags_pos, lis_align::forward_list<lis_align::element<double> >& L,
                            double a, double b, const int max_mer_count) {
  while(parser.next()) { // Process each k-mer
    const bool is_canonical = parser.m < parser.rm;
    auto it = ary.find_pos(is_canonical ? parser.m : parser.rm);
    if(max_mer_count && std::distance(it, ary.pos_end()) > max_mer_count) continue;
    for( ; it != ary.pos_end(); ++it) {
      mer_lists& ml = frags_pos[it->frag->name];
      ml.frag       = it->frag;
      const int offset = is_canonical ? it->offset : -it->offset;
      if(offset > 0)
        ml.fwd.offsets.push_back(pb_sr_offsets(parser.offset + 1, offset));
      else
        ml.bwd.offsets.push_back(pb_sr_offsets(parser.offset + 1, offset));
    }
  }
  if(frags_pos.empty()) return;

  // Compute LIS forward and backward on every super reads.
  for(auto it = frags_pos.begin(); it != frags_pos.end(); ++it) {
    mer_lists& mer_list = it->second;
    mer_list.fwd.lis.clear();
    mer_list.bwd.lis.clear();
    L.clear();
    lis_align::indices(mer_list.fwd.offsets.cbegin(), mer_list.fwd.offsets.cend(),
                       L, mer_list.fwd.lis, a, b);
    L.clear();
    lis_align::indices(mer_list.bwd.offsets.cbegin(), mer_list.bwd.offsets.cend(),
                       L, mer_list.bwd.lis, a, b);
  }
}


void align_pb_reads(int threads, mer_pos_hash_type& hash, double stretch_const, double stretch_factor,
                    bool compress, bool forward, bool duplicated,
                    const char* pb_path,
                    const char* coords_path, const char* details_path,
                    const int* lengths, const int nb_unitigs, const unsigned int k_len,
                    double matching_mers, double matching_bases) {
  output_file details, coords;
  if(coords_path) coords.open(coords_path, threads);
  if(details_path) details.open(details_path, threads);
  file_vector files;
  files.push_back(pb_path);
  stream_manager streams(files.cbegin(), files.cend());
  align_pb aligner(threads, hash, streams, stretch_const, stretch_factor,
                   forward, compress, duplicated,
                   matching_mers, matching_bases);
  if(coords_path) aligner.coords_multiplexer(coords.multiplexer(), true);
  if(details_path) aligner.details_multiplexer(details.multiplexer());
  std::vector<int> unitigs_lengths;
  if(nb_unitigs && lengths && k_len) {
    unitigs_lengths.insert(unitigs_lengths.begin(), lengths, lengths + nb_unitigs);
    aligner.unitigs_lengths(&unitigs_lengths, k_len);
  }
  aligner.exec_join(threads);
}

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
