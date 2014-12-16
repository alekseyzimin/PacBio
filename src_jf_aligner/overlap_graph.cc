#include <src_jf_aligner/overlap_graph.hpp>
#include <boost/icl/interval_set.hpp>

void overlap_graph::traverse(const std::vector<int>& sort_array, const align_pb::coords_info_type& coords,
                             std::vector<node_info>& nodes, std::ostream* dot) const {
  for(size_t i = 0; i != sort_array.size(); ++i) {
    const size_t it_i     = sort_array[i];
    auto&        node_i   = nodes[it_i];
    const auto&  coords_i = coords[it_i];

    if(node_i.imp_e >= coords_i.rl) continue; // Hanging off 3' end of PacBio read: no overlap
    for(size_t j = i + 1; j != sort_array.size(); ++j) {
      const size_t it_j     = sort_array[j];
      auto&        node_j   = nodes[it_j];
      const auto&  coords_j = coords[it_j];

      if(node_j.imp_s <= 1) continue; // Hanging off 5' end of PacBio read: no overlap

      const double position_len = node_i.imp_e - node_j.imp_s;
      if(position_len * overlap_play < k_len) break; // maximum implied overlap is less than a k-mer length
      const int nb_u_overlap = coords_i.name_u->unitigs.overlap(coords_j.name_u->unitigs);
      if(!nb_u_overlap) continue; // No overlap according to unitig names
      if(coords_i.name_u->unitigs == coords_j.name_u->unitigs) continue;
      int u_overlap_len  = 0;
      int common_overlap = 0;
      for(int u = 0; u < nb_u_overlap; ++u) {
        u_overlap_len    += unitigs_lengths[coords_j.name_u->unitigs.unitig_id(u)];
        common_overlap   += maximize_bases ? coords_j.bases_info[2 * u] : coords_j.kmers_info[2 * u];
        if(u > 0)
          common_overlap -= maximize_bases ? coords_j.bases_info[2 * u - 1] : coords_j.kmers_info[2 * u - 1];
      }
      u_overlap_len -= (nb_u_overlap - 1) * (k_len - 1);
      double error = nb_errors * (coords_i.avg_err + coords_j.avg_err);
      if(u_overlap_len > overlap_play * position_len + error || position_len > overlap_play * u_overlap_len + error)
        continue; // Overlap lengths (position, unitigs) do not agree

      // We have an overlap between nodes i and j
      node_i.end_node    = false;
      node_j.start_node  = false;
      node_i.component  |= node_j.component;

      // Update longest path
      int nlpath = node_i.lpath + (maximize_bases ? coords_j.sr_cover : coords_j.nb_mers) - common_overlap;
      if(nlpath > node_j.lpath ||
         (nlpath == node_j.lpath && (node_j.lstart == -1 || node_i.l_start_node(nodes).imp_s > node_j.l_start_node(nodes).imp_s))) {
        node_j.lpath  = nlpath;
        node_j.lstart = node_i.lstart == -1 ? it_i : node_i.lstart;
        node_j.lprev  = it_i;
        node_j.lunitigs = node_i.lunitigs + coords_j.name_u->unitigs.size() - nb_u_overlap;
      }
      if(dot)
        *dot << "n" << it_i << " -> n" << it_j << " [tooltip=\"" << "..." << "\", label=\"" << common_overlap << "\"];\n";
    }
  }
}

mega_read_info mega_read_info::make(int i, const std::vector<node_info>& nodes,
                                    const align_pb::coords_info_type&  coords) {
  mega_read_info res;
  res.start_node   = nodes[i].lstart == -1 ? i : nodes[i].lstart;
  res.end_node     = i;
  res.start_unitig = 0;
  res.nb_unitigs   = nodes[i].lunitigs;
  res.end_unitig   = coords[i].kmers_info.size() / 2;
  res.tiling_start = coords[res.start_node].rs;
  res.tiling_end   = coords[i].re;
  return res;
}

void overlap_graph::update_mr_trim(mega_read_info& mr, std::vector<node_info>& nodes,
                                   const align_pb::coords_info_type& coords) const {
  {
    const auto& coord = coords[mr.start_node];
    int offset        = 0;
    for(mr.start_unitig = 0; mr.start_unitig < (int)coord.kmers_info.size(); mr.start_unitig += 2) {
      if(coord.kmers_info[mr.start_unitig])
        break;
      offset += unitigs_lengths[mr.start_unitig / 2];
    }
    mr.start_unitig /= 2;
    mr.nb_unitigs   -= mr.start_unitig;
    mr.tiling_start  = coord.stretch * (coord.qs + offset) + coord.offset;
  }

  {
    const auto& coord  = coords[mr.end_node];
    int         offset = 0;
    for(mr.end_unitig = coord.kmers_info.size() - 1; mr.end_unitig >= 0; mr.end_unitig -= 2) {
      if(coord.kmers_info[mr.end_unitig])
        break;
      offset += unitigs_lengths[mr.end_unitig / 2];
    }
    mr.end_unitig /= 2;
    mr.nb_unitigs += mr.end_unitig - (int)(coord.kmers_info.size() / 2);
    mr.tiling_end  = coord.stretch * (coord.qe - offset) + coord.offset;
  }
}

void overlap_graph::mega_reads_per_comp(const int n, size_t pb_size, std::vector<node_info>& nodes,
                                        const align_pb::coords_info_type& coords, comp_to_path& components,
                                        double min_density, double min_len, std::ostream* dot) const {
  // For each connected component, keep the index of the terminal node
  // of the longest path found
  for(int i = 0; i < n; ++i) {
    auto&        node    = nodes[i];
    const double imp_len = std::min((double)pb_size + 0.5, node.imp_e) - std::max(0.5, node.l_start_node(nodes).imp_s);
    auto         mr      = mega_read_info::make(i, nodes, coords);
    update_mr_trim(mr, nodes, coords);
    mr.density = (double)node.lpath / (mr.tiling_end - mr.tiling_start);
    //    double       density = (double)node.lpath / imp_len;
    if(dot) {
      const char* color = "";
      if(node.start_node) {
        color = ", color=\"blue\"";
      } else if(node.end_node) {
        color = ", color=\"green\"";
      }
      const auto& ci = coords[i];
      *dot << std::fixed << "n" << i << " [label=\"" << i << " L" << ci.ql << " #" << ci.nb_mers
           << "\\nP(" << ci.rs << ',' << ci.re << ") S(" << ci.qs << ',' << ci.qe << ")"
           << "\\nI(" << node.imp_s << ',' << node.imp_e << ")"
           << "\\nLP #" << node.lpath << " L" << std::setprecision(1) << imp_len
           << " d" << std::setprecision(2) << mr.density << "\""
           << color << "];\n";
    }
    if(!node.end_node || mr.density < min_density || (mr.tiling_end - mr.tiling_start) - min_len) continue;
    // || node.ldensity < min_density || imp_len < min_len) continue;

    auto comp_root = node.component.root();
    auto comp_it   = components.find(comp_root);
    if(comp_it == components.end()) {
      components.insert(std::make_pair(comp_root, mr));
    } else {
      const auto& onode = nodes[comp_it->second.end_node]; // current terminal node of longest path
      if(node.lpath > onode.lpath || (node.lpath == onode.lpath && mr.density > comp_it->second.density))
        comp_it->second = mr;
    }
  }
}

typedef boost::icl::right_open_interval<double> pos_interval;
typedef boost::icl::interval_set<double, std::less, pos_interval> pos_set;
int overlap_graph::tile_greedy(const std::vector<int>& sort_array,
                               const std::vector<const mega_read_info*>& mega_reads,
                               const std::vector<node_info>& nodes,
                               std::vector<int>& res,
                               size_t at_most) const {
  pos_set                   covered;
  std::vector<pos_interval> placed;
  int                       score = 0;

  for(const int it_i : sort_array) {
    const auto& mr     = *mega_reads[it_i];
    pos_interval pos_i(mr.tiling_start, mr.tiling_end);
    const double max_overlap       = std::min(k_len * overlap_play, boost::icl::length(pos_i));
    const auto   overlaps          = covered & pos_i;
    const bool   has_large_overlap =
      std::any_of(overlaps.begin(), overlaps.end(),
                  [=](const pos_interval& x) { return boost::icl::length(x) >= max_overlap; });

    if(has_large_overlap) continue;
    const bool contains =
      std::any_of(placed.begin(), placed.end(),
                  [&](const pos_interval& x) { return boost::icl::contains(pos_i, x); });
    if(contains) continue;

    covered += pos_i;
    placed.push_back(pos_i);
    score   += nodes[it_i].lpath;
    res.push_back(it_i);
    if(res.size() >= at_most) break;
  }
  return score;
  return 0;
}

struct max_tile_info {
  int    score;
  double pos;                 // last base position
  int    node;                // last mega read
  int    previous;            // back pointer
  int    length;              // number of mega reads in tiling
};
std::ostream& operator<<(std::ostream& os, const max_tile_info& i) {
  return os << "{score:" << i.score << ", pos:" << i.pos << ", node:" << i.node << ", prev:" << i.previous << ", len:" << i.length << "}";
}

//static max_tile_info mtig = { 0, std::numeric_limits<double>::min(), -1, -1, 0 };

int overlap_graph::tile_maximal(const std::vector<int>& sort_array,
                                const std::vector<const mega_read_info*>& mega_reads,
                                const std::vector<node_info>& nodes,
                                std::vector<int>& res) const {
  std::vector<max_tile_info> info;
  info.reserve(sort_array.size());

  auto       it  = sort_array.cbegin();
  const auto end = sort_array.cend();
  if(it == end) return 0;
  info.push_back({nodes[mega_reads[*it]->end_node].lpath, mega_reads[*it]->tiling_end, *it, -1, 1 });

  for(++it; it != end; ++it) {
    const double lpath_start = mega_reads[*it]->tiling_start;
    const auto lb = std::upper_bound(info.cbegin(), info.cend(),
                                     std::min(lpath_start + k_len * overlap_play, mega_reads[*it]->tiling_end),
                                     [](const double x, const max_tile_info& y) { return x < y.pos; });
    int i = std::distance(info.cbegin(), lb) - 1;
    while(i >= 0 && mega_reads[info[i].node]->tiling_start >= lpath_start)
      i = info[i].previous;

    const int nscore = (i >= 0 ? info[i].score : 0) + nodes[mega_reads[*it]->end_node].lpath;
    if(nscore > info.back().score)
      info.push_back({ nscore, mega_reads[*it]->tiling_end, *it, i, (i >= 0 ? info[i].length : 0) + 1 });
  }

  res.resize(info.back().length);
  int ptr = info.size() - 1;
  for(auto it = res.rbegin(); it != res.rend(); ++it) {
    assert(ptr >= 0);
    *it = info[ptr].node;
    int optr = ptr;
    ptr = info[ptr].previous;
    assert(ptr < optr);
  }
  assert(ptr < 0);
  return info.back().score;
  return 0;
}

void overlap_graph::print_mega_reads(std::ostream& output, const std::vector<int>& sort_array,
                                     const std::vector<const mega_read_info*>& mega_reads,
                                     const align_pb::coords_info_type& coords,
                                     const std::vector<node_info>& nodes,
                                     const std::vector<std::string>* unitigs_sequences,
                                     std::ostream* dot) const {
  for(const int cmr : sort_array) {
    const auto& mr      = *mega_reads[cmr];
    const auto& end_n   = nodes[mr.end_node];
    //    const auto& start_n = nodes[mr.start_node];
    const auto& end_c   = coords[mr.end_node];
    const auto& start_c = coords[mr.start_node];

    super_read_name        sr(mr.nb_unitigs);
    size_t                 offset = sr.prepend(end_c.name_u->unitigs, mr.end_unitig);
    int                    node_j = mr.end_node;
    int                    node_i = end_n.lprev;
    while(node_i >= 0) {
      const size_t overlap = nodes[node_i].lunitigs + coords[node_j].name_u->unitigs.size() - nodes[node_j].lunitigs;
      offset = sr.prepend(coords[node_i].name_u->unitigs, overlap, offset);
      if(dot) *dot << "n" << node_i << " -> n" << node_j << " [color=\"red\"];\n";
      node_j = node_i;
      node_i = nodes[node_i].lprev;
    }

    int sr_len = 0;
    for(auto unitig : sr.unitigs())
      sr_len += unitigs_lengths[unitig.id()];
    sr_len -= (sr.size() - 1) * (k_len - 1);

    // TODO: should the implied start and end of the mega read be
    // computed and kept in the mega_read_info, instead of being
    // computed when printing out results?
    const std::vector<std::string>& us = *unitigs_sequences;
    double imp_s = start_c.qs;
    if(mr.start_unitig) {
      for(int i = 0; i < mr.start_unitig; ++i)
        imp_s += us[start_c.name_u->unitigs[i].id_].size();
      imp_s -= (k_len - 1) * (mr.start_unitig - 1);
    }
    imp_s = imp_s * start_c.stretch + start_c.offset;

    double imp_e = end_c.qe;
    if(mr.end_unitig != (int)end_c.name_u->unitigs.size()) {
      for(int i = end_c.name_u->unitigs.size() - 1; i >= 0; --i)
        imp_e -= us[end_c.name_u->unitigs[i].id_].size();
      imp_e += (k_len - 1) * (end_c.name_u->unitigs.size() - mr.end_unitig - 1);
    }
    imp_e = imp_e * end_c.stretch + end_c.offset;

    output << std::fixed << std::setprecision(2)
           << imp_s << ' ' << imp_e << ' '
           << start_c.rs << ' ' << end_c.re << ' '
           << start_c.qs << ' ' << (sr_len - (end_c.ql - end_c.qe)) << ' '
           << end_n.lpath << ' ' << std::setprecision(4) << mr.density
           << ' ' << sr << ' ' << sr_len;

    if(unitigs_sequences) {
      output << ' ';
      sr.print_sequence(output, *unitigs_sequences, k_len);
    }

    output << '\n';
  }
}
