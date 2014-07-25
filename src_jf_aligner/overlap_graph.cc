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
      const int nb_u_overlap = coords_i.unitigs.overlap(coords_j.unitigs);
      if(!nb_u_overlap) continue; // No overlap according to unitig names
      if(coords_i.unitigs == coords_j.unitigs) continue;
      int u_overlap_len  = 0;
      int common_overlap = 0;
      for(int u = 0; u < nb_u_overlap; ++u) {
        u_overlap_len    += unitigs_lengths[coords_j.unitigs.unitig_id(u)];
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
        node_j.lunitigs = node_i.lunitigs + coords_j.unitigs.size() - nb_u_overlap;
      }
      if(dot)
        *dot << "n" << it_i << " -> n" << it_j << " [tooltip=\"" << "..." << "\", label=\"" << common_overlap << "\"];\n";
    }
  }
}

void overlap_graph::term_node_per_comp(const int n, size_t pb_size, std::vector<node_info>& nodes,
                                       const align_pb::coords_info_type& coords, comp_to_path& components,
                                       double min_density, double min_len, std::ostream* dot) const {
  // For each connected component, keep the index of the terminal node
  // of the longest path found
  for(int i = 0; i < n; ++i) {
    auto& node = nodes[i];
    const double imp_len = std::min((double)pb_size + 0.5, node.imp_e) - std::max(0.5, node.l_start_node(nodes).imp_s);
    node.ldensity        = (double)node.lpath / imp_len;
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
           << " d" << std::setprecision(2) << node.ldensity << "\""
           << color << "];\n";
    }
    if(!node.end_node || node.ldensity < min_density || imp_len < min_len) continue;

    auto comp_root = node.component.root();
    auto comp_it   = components.find(comp_root);
    if(comp_it == components.end()) {
      components.insert(std::make_pair(comp_root, i));
      continue;
    }
    const auto& onode = nodes[comp_it->second]; // current terminal node of longest path
    if(node.lpath > onode.lpath || (node.lpath == onode.lpath &&  node.ldensity > onode.ldensity))
      comp_it->second = i;
  }
}

typedef boost::icl::right_open_interval<double> pos_interval;
typedef boost::icl::interval_set<double, std::less, pos_interval> pos_set;
int overlap_graph::tile_greedy(const std::vector<int>& sort_array,
                               const std::vector<node_info>& nodes, std::vector<int>& res,
                               size_t at_most) const {
  pos_set                   covered;
  std::vector<pos_interval> placed;
  int                       score = 0;

  for(const int it_i : sort_array) {
    const auto& node_i = nodes[it_i];
    pos_interval pos_i(node_i.l_start_node(nodes).imp_s, node_i.imp_e);
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
}

struct max_tile_info {
  int    score;
  double pos;                 // last base position
  int    node;                // last node
  int    previous;            // back pointer
  int    length;              // number of mega reads in tiling
};
std::ostream& operator<<(std::ostream& os, const max_tile_info& i) {
  return os << "{score:" << i.score << ", pos:" << i.pos << ", node:" << i.node << ", prev:" << i.previous << ", len:" << i.length << "}";
}

static max_tile_info mtig = { 0, std::numeric_limits<double>::min(), -1, -1, 0 };

int overlap_graph::tile_maximal(const std::vector<int>& sort_array,
                                const std::vector<node_info>& nodes, std::vector<int>& res) const {
  std::vector<max_tile_info> info;
  info.reserve(sort_array.size());

  auto       it  = sort_array.cbegin();
  const auto end = sort_array.cend();
  if(it == end) return 0;
  info.push_back({nodes[*it].lpath, nodes[*it].imp_e, *it, -1, 1 });

  for(++it; it != end; ++it) {
    const double lpath_start = nodes[*it].l_start_node(nodes).imp_s;
    const auto lb = std::upper_bound(info.cbegin(), info.cend(),
                                     std::min(lpath_start + k_len * overlap_play, nodes[*it].imp_e),
                                     [](const double x, const max_tile_info& y) { return x < y.pos; });
    int i = std::distance(info.cbegin(), lb) - 1;
    while(i >= 0 && nodes[info[i].node].l_start_node(nodes).imp_s >= lpath_start)
      i = info[i].previous;

    const int nscore = (i >= 0 ? info[i].score : 0) + nodes[*it].lpath;
    if(nscore > info.back().score)
      info.push_back({ nscore, nodes[*it].imp_e, *it, i, (i >= 0 ? info[i].length : 0) + 1 });
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
}

void overlap_graph::print_mega_reads(std::ostream& output, const std::vector<int>& sort_array,
                                     const align_pb::coords_info_type& coords,
                                     const std::vector<node_info>& nodes,
                                     const std::vector<std::string>* unitigs_sequences,
                                     std::ostream* dot) const {
  for(const int cnode : sort_array) {
    const auto& end_n   = nodes[cnode];
    const auto& start_n = end_n.l_start_node(nodes);
    output << std::fixed << std::setprecision(2)
           << start_n.imp_s << ' ' << end_n.imp_e << ' '
           << end_n.lpath << ' ' << std::setprecision(4) << end_n.ldensity;
    const super_read_name *asr;
    super_read_name sr(end_n.lunitigs);
    if(end_n.lstart == -1) {
      asr = &coords[cnode].unitigs;
    } else {
      size_t offset = sr.prepend(coords[cnode].unitigs);
      int    node_j = cnode;
      int    node_i = end_n.lprev;
      while(node_i >= 0) {
        const size_t overlap = nodes[node_i].lunitigs + coords[node_j].unitigs.size() - nodes[node_j].lunitigs;
        offset = sr.prepend(coords[node_i].unitigs, overlap, offset);
        if(dot) *dot << "n" << node_i << " -> n" << node_j << " [color=\"red\"];\n";
        node_j = node_i;
        node_i = nodes[node_i].lprev;
      }
      asr = &sr;
    }

    int sr_len = 0;
    for(auto unitig : asr->unitigs())
      sr_len += unitigs_lengths[unitig.id()];
    sr_len -= (sr.size() - 1) * (k_len - 1);
    output << ' ' << *asr << ' ' << sr_len;

    if(unitigs_sequences) {
      output << ' ';
      asr->print_sequence(output, *unitigs_sequences, k_len);
    }

    output << '\n';
  }
}
