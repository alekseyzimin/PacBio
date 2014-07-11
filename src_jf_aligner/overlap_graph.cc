#include <src_jf_aligner/overlap_graph.hpp>

///////////////////////
// Helper functions. //
///////////////////////

// Return the starting node of the longest path ending at node n
static node_info& l_node_start(node_info& n, std::vector<node_info>& nodes) {
  return n.lstart == -1 ? n : nodes[n.lstart];
}

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
         (nlpath == node_j.lpath && (node_j.lstart == -1 || l_node_start(node_i, nodes).imp_s > l_node_start(node_j, nodes).imp_s))) {
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

void overlap_graph::comp_mega_reads(const int n, size_t pb_size, std::vector<node_info>& nodes,
                                    const align_pb::coords_info_type& coords, comp_to_path& components,
                                    std::ostream* dot) const {
  // For each connected component, keep the index of the terminal node
  // of the longest path found
  for(int i = 0; i < n; ++i) {
    auto& node = nodes[i];
    const double imp_len = std::min((double)pb_size, node.imp_e) - std::max(1.0, l_node_start(node, nodes).imp_s);
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
           << "\\nP(" << ci.rs << ',' << ci.re << ")\\nS(" << ci.qs << ',' << ci.qe << ")"
           << "\\nLP #" << node.lpath << " L" << std::setprecision(1) << imp_len
           << " d" << std::setprecision(2) << node.ldensity << "\""
           << color << "];\n";
    }
    if(!node.end_node) continue;

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

void overlap_graph::print_mega_reads(std::ostream& output, const comp_to_path& mega_reads,
                                     const align_pb::coords_info_type& coords,
                                     const std::vector<node_info>& nodes,
                                     std::ostream* dot) const {
  for(const auto comp : mega_reads) {
    int         node_j  = comp.second;
    const auto& end_n   = nodes[node_j];
    const auto& start_n = end_n.lstart == -1 ? end_n : nodes[end_n.lstart];
    output << std::fixed << std::setprecision(2)
           << start_n.imp_s << ' ' << end_n.imp_e << ' '
           << end_n.lpath << ' ' << std::setprecision(4) << end_n.ldensity;
    const super_read_name *asr;
    super_read_name sr(end_n.lunitigs);
    if(end_n.lstart == -1) {
      asr = &coords[node_j].unitigs;
    } else {
      size_t          offset = sr.prepend(coords[node_j].unitigs);
      int             node_i = end_n.lprev;
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
    output << ' ' << *asr << ' ' << sr_len << '\n';
  }
}
