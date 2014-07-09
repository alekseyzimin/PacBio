#include <src_jf_aligner/overlap_graph.hpp>

///////////////////////
// Helper functions. //
///////////////////////

// Return the starting node of the longest path ending at node n
static node_info& l_node_start(node_info& n, std::vector<node_info>& nodes) {
  return n.lstart == -1 ? n : nodes[n.lstart];
}

static int implied_len(const node_info& n, long pb_size) {
  return std::min(pb_size, lrint(n.imp_e)) - std::max(lrint(n.imp_s), (long)1);
}

void overlap_graph::traverse(const std::vector<int>& sort_array, const align_pb::coords_info_type& coords,
                             std::vector<node_info>& nodes) const {
  for(size_t it_i = 0; it_i != sort_array.size(); ++it_i) {
    auto&       node_i   = nodes[it_i];
    const auto& coords_i = coords[it_i];
    if(node_i.imp_e >= coords_i.rl)
      continue; // Hanging off 3' end of PacBio read: no overlap
    for(size_t it_j = it_i + 1; it_j != sort_array.size(); ++it_j) {
      auto&       node_j   = nodes[it_j];
      const auto& coords_j = coords[it_j];
      if(node_j.imp_s <= 1)
        continue; // Hanging off 5' end of PacBio read: no overlap
      const int position_len = lrint(node_i.imp_e) - lrint(node_j.imp_s) + 1;
      if(position_len * overlap_play < k_len)
        break; // maximum implied overlap is less than a k-mer length
      const int nb_u_overlap = coords_i.unitigs.overlap(coords_j.unitigs);
      if(!nb_u_overlap)
        continue; // No overlap according to unitig names
      int u_overlap_len  = 0;
      int common_overlap = 0;
      for(int i = 0; i < nb_u_overlap; ++i) {
        u_overlap_len    += unitigs_lengths[coords_j.unitigs.unitig_id(i)];
        common_overlap   += maximize_bases ? coords_j.bases_info[2 * i] : coords_j.kmers_info[2 * i];
        if(i > 0)
          common_overlap -= maximize_bases ? coords_j.bases_info[2 * i - 1] : coords_j.kmers_info[2 * i - 1];
      }
      u_overlap_len -= (nb_u_overlap - 1) * k_len;
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
         (nlpath == node_j.lpath && l_node_start(node_i, nodes).imp_s > l_node_start(node_j, nodes).imp_s)) {
        node_j.lpath  = nlpath;
        node_j.lstart = node_i.lstart;
        node_j.lprev  = it_i;
        node_j.lunitigs = node_i.lunitigs + coords_j.unitigs.size() - nb_u_overlap;
      }
    }
  }
}

void overlap_graph::comp_mega_reads(const int n, size_t pb_size, std::vector<node_info>& nodes, comp_to_path& components) const {
  // For each connected component, keep the index of the terminal node
  // of the longest path found
  for(int i = 0; i < n; ++i) {
    auto& node = nodes[i];
    if(!node.end_node) continue;
    auto comp_root = node.component.root();
    auto comp_it   = components.find(comp_root);
    if(comp_it == components.end()) {
      components.insert(std::make_pair(comp_root, i));
      continue;
    }
    const auto& onode = nodes[comp_it->second]; // current terminal node of longest path
    const int imp_len = implied_len(node, pb_size);
    if(node.lpath > onode.lpath ||
       (node.lpath == onode.lpath &&  imp_len < implied_len(onode, pb_size))) {
      node.ldensity = (double)node.lpath / (double)imp_len;
      comp_it->second = i;
    }
  }
}

void overlap_graph::print_mega_reads(std::ostream& output, const comp_to_path& mega_reads,
                                     const align_pb::coords_info_type& coords,
                                     const std::vector<node_info>& nodes) const {
  for(const auto comp : mega_reads) {
    int         node_j  = comp.second;
    const auto& end_n   = nodes[node_j];
    const auto& start_n = end_n.lstart == -1 ? end_n : nodes[end_n.lstart];
    output << start_n.imp_s << " " << end_n.imp_e << " "
           << end_n.lpath << " " << end_n.ldensity;
    if(end_n.lstart == -1) {
      output << coords[node_j].unitigs.name() << "\n";
    } else {
      super_read_name sr(end_n.lunitigs);
      size_t          offset = sr.prepend(coords[node_j].unitigs);
      int             node_i = end_n.lprev;
      while(node_i >= 0) {
        const size_t overlap = nodes[node_i].lunitigs + coords[node_j].unitigs.size() - nodes[node_j].lunitigs;
        offset = sr.prepend(coords[node_i].unitigs, overlap, offset);
        node_i = nodes[node_i].lprev;
      }
    }
  }

}
