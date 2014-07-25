#ifndef _OVERLAP_GRAPH_H_
#define _OVERLAP_GRAPH_H_

#include <map>
#include <src_jf_aligner/union_find.hpp>
#include <src_jf_aligner/pb_aligner.hpp>

struct node_info {
  bool            start_node;
  bool            end_node;
  double          imp_s, imp_e;
  union_find::set component;
  int             lstart; // id of node starting longest path (-1 if first node in path)
  int             lprev;  // id of previous node in longest path (-1 if first)
  int             lpath;  // number of k-mers/bases in longest path
  int             lunitigs; // number of unitigs in longest path
  double          ldensity; // density of longest path (updated only at terminal node of longest paths)

  node_info() = default;

  void reset(const align_pb::coords_info& coords, bool bases) {
    start_node = true;
    end_node   = true;
    imp_s      = coords.stretch + coords.offset;
    imp_e      = coords.stretch * coords.ql + coords.offset;
    component.reset();
    lstart     = -1;
    lprev      = -1;
    lpath      = bases ? coords.sr_cover : coords.nb_mers;
    lunitigs   = coords.unitigs.size();
  }

  node_info& l_start_node(std::vector<node_info>& nodes) {
    return lstart == -1 ? *this : nodes[lstart];
  }
  const node_info& l_start_node(const std::vector<node_info>& nodes) const {
    return lstart == -1 ? *this : nodes[lstart];
  }

};

struct overlap_graph {
  const double            overlap_play;
  const unsigned int      k_len;
  const std::vector<int>& unitigs_lengths;
  const double            nb_errors;
  const bool              maximize_bases;

  // Used to traverse a graph of overlaps between super reads. They
  // must overlap according to their position (e.g. alignment to a PB
  // read) and according to the unitigs sequence.
  //
  // play: ratio of difference in length between position overlap and
  // unitig overlap. E.g. 1.2 allows for 20% of play.
  //
  // len: length of k-mer used for creating k-unitigs
  //
  // lengths: vector of unitigs lengths in bases
  //
  // errors: number of average errors from best linear fit to use for
  // computing position.
  //
  // bases: maximize number of bases in path instead of number of
  // mers.
  overlap_graph(double play, unsigned int len, const std::vector<int>& lengths, double errors,
                bool bases) :
    overlap_play(play), k_len(len), unitigs_lengths(lengths), nb_errors(errors),
    maximize_bases(bases)
  { }

  // Traverse the graph and mark the longest paths. At the same time,
  // computes the connected components.
  void traverse(const std::vector<int>& sort_array, const align_pb::coords_info_type& coords,
                std::vector<node_info>& nodes, std::ostream* dot = 0) const;

  // Given the information created in 'traverse', for each connected
  // component find the end node of the longest path.
  typedef std::map<union_find::set*, int> comp_to_path;
  void term_node_per_comp(const int n, size_t pb_size, std::vector<node_info>& nodes,
                          const align_pb::coords_info_type& coords, comp_to_path& res,
                          double min_density = 0, double min_len = 0,
                          std::ostream* dot = 0) const;
  comp_to_path term_node_per_comp(const int n, const size_t pb_size, std::vector<node_info>& nodes,
                                  const align_pb::coords_info_type& coords,
                                  double min_density = 0, double min_len = 0,
                                  std::ostream* dot = 0) const {
    comp_to_path res;
    term_node_per_comp(n, pb_size, nodes, coords, res, min_density, min_len, dot);
    return res;
  }

  // Tile the mega reads in a greedy fashion. Expect sort_array to be
  // the indices of the last nodes in the mega reads, sorted by size
  // of the mega reads (e.g. number of aligned k-mers or number of
  // bases).
  int tile_greedy(const std::vector<int>& sort_array,
                   const std::vector<node_info>& nodes,
                   std::vector<int>& res, size_t at_most = std::numeric_limits<size_t>::max()) const;

  std::pair<int, std::vector<int>> tile_greedy(const std::vector<int>& sort_array,
                                               const std::vector<node_info>& nodes,
                                               size_t at_most = std::numeric_limits<size_t>::max()) const {
    std::vector<int> res;
    int score = tile_greedy(sort_array, nodes, res, at_most);
    return std::make_pair(score, std::move(res));
  }

  // Tile the mega reads maximizing the number of aligned
  // k-mers. Expect sort_array to be the indices of the last nodes in
  // the mega reads, sorted by the number of mers aligned.
  int tile_maximal(const std::vector<int>& sort_array,
                    const std::vector<node_info>& nodes,
                    std::vector<int>& res) const;

  std::pair<int, std::vector<int>> tile_maximal(const std::vector<int>& sort_array,
                                                const std::vector<node_info>& nodes) const {
    std::vector<int> res;
    int score = tile_maximal(sort_array, nodes, res);
    return std::make_pair(score, std::move(res));
  }

  // Given the terminal nodes found by term_node_per_comp, print the mega reads
  void print_mega_reads(std::ostream& os, const std::vector<int>& sort_array,
                        const align_pb::coords_info_type& coords,
                        const std::vector<node_info>& nodes,
                        const std::vector<std::string>* unitigs_sequences,
                        std::ostream* dot = 0) const;


  struct thread {
    const overlap_graph&              og_;
    std::vector<int>                  sort_nodes_, sort_mega_reads_, tiling_;
    std::vector<node_info>            nodes_;
    const align_pb::coords_info_type* coords_;
    comp_to_path                      mega_reads_;
    std::ostream*                     dot_;

    thread(const overlap_graph& og) : og_(og) { }

    void reset(const align_pb::coords_info_type& coords, const std::string& pb_name, std::ostream* dot = 0) {
      coords_     = &coords;
      const int n = coords_->size();
      sort_nodes_.resize(n);
      nodes_.resize(n);
      for(int i = 0; i < n; ++i) {
        sort_nodes_[i] = i;
        nodes_[i].reset(coords[i], og_.maximize_bases);
      }
      std::sort(sort_nodes_.begin(), sort_nodes_.end(),
                [&] (int i, int j) { return nodes_[i].imp_s < nodes_[j].imp_s || (nodes_[i].imp_s == nodes_[j].imp_s &&
                                                                                  nodes_[i].imp_e < nodes_[j].imp_e); });
      dot_ = dot;
      if(dot_) {
        *dot << "digraph \"" << pb_name << "\" {\nnode [fontsize=\"10\"];\n";
        for(size_t i = 0; i < sort_nodes_.size(); ++i) {
          const size_t it_i = sort_nodes_[i];
          *dot << "n" << it_i << "[tooltip=\"" << coords[it_i].unitigs << "\"];\n";
        }
      }
    }

    void traverse() { og_.traverse(sort_nodes_, *coords_, nodes_, dot_); }
    void term_node_per_comp(size_t pb_size, double min_density = 0, double min_len = 0) {
      mega_reads_.clear();
      og_.term_node_per_comp(coords_->size(), pb_size, nodes_, *coords_, mega_reads_, min_density, min_len, dot_);
      sort_mega_reads_.clear();
      for(const auto& comp : mega_reads_)
        sort_mega_reads_.push_back(comp.second);
      tiling_.clear();
    }

    void tile_greedy(size_t at_most = std::numeric_limits<size_t>::max()) {
      std::sort(sort_mega_reads_.begin(), sort_mega_reads_.end(),
                [&](int i, int j) { return nodes_[j].lpath < nodes_[i].lpath; });
      og_.tile_greedy(sort_mega_reads_, nodes_, tiling_, at_most);
      std::sort(tiling_.begin(), tiling_.end(),
                [&](int i, int j) -> bool {
                  auto st_i = nodes_[i].l_start_node(nodes_).imp_s, st_j = nodes_[j].l_start_node(nodes_).imp_s;
                  return st_i < st_j || (st_i == st_j && nodes_[i].imp_e < nodes_[j].imp_e);
                });
    }

    void tile_maximal() {
      std::sort(sort_mega_reads_.begin(), sort_mega_reads_.end(),
                [&](int i, int j) { return nodes_[i].imp_e < nodes_[j].imp_e; });
      og_.tile_maximal(sort_mega_reads_, nodes_, tiling_);
      std::sort(tiling_.begin(), tiling_.end(),
                [&](int i, int j) -> bool {
                  const auto st_i = nodes_[i].l_start_node(nodes_).imp_s;
                  const auto st_j = nodes_[j].l_start_node(nodes_).imp_s;
                  return st_i < st_j || (st_i == st_j && nodes_[j].imp_e < nodes_[i].imp_e);
                });
    }


    void print_mega_reads(std::ostream& os, const std::string& name,
                          const std::vector<std::string>* unitigs_sequences = 0) const {
      if(!mega_reads_.empty()) {
        //        os << '>' << name << ' ' << '(' << sort_mega_reads_.size() << ',' << tiling_.size() << ')' << '\n';
        os << '>' << name << '\n';
        og_.print_mega_reads(os, tiling_.empty() ? sort_mega_reads_ : tiling_, *coords_, nodes_, unitigs_sequences, dot_);
        if(dot_)
          *dot_ << "}\n";
      }
    }
  };
};

#endif /* _OVERLAP_GRAPH_H_ */
