#ifndef _OVERLAP_GRAPH_H_
#define _OVERLAP_GRAPH_H_

#include <map>
#include <src_jf_aligner/union_find.hpp>
#include <src_jf_aligner/pb_aligner.hpp>

struct node_info {
  bool            start_node;
  bool            end_node;
  double          imp_s, imp_e; // coordinates used for longest path
  union_find::set component;
  int             lstart; // id of node starting longest path (-1 if first node in path)
  int             lprev;  // id of previous node in longest path (-1 if first)
  int             lpath;  // number of k-mers/bases in longest path
  int             lunitigs; // number of unitigs in longest path
  //  double          ldensity; // density of longest path

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
    lunitigs   = coords.name_u->unitigs.size();
  }

  node_info& l_start_node(std::vector<node_info>& nodes) {
    return lstart == -1 ? *this : nodes[lstart];
  }
  const node_info& l_start_node(const std::vector<node_info>& nodes) const {
    return lstart == -1 ? *this : nodes[lstart];
  }
};

// A mega-read is the index of two nodes, where the start node is on
// the backward longest path from the end node. In addition, we keep
// the offset in the unitig list of the start (resp. end) unitig to
// use.
struct mega_read_info {
  int    start_node, end_node;
  int    start_unitig, end_unitig;
  int    start_offset, end_offset;
  int    nb_unitigs;
  double imp_s, imp_e;
  double tiling_start, tiling_end; // coordinates used for tiling
  double density;
  static mega_read_info make(int i, const std::vector<node_info>& nodes,
                             const align_pb::coords_info_type&  coords);
};


struct overlap_graph {
  const double            overlap_play;
  const unsigned int      k_len;
  const std::vector<int>& unitigs_lengths;
  const double            nb_errors;
  const bool              maximize_bases;

  enum trim_action { NONE, MATCH, BRANCH };

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
  // component find the candidate mega-reads in each connected component.
  typedef std::map<union_find::set*, mega_read_info> comp_to_path;
  void mega_reads_per_comp(const int n, size_t pb_size, std::vector<node_info>& nodes,
                           const align_pb::coords_info_type& coords, comp_to_path& res,
                           double min_density = 0, double min_len = 0,
                           trim_action trim = NONE,
                          std::ostream* dot = nullptr) const;
  comp_to_path mega_reads_per_comp(const int n, const size_t pb_size, std::vector<node_info>& nodes,
                                   const align_pb::coords_info_type& coords,
                                   double min_density = 0, double min_len = 0,
                                   trim_action trim = NONE,
                                   std::ostream* dot = nullptr) const {
    comp_to_path res;
    mega_reads_per_comp(n, pb_size, nodes, coords, res, min_density, min_len, trim, dot);
    return res;
  }

  void trim_match(mega_read_info& mr, std::vector<node_info>& nodes,
                  const align_pb::coords_info_type& coords) const;

  // Tile the mega reads in a greedy fashion. Expect sort_array to be
  // the indices of the last nodes in the mega reads, sorted by size
  // of the mega reads (e.g. number of aligned k-mers or number of
  // bases).
  int tile_greedy(const std::vector<int>& sort_array,
                  const std::vector<const mega_read_info*>& mega_reads,
                  const std::vector<node_info>& nodes,
                  std::vector<int>& res, size_t at_most = std::numeric_limits<size_t>::max()) const;

  std::pair<int, std::vector<int>> tile_greedy(const std::vector<int>& sort_array,
                                               const std::vector<const mega_read_info*>& mega_reads,
                                               const std::vector<node_info>& nodes,
                                               size_t at_most = std::numeric_limits<size_t>::max()) const {
    std::vector<int> res;
    int score = tile_greedy(sort_array, mega_reads, nodes, res, at_most);
    return std::make_pair(score, std::move(res));
  }

  // Tile the mega reads maximizing the number of aligned
  // k-mers. Expect sort_array to be the indices of the last nodes in
  // the mega reads, sorted by the number of mers aligned.
  int tile_maximal(const std::vector<int>& sort_array,
                   const std::vector<const mega_read_info*>& mega_reads,
                   const std::vector<node_info>& nodes,
                   std::vector<int>& res) const;

  std::pair<int, std::vector<int>> tile_maximal(const std::vector<int>& sort_array,
                                                const std::vector<const mega_read_info*>& mega_reads,
                                                const std::vector<node_info>& nodes) const {
    std::vector<int> res;
    int score = tile_maximal(sort_array, mega_reads, nodes, res);
    return std::make_pair(score, std::move(res));
  }

  // Given the terminal nodes found by term_node_per_comp, print the mega reads
  void print_mega_reads(std::ostream& os, const std::vector<int>& sort_array,
                        const std::vector<const mega_read_info*>& mega_reads,
                        const align_pb::coords_info_type& coords,
                        const std::vector<node_info>& nodes,
                        const std::vector<std::string>* unitigs_sequences,
                        std::ostream* dot = 0) const;


  struct thread {
    const overlap_graph&               og_;
    std::vector<int>                   sort_nodes_, sort_tiling_, tiled_mr_;
    std::vector<node_info>             nodes_;
    const align_pb::coords_info_type*  coords_;
    comp_to_path                       comp_mega_reads_;
    std::vector<const mega_read_info*> mega_reads_;
    std::ostream*                      dot_;
    overlap_graph::trim_action         trim_;

    thread(const overlap_graph& og) : og_(og), dot_(nullptr), trim_(overlap_graph::NONE) { }

    void dot(std::ostream* d) { dot_ = d; }
    void trim_match() {trim_ = overlap_graph::MATCH; }

    void reset(const align_pb::coords_info_type& coords, const std::string& pb_name) {
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
      if(dot_) {
        *dot_ << "digraph \"" << pb_name << "\" {\nnode [fontsize=\"10\"];\n";
        for(size_t i = 0; i < sort_nodes_.size(); ++i) {
          const size_t it_i = sort_nodes_[i];
          *dot_ << "n" << it_i << "[tooltip=\"" << coords[it_i].name_u->unitigs << "\"];\n";
        }
      }
    }

    void traverse() { og_.traverse(sort_nodes_, *coords_, nodes_, dot_); }
    void term_node_per_comp(size_t pb_size, double min_density = 0, double min_len = 0) {
      comp_mega_reads_.clear();
      og_.mega_reads_per_comp(coords_->size(), pb_size, nodes_, *coords_, comp_mega_reads_, min_density, min_len,
                              trim_, dot_);
      mega_reads_.clear();
      sort_tiling_.clear();
      for(const auto& comp : comp_mega_reads_) {
        sort_tiling_.push_back(mega_reads_.size());
        mega_reads_.push_back(&comp.second);
      }
    }

    void tile_greedy(size_t at_most = std::numeric_limits<size_t>::max()) {
      std::sort(sort_tiling_.begin(), sort_tiling_.end(),
                [&](int i, int j) { return nodes_[mega_reads_[j]->end_node].lpath < nodes_[mega_reads_[i]->end_node].lpath; });
      tiled_mr_.clear();
      og_.tile_greedy(sort_tiling_, mega_reads_, nodes_, tiled_mr_, at_most);
      std::sort(tiled_mr_.begin(), tiled_mr_.end(),
                [&](int i, int j) -> bool {
                  const auto st_i = mega_reads_[i]->imp_s, st_j = mega_reads_[j]->imp_s;
                  return st_i < st_j || (st_i == st_j && mega_reads_[i]->imp_e < mega_reads_[j]->imp_e);
                });
    }

    void tile_maximal() {
      std::sort(sort_tiling_.begin(), sort_tiling_.end(),
                [&](int i, int j) { return mega_reads_[i]->tiling_end < mega_reads_[j]->tiling_end; });
      og_.tile_maximal(sort_tiling_, mega_reads_, nodes_, tiled_mr_);
      std::sort(tiled_mr_.begin(), tiled_mr_.end(),
                [&](int i, int j) -> bool {
                  const auto st_i = mega_reads_[i]->imp_s, st_j = mega_reads_[j]->imp_s;
                  return st_i < st_j || (st_i == st_j && mega_reads_[i]->imp_e < mega_reads_[j]->imp_e);
                });
    }


    void print_mega_reads(std::ostream& os, const std::string& name,
                          const std::vector<std::string>* unitigs_sequences = 0) const {
      if(!comp_mega_reads_.empty()) {
        //        os << '>' << name << ' ' << '(' << sort_mega_reads_.size() << ',' << tiling_.size() << ')' << '\n';
        os << '>' << name << '\n';
        og_.print_mega_reads(os, tiled_mr_.empty() ? sort_tiling_ : tiled_mr_, mega_reads_, *coords_, nodes_, unitigs_sequences, dot_);
        if(dot_)
          *dot_ << "}\n";
      }
    }
  };
};

#endif /* _OVERLAP_GRAPH_H_ */
