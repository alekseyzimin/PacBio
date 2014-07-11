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

  void traverse(const std::vector<int>& sort_array, const align_pb::coords_info_type& coords,
                std::vector<node_info>& nodes, std::ostream* dot = 0) const;
  typedef std::map<union_find::set*, int> comp_to_path;
  void comp_mega_reads(const int n, size_t pb_size, std::vector<node_info>& nodes,
                       const align_pb::coords_info_type& coords, comp_to_path& res,
                       std::ostream* dot = 0) const;
  comp_to_path comp_mega_reads(const int n, const size_t pb_size, std::vector<node_info>& nodes,
                               const align_pb::coords_info_type& coords,
                               std::ostream* dot = 0) const {
    comp_to_path res;
    comp_mega_reads(n, pb_size, nodes, coords, res);
    return res;
  }

  void print_mega_reads(std::ostream& os, const comp_to_path& mega_reads,
                        const align_pb::coords_info_type& coords,
                        const std::vector<node_info>& nodes,
                        std::ostream* dot = 0) const;


  struct thread {
    const overlap_graph&              og_;
    std::vector<int>                  sort_array_;
    std::vector<node_info>            nodes_;
    const align_pb::coords_info_type* coords_;
    comp_to_path                      mega_reads_;
    std::ostream*                     dot_;

    thread(const overlap_graph& og) : og_(og) { }

    void reset(const align_pb::coords_info_type& coords, std::ostream* dot = 0) {
      coords_     = &coords;
      const int n = coords_->size();
      sort_array_.resize(n);
      if((int)nodes_.size() < n)
        nodes_.resize(n);
      for(int i = 0; i < n; ++i) {
        sort_array_[i] = i;
        nodes_[i].reset(coords[i], og_.maximize_bases);
      }
      std::sort(sort_array_.begin(), sort_array_.end(),
                [&] (int i, int j) { return nodes_[i].imp_s < nodes_[j].imp_s || (nodes_[i].imp_s == nodes_[j].imp_s &&
                                                                                  nodes_[i].imp_e < nodes_[j].imp_e); });
      dot_ = dot;
      if(dot_) {
        for(size_t i = 0; i < sort_array_.size(); ++i) {
          const size_t it_i = sort_array_[i];
          *dot_ << "n" << it_i << "[tooltip=\"" << coords[it_i].qname << "\"];\n";
        }
      }
    }

    void traverse() { og_.traverse(sort_array_, *coords_, nodes_, dot_); }
    void compute_mega_reads(size_t pb_size) {
      mega_reads_.clear();
      og_.comp_mega_reads(coords_->size(), pb_size, nodes_, *coords_, mega_reads_, dot_);
    }
    void print_mega_reads(std::ostream& os) const {
      og_.print_mega_reads(os, mega_reads_, *coords_, nodes_, dot_);
    }
  };
};

#endif /* _OVERLAP_GRAPH_H_ */
