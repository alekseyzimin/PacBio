#include <gtest/gtest.h>
#include <boost/icl/interval_set.hpp>
#include <src_jf_aligner/overlap_graph.hpp>
#include <random>

namespace {
std::vector<int> range(int n) {
  std::vector<int> res;
  res.reserve(n);
  for(int i = 0; i < n; ++i)
    res.push_back(i);
  return res;
}

typedef boost::icl::right_open_interval<double> pos_interval;
typedef boost::icl::interval_set<double, std::less, pos_interval> pos_set;
void check_no_overlap(std::vector<int> tiling, std::vector<node_info> nodes, double max_ovl) {
  pos_set covered;
  for(auto it : tiling) {
    pos_interval pos(nodes[it].l_start_node(nodes).imp_s, nodes[it].imp_e);
    const auto   overlaps           = covered & pos;
    covered                        += pos;
    const double max                = std::min(max_ovl, boost::icl::length(pos));
    const bool   has_large_overlap  = std::any_of(overlaps.begin(), overlaps.end(),
                                                  [=](const pos_interval& x) { return boost::icl::length(x) >= max; });
    EXPECT_FALSE(has_large_overlap);
  }
}

int compute_score(const std::vector<int> tiling, const std::vector<node_info>& nodes) {
  int score = 0;
  for(auto it : tiling) score += nodes[it].lpath;
  return score;
}

TEST(Tiling,Uniform) {
  static const int pb_len = 1000;
  static const int nb_nodes = 100;

  const unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine rng(seed1);
  std::uniform_real_distribution<double> pos_dist(0, pb_len);
  std::uniform_real_distribution<double> density_dist(0.02, 0.2);

  std::vector<node_info>             nodes(nb_nodes);
  std::vector<mega_read_info>        mega_reads_data(nb_nodes);
  std::vector<const mega_read_info*> mega_reads(nb_nodes);
  for(int i = 0; i < nb_nodes; ++i) {
    auto& n  = nodes[i];
    auto& mr = mega_reads_data[i];
    n.imp_s = pos_dist(rng);
    n.imp_e = pos_dist(rng);
    if(n.imp_e < n.imp_s) std::swap(n.imp_s, n.imp_e);
    mr.tiling_start = n.imp_s;
    mr.tiling_end   = n.imp_e;
    mr.density = density_dist(rng);
    n.lpath = mr.density * (n.imp_e - n.imp_s);
    n.lstart = -1;
    mr.start_node = mr.end_node = i;
    mega_reads[i] = &mr;
  }

  auto sort_greedy  = range(nb_nodes);
  auto sort_maximal = range(nb_nodes);
  std::sort(sort_greedy.begin(), sort_greedy.end(), [&](int i, int j) { return nodes[i].lpath < nodes[j].lpath; });
  std::sort(sort_maximal.begin(), sort_maximal.end(), [&](int i, int j) { return nodes[i].imp_e < nodes[j].imp_e; });

  overlap_graph graph(1.2, 70, {}, 0, false);
  const auto tile_greedy = graph.tile_greedy(sort_greedy, mega_reads, nodes);
  const auto tile_maximal = graph.tile_maximal(sort_maximal, mega_reads, nodes);
  EXPECT_LE(tile_greedy.first, tile_maximal.first);
  EXPECT_EQ(compute_score(tile_greedy.second, nodes), tile_greedy.first);
  EXPECT_EQ(compute_score(tile_maximal.second, nodes), tile_maximal.first);
  check_no_overlap(tile_greedy.second, nodes, graph.overlap_play * graph.k_len);
  check_no_overlap(tile_maximal.second, nodes, graph.overlap_play * graph.k_len);
}
}
