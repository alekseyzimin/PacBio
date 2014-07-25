#include <iostream>
#include <fstream>

#include <src_jf_aligner/tile_mega_reads_cmdline.hpp>
#include <src_jf_aligner/overlap_graph.hpp>

typedef tile_mega_reads_cmdline cmdline_args;
cmdline_args args;

void tile_mega_reads(std::ostream& out, const overlap_graph& og, const std::string& header, const std::vector<node_info>& nodes) {
  std::vector<int> sort(nodes.size());
  for(int i = 0, s = nodes.size(); i < s; ++i) sort[i] = i;
  std::vector<int> res;
  if(args.maximal_tiling_flag)
    og.tile_maximal(sort, nodes, res);
  else
    og.tile_greedy(sort, nodes, res);

  // out << header << '\n';
  for(auto i : res) {
    const auto& n = nodes[i];
    out << n.imp_s << ' ' << n.imp_e << ' ' << n.lpath << ' ' << n.ldensity
        << '\n';
  }
}

int main(int argc, char *argv[]) {
  args.parse(argc, argv);

  std::ifstream          input(args.mega_reads_arg);
  std::string            header;
  std::string            unitigs, sequence;
  long                   length;
  std::vector<node_info> nodes;
  std::vector<int>       uls;
  int id;

  overlap_graph og(args.overlap_play_arg, args.k_mer_arg, uls, 0, false);

  while(true) {
    // if(input.peek() == '>' || input.peek() == EOF) {
    //   if(!nodes.empty())
    //     tile_mega_reads(std::cout, og, header, nodes);
    //   std::getline(input, header);
    //   if(!input.good()) break;
    //   nodes.clear();
    //   continue;
    // }

    node_info n;
    n.start_node = n.end_node = true;
    n.lstart = n.lprev = -1;
    input >> id >> n.imp_s >> n.imp_e >> n.lpath >> n.ldensity >> unitigs >> length >> sequence;
    if(!input.good()) break;
    n.lunitigs = std::count(unitigs.cbegin(), unitigs.cend(), '_') + 1;
    nodes.push_back(std::move(n));
  }
  tile_mega_reads(std::cout, og, "", nodes);

  return 0;
}
