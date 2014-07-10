#include <cerrno>
#include <climits>
#include <ios>
#include <thread>

#include <src_jf_aligner/longest_path_overlap_graph2_cmdline.hpp>
#include <src_jf_aligner/output_file.hpp>
#include <src_jf_aligner/pb_aligner.hpp>
#include <src_jf_aligner/overlap_graph.hpp>
#include <src_jf_aligner/coords_parsing.hpp>


typedef longest_path_overlap_graph2_cmdline cmdline_args;
cmdline_args args;

void fill_coords(std::vector<std::string>& lines, align_pb::coords_info_type& coords) {
  coords.resize(lines.size());
  for(size_t i = 0; i < lines.size(); ++i) {
    std::istringstream is(lines[i]);
    is >> coords[i];
  }
}

void create_mega_reads(coords_parser* parser, Multiplexer* output_m, overlap_graph* graph_walker) {
  align_pb::coords_info_type coords;
  overlap_graph::thread      graph(*graph_walker);
  Multiplexer::ostream       output(output_m);

  for(coords_parser::stream coords_stream(*parser); coords_stream; ++coords_stream) {
    fill_coords(coords_stream->lines, coords);
    graph.reset(coords);
    graph.traverse();
    graph.compute_mega_reads(coords[0].rl);
    output << ">" << coords_stream->header << "\n";
    graph.print_mega_reads(output);
  }
}

int main(int argc, char* argv[]) {
  args.parse(argc, argv);
  std::ios::sync_with_stdio(false);

  // Open output file for early error reporting
  output_file output;
  if(args.output_given) {
    output.open(args.output_arg, args.threads_arg);
  } else {
    output.set(std::cout, args.threads_arg);
  }

  // Read k-unitig lengths
  std::vector<int> unitigs_lengths;
  {
    std::ifstream is(args.unitigs_lengths_arg);
    if(!is.good())
      cmdline_args::error() << "Failed to open unitig lengths map file '" << args.unitigs_lengths_arg << "'";
    std::string unitig;
    unsigned int len;
    while(is.good()) {
      is >> unitig >> len;
      unitigs_lengths.push_back(len);
    }
  }

  coords_parser parser(args.threads_arg, args.coords_arg);
  parser.start_parsing();

  overlap_graph graph_walker(args.overlap_play_arg, args.k_mer_arg, unitigs_lengths, args.errors_arg, args.bases_flag);
  std::vector<std::thread> threads;
  for(unsigned int i = 0; i < args.threads_arg; ++i)
    threads.push_back(std::thread(create_mega_reads, &parser, output.multiplexer(), &graph_walker));
  for(unsigned int i = 0; i < args.threads_arg; ++i)
    threads[i].join();


  return 0;
}
