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

void create_mega_reads(coords_parser* parser, Multiplexer* output_m, overlap_graph* graph_walker, Multiplexer* dot_m) {
  align_pb::coords_info_type            coords;
  overlap_graph::thread                 graph(*graph_walker);
  Multiplexer::ostream                  output(output_m);
  std::unique_ptr<Multiplexer::ostream> dot(dot_m ? new Multiplexer::ostream(dot_m) : 0);

  for(coords_parser::stream coords_stream(*parser); coords_stream; ++coords_stream) {
    fill_coords(coords_stream->lines, coords);

    graph.reset(coords, coords_stream->header, dot.get());
    graph.traverse();
    graph.term_node_per_comp(coords[0].rl, args.density_arg, args.min_length_arg);
    graph.print_mega_reads(output, coords_stream->header);
    output.end_record();
    if(dot)
      dot->end_record();
  }
}

inline static std::istream& skip_header(std::istream& is) {
  return is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

int main(int argc, char* argv[]) {
  args.parse(argc, argv);
  std::ios::sync_with_stdio(false);

  if(!args.unitigs_lengths_given && !args.unitigs_sequences_given)
    cmdline_args::error() << "One of --unitigs-lengths or --unitigs-sequences is required.";

  // Open output file for early error reporting
  output_file output;
  if(args.output_given) {
    output.open(args.output_arg, args.threads_arg);
  } else {
    output.set(std::cout, args.threads_arg);
  }
  output_file dot;
  if(args.dot_given)
    dot.open(args.dot_arg, args.threads_arg);

  // Read k-unitig lengths
  std::vector<int> unitigs_lengths;
  std::vector<std::string> sequences;
  {
    if(args.unitigs_lengths_given) { // File with lengths
      std::ifstream is(args.unitigs_lengths_arg);
      if(!is.good())
        cmdline_args::error() << "Failed to open unitig lengths map file '" << args.unitigs_lengths_arg << "'";
      std::string unitig;
      unsigned int len;
      is >> unitig >> len;
      while(is.good()) {
        unitigs_lengths.push_back(len);
        is >> unitig >> len;
      }
    } else { // Sequence in fasta file given
      std::ifstream is(args.unitigs_sequences_arg);
      if(!is.good())
        cmdline_args::error() << "Failed to open unitigs sequence file '" << args.unitigs_sequences_arg << "'";
      while(skip_header(is)) {
        sequences.push_back("");
        std::getline(is, sequences.back());
        unitigs_lengths.push_back(sequences.back().size());
      }
    }
  }

  coords_parser parser(args.threads_arg, args.coords_arg);
  parser.start_parsing();

  overlap_graph graph_walker(args.overlap_play_arg, args.k_mer_arg, unitigs_lengths, args.errors_arg, args.bases_flag);
  std::vector<std::thread> threads;
  for(unsigned int i = 0; i < args.threads_arg; ++i)
    threads.push_back(std::thread(create_mega_reads, &parser, output.multiplexer(), &graph_walker, dot.multiplexer()));
  for(auto& th : threads)
    th.join();

  return 0;
}
