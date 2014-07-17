#include <memory>
#include <ios>
#include <thread>

#include <src_jf_aligner/create_mega_reads_cmdline.hpp>
#include <src_jf_aligner/superread_parser.hpp>
#include <src_jf_aligner/pb_aligner.hpp>
#include <src_jf_aligner/overlap_graph.hpp>

create_mega_reads_cmdline args;

void create_mega_reads(read_parser* reads, Multiplexer* output_m, const align_pb* align_data,
                       overlap_graph* graph_walker, Multiplexer* dot_m) {
  parse_sequence                        parser;
  align_pb::thread                      aligner(*align_data);
  overlap_graph::thread                 graph(*graph_walker);
  std::string                           name;
  Multiplexer::ostream                  output(output_m);
  std::unique_ptr<Multiplexer::ostream> dot(dot_m ? new Multiplexer::ostream(dot_m) : 0);

  while(true) {
    read_parser::job job(*reads);
    if(job.is_empty()) break;

    for(size_t i = 0; i < job->nb_filled; ++i) { // Process each PB read
      auto name_end = job->data[i].header.find_first_of(" \t\n\v\f\r");
      name = job->data[i].header.substr(0, name_end);
      const size_t pb_size = job->data[i].seq.size();
      parser.reset(job->data[i].seq);
      aligner.align_sequence_max(parser, pb_size);

      const auto& coords = aligner.coords();
      graph.reset(coords, name, dot.get());
      graph.traverse();
      graph.term_node_per_comp(pb_size);
      output << ">" << name << "\n";
      graph.print_mega_reads(output);
      output.end_record();
      if(dot)
        dot->end_record();
    }
  }

}

int main(int argc, char *argv[])
{
  args.parse(argc, argv);
  mer_dna::k(args.mer_arg);
  std::ios::sync_with_stdio(false);

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
  {
    std::ifstream is(args.unitigs_lengths_arg);
    if(!is.good())
      create_mega_reads_cmdline::error() << "Failed to open unitig lengths map file '" << args.unitigs_lengths_arg << "'";
    std::string unitig;
    unsigned int len;
    while(is.good()) {
      is >> unitig >> len;
      unitigs_lengths.push_back(len);
    }
  }

  // Read the super reads
  mer_pos_hash_type hash(args.size_arg);
  frag_lists names(args.threads_arg);
  superread_parse(args.threads_arg, hash, names, args.superreads_arg.cbegin(), args.superreads_arg.cend());

  // Prepare I/O
  stream_manager streams(args.pacbio_arg.cbegin(), args.pacbio_arg.cend());
  read_parser    reads(4 * args.threads_arg, 100, 1, streams);

  // Create aligner
  align_pb align_data(hash, args.stretch_constant_arg, args.stretch_factor_arg,
                      true /* forward */, args.max_match_flag,
                      args.max_count_arg ? args.max_count_arg : std::numeric_limits<int>::max(),
                      args.mers_matching_arg / 100.0, args.bases_matching_arg / 100.0);
  align_data.unitigs_lengths(&unitigs_lengths, args.k_mer_arg);

  // Output candidate mega_reads
  overlap_graph graph_walker(args.overlap_play_arg, args.k_mer_arg, unitigs_lengths, args.errors_arg, args.bases_flag);
  std::vector<std::thread> threads;
  for(unsigned int i = 0; i < args.threads_arg; ++i)
    threads.push_back(std::thread(create_mega_reads, &reads, output.multiplexer(), &align_data,
                                  &graph_walker, dot.multiplexer()));
  for(unsigned int i = 0; i < args.threads_arg; ++i)
    threads[i].join();

  return 0;
}
