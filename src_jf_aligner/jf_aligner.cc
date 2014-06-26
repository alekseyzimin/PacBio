#include <memory>
#include <ios>
#include <thread>

#include <src_jf_aligner/jf_aligner_cmdline.hpp>
#include <src_jf_aligner/superread_parser.hpp>
#include <src_jf_aligner/pb_aligner.hpp>

jf_aligner_cmdline args;

/**
 * A unique ptr to a multiplexed output stream with a convenient
 * constructor.
 */
class mstream : public std::unique_ptr<Multiplexer::ostream> {
public:
  mstream(Multiplexer* m) :
    std::unique_ptr<Multiplexer::ostream>(m ? new Multiplexer::ostream(m) : 0)
  { }
};

void print_alignments(read_parser* reads, Multiplexer* details_m, Multiplexer* coords_m, const align_pb* align_data) {
  mer_dna                  tmp_m;
  parse_sequence           parser;
  align_pb::thread         aligner(*align_data);

  mstream     details_io(details_m);
  mstream     coords_io(coords_m);
  std::string name;

  while(true) {
    read_parser::job job(*reads);
    if(job.is_empty()) break;

    for(size_t i = 0; i < job->nb_filled; ++i) { // Process each read
      auto name_end = job->data[i].header.find_first_of(" \t\n\v\f\r");
      name = job->data[i].header.substr(0, name_end);
      const size_t pb_size = job->data[i].seq.size();
      parser.reset(job->data[i].seq);
      aligner.align_sequence_max(parser, pb_size);

      align_pb::print_coords(*coords_io, name, pb_size, args.compact_flag, aligner.coords());
      if(details_io) align_pb::print_details(*details_io, name, aligner.frags_pos());
    }
  }

}

int main(int argc, char *argv[])
{
  args.parse(argc, argv);
  mer_dna::k(args.mer_arg);
  std::ios::sync_with_stdio(false);

  if(!args.details_given && !args.coords_given)
    jf_aligner_cmdline::error() << "No output file given. Doing nothing ungracefully.";

  // Open output file for early error reporting
  output_file details;
  output_file coords;
  if(args.details_given) details.open(args.details_arg, args.threads_arg);
  if(args.coords_given) {
    coords.open(args.coords_arg, args.threads_arg);
  } else {
    coords.set(std::cout, args.threads_arg);
  }

  // Read k-unitig lengths
  std::unique_ptr<std::vector<int> > unitigs_lengths;
  if(args.unitigs_lengths_given) {
    if(!args.k_mer_given)
      jf_aligner_cmdline::error() << "The mer length used for generating the k-unitigs (-k, --k-mer) is required if the unitig lengths switch (-l, --unitig-lengths) is passed.";
    unitigs_lengths.reset(new std::vector<int>);
    std::ifstream is(args.unitigs_lengths_arg);
    if(!is.good())
      jf_aligner_cmdline::error() << "Failed to open unitig lengths map file '" << args.unitigs_lengths_arg << "'";
    std::string unitig;
    unsigned int len;
    while(is.good()) {
      is >> unitig >> len;
      unitigs_lengths->push_back(len);
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
                      args.forward_flag, args.max_match_flag, args.max_count_arg,
                      args.mers_matching_arg / 100.0, args.bases_matching_arg / 100.0);
  if(unitigs_lengths) align_data.unitigs_lengths(unitigs_lengths.get(), args.k_mer_arg);

  // Output header if necessary
  if(!args.no_header_flag)
    align_pb::print_coords_header(coords.multiplexer(), args.compact_flag);

  // Output alignements
  std::vector<std::thread> threads;
  for(unsigned int i = 0; i < args.threads_arg; ++i)
    threads.push_back(std::thread(print_alignments, &reads, details.multiplexer(), coords.multiplexer(),
                                  &align_data));
  for(unsigned int i = 0; i < args.threads_arg; ++i)
    threads[i].join();

  return 0;
}
