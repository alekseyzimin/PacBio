#include <memory>
#include <ios>

#include <src_jf_aligner/jf_aligner_cmdline.hpp>
#include <src_jf_aligner/superread_parser.hpp>
#include <src_jf_aligner/pb_aligner.hpp>

jf_aligner_cmdline args;

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
  if(args.coords_given) coords.open(args.coords_arg, args.threads_arg);

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
  superread_parse(args.threads_arg, hash, names, args.superreads_arg.cbegin(), args.superreads_arg.cend(),
                  args.compress_flag);

  // Create aligner
  stream_manager streams(args.pacbio_arg.cbegin(), args.pacbio_arg.cend());
  align_pb aligner(args.threads_arg, hash, streams, args.stretch_constant_arg, args.stretch_factor_arg,
                   args.consecutive_arg, args.nmers_arg, args.forward_flag, args.compress_flag, args.duplicated_flag);
  aligner
    .max_mer_count(args.max_count_arg)
    .compact_format(args.compact_flag);
  if(args.details_given) aligner.details_multiplexer(details.multiplexer());
  if(args.coords_given) aligner.coords_multiplexer(coords.multiplexer(), !args.no_header_flag);
  if(unitigs_lengths) aligner.unitigs_lengths(unitigs_lengths.get(), args.k_mer_arg);

  // Output matches
  aligner.exec_join(args.threads_arg);

  return 0;
}
