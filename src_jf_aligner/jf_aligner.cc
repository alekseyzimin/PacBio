#include <memory>

#include <jflib/multiplexed_io.hpp>
#include <src_jf_aligner/jf_aligner_cmdline.hpp>
#include <src_jf_aligner/superread_parser.hpp>
#include <src_jf_aligner/pb_aligner.hpp>


int main(int argc, char *argv[])
{
  jf_aligner_cmdline args(argc, argv);
  mer_dna::k(args.mer_arg);

  // Open output file for early error reporting
  std::unique_ptr<std::ostream> out;
  if(args.output_given) {
    out.reset(new std::ofstream(args.output_arg));
    if(!out->good())
      jf_aligner_cmdline::error() << "Failed to open output file '" << args.output_arg << "'";
  } else {
    out.reset(new std::ostream(std::cout.rdbuf()));
  }

  // Read the super reads
  mer_pos_hash_type hash(args.size_arg);
  name_lists names(args.threads_arg);
  superread_parse(args.threads_arg, hash, names, args.superreads_arg.cbegin(), args.superreads_arg.cend());

  // Output matches
  stream_manager streams(args.pacbio_arg.cbegin(), args.pacbio_arg.cend());
  jflib::o_multiplexer multiplexer(out.get(), 4 * args.threads_arg, 4096);
  align_pb aligner(args.threads_arg, hash, streams, multiplexer);
  aligner.exec_join(args.threads_arg);

  return 0;
}
