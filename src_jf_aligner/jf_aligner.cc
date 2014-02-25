#include <memory>
#include <ios>

#include <jflib/multiplexed_io.hpp>
#include <src_jf_aligner/jf_aligner_cmdline.hpp>
#include <src_jf_aligner/superread_parser.hpp>
#include <src_jf_aligner/pb_aligner.hpp>

jf_aligner_cmdline args;


class output_file {
  std::unique_ptr<std::ostream>  file_;
  std::unique_ptr<o_multiplexer> multiplexer_;

public:
  output_file() = default;
  output_file(const char* path) {
    open(path);
  }
  void open(const char* path) {
    file_.reset(new std::ofstream(path));
    if(!file_->good())
      jf_aligner_cmdline::error() << "Failed to open file '" << path << "'";
    multiplexer_.reset(new o_multiplexer(file_.get(), 4 * args.threads_arg, 4096));
  }
  std::ostream& file() { return *file_; }
  o_multiplexer* multiplexer() { return multiplexer_.get(); }
};

int main(int argc, char *argv[])
{
  args.parse(argc, argv);
  mer_dna::k(args.mer_arg);
  std::ios::sync_with_stdio(false);

  // Open output file for early error reporting
  output_file details;
  output_file coords;
  if(args.details_given) details.open(args.details_arg);
  if(args.coords_given) coords.open(args.coords_arg);

  // Read the super reads
  mer_pos_hash_type hash(args.size_arg);
  name_lists names(args.threads_arg);
  stream_manager streams(args.pacbio_arg.cbegin(), args.pacbio_arg.cend());
  superread_parse(args.threads_arg, hash, names, args.superreads_arg.cbegin(), args.superreads_arg.cend());

  // Create aligner
  align_pb aligner(args.threads_arg, hash, streams, args.stretch_constant_arg, args.stretch_factor_arg);
  if(args.details_given) aligner.details_multiplexer(details.multiplexer());
  if(args.coords_given) {
    coords.file() << "S1 E1 S2 E2 N L1 L2 Q R" << std::endl;
    aligner.coords_multiplexer(coords.multiplexer());
  }

  // Output matches
  aligner.exec_join(args.threads_arg);

  return 0;
}
