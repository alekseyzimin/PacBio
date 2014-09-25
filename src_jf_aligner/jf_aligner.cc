#include <memory>
#include <ios>
#include <thread>

#include <src_jf_aligner/jf_aligner_cmdline.hpp>
#include <src_jf_aligner/superread_parser.hpp>
#include <src_jf_aligner/pb_aligner.hpp>
#include <debug.hpp>

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

void print_coords_header(Multiplexer* m, bool compact) {
  Multiplexer::ostream o(m); // Write header
  o << "Rstart Rend Qstart Qend Nmers Rcons Qcons Rcover Qcover Rlen Qlen Stretch Offset Err";
  if(!compact)
    o << " Rname";
  o << " Qname\n";
  o.end_record();
}

void print_coords(Multiplexer::ostream& out, const std::string& pb_name, const size_t pb_size,
                  const bool compact, const align_pb::coords_info_type& coords,
                            const std::vector<int>& order) {
  const size_t nb_lines = coords.size();
  if(nb_lines == 0) return;

  if(compact)
    out << ">" << nb_lines << " " << pb_name << "\n";
  for(size_t i = 0; i < nb_lines; ++i) {
    const auto& it = coords[order[i]];
    out << it.rs << " " << it.re << " " << it.qs << " " << it.qe << " "
        << it.nb_mers << " "
        << it.pb_cons << " " << it.sr_cons << " "
        << it.pb_cover << " " << it.sr_cover << " "
        << pb_size << " " << it.ql
        << " " << it.stretch << " " << it.offset << " " << it.avg_err;
    if(!compact)
      out << " " << pb_name;
    out << " " << it.unitigs.name();
    auto mit = it.kmers_info.cbegin();
    auto bit = it.bases_info.cbegin();
    for( ; mit != it.kmers_info.cend(); ++mit, ++bit)
      out << " " << *mit << ":" << *bit;
    out << "\n";
  }
  out.end_record();
}

void print_details(Multiplexer::ostream& out, const std::string& pb_name, const align_pb::frags_pos_type& frags_pos) {
  for(auto it = frags_pos.cbegin(); it != frags_pos.cend(); ++it) {
    out << pb_name << " " << it->first;
    const align_pb::mer_lists& ml              = it->second;
    const bool                       fwd_align =
      std::distance(ml.fwd.lis.cbegin(), ml.fwd.lis.cend()) > std::distance(ml.bwd.lis.cbegin(), ml.bwd.lis.cend());
    auto                       lisit           = fwd_align ? ml.fwd.lis.cbegin() : ml.bwd.lis.cbegin();
    const auto                 lisend          = fwd_align ? ml.fwd.lis.cend() : ml.bwd.lis.cend();
    auto                       fwd_offit       = ml.fwd.offsets.cbegin();
    const auto                 fwd_offbegin    = fwd_offit;
    const auto                 fwd_offend      = ml.fwd.offsets.cend();
    auto                       bwd_offit       = ml.bwd.offsets.cbegin();
    const auto                 bwd_offbegin    = bwd_offit;
    const auto                 bwd_offend      = ml.bwd.offsets.cend();

    while(fwd_offit != fwd_offend || bwd_offit != bwd_offend) {
      std::pair<int, int> pos;
      bool                part_of_lis = false;
      if(fwd_offit != fwd_offend && (bwd_offit == bwd_offend || fwd_offit->first <= bwd_offit->first)) {
        pos = *fwd_offit;
        part_of_lis = fwd_align && (lisit < lisend) && (*lisit == std::distance(fwd_offbegin, fwd_offit));
        ++fwd_offit;
      } else if(bwd_offit != bwd_offend && (fwd_offit == fwd_offend || bwd_offit->first < fwd_offit->first)) {
        pos = *bwd_offit;
        part_of_lis = !fwd_align && (lisit < lisend) && (*lisit == std::distance(bwd_offbegin, bwd_offit));
        ++bwd_offit;
      }
      out << " " <<(part_of_lis ? "[" : "")
          << pos.first << ":" << pos.second
          << (part_of_lis ? "]" : "");
      if(part_of_lis)
        ++lisit;
    }
    out << "\n";
  }
  out.end_record();
}

void print_alignments(read_parser* reads, Multiplexer* details_m, Multiplexer* coords_m, const align_pb* align_data) {
  parse_sequence           parser;
  align_pb::thread         aligner(*align_data);

  mstream          details_io(details_m);
  mstream          coords_io(coords_m);
  std::string      name;
  std::vector<int> sort_array;

  while(true) {
    read_parser::job job(*reads);
    if(job.is_empty()) break;

    for(size_t i = 0; i < job->nb_filled; ++i) { // Process each read
      auto name_end = job->data[i].header.find_first_of(" \t\n\v\f\r");
      name = job->data[i].header.substr(0, name_end);
      const size_t pb_size = job->data[i].seq.size();
      parser.reset(job->data[i].seq);
      aligner.align_sequence_max(parser, pb_size);

      const auto& coords = aligner.coords();
      const int n = coords.size();
      if((int)sort_array.size() < n)
        sort_array.resize(n);
      for(int i = 0; i < n; ++i)
        sort_array[i] = i;
      std::sort(sort_array.begin(), sort_array.begin() + n, [&coords] (int i, int j) { return coords[i] < coords[j]; });
      print_coords(*coords_io, name, pb_size, args.compact_flag, coords, sort_array);
      if(details_io) print_details(*details_io, name, aligner.frags_pos());
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
    is >> unitig >> len;
    while(is.good()) {
      unitigs_lengths->push_back(len);
      is >> unitig >> len;
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
  align_pb align_data(hash, args.stretch_factor_arg, args.stretch_constant_arg, args.window_size_arg,
                      args.forward_flag, args.max_match_flag,
                      args.max_count_arg ? args.max_count_arg : std::numeric_limits<int>::max(),
                      args.mers_matching_arg / 100.0, args.bases_matching_arg / 100.0);
  if(unitigs_lengths) align_data.unitigs_lengths(unitigs_lengths.get(), args.k_mer_arg);

  // Output header if necessary
  if(!args.no_header_flag)
    print_coords_header(coords.multiplexer(), args.compact_flag);

  // Output alignements
  std::vector<std::thread> threads;
  for(unsigned int i = 0; i < args.threads_arg; ++i)
    threads.push_back(std::thread(print_alignments, &reads, details.multiplexer(), coords.multiplexer(),
                                  &align_data));
  for(auto& th : threads)
    th.join();

  return 0;
}
