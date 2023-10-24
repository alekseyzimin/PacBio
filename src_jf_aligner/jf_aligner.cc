/******************************************
Copyright University of Maryland 2015
******************************************/
#include <memory>
#include <ios>
#include <thread>

#include <src_jf_aligner/superread_parser.hpp>
#include <src_jf_aligner/coarse_aligner.hpp>
#include <src_jf_aligner/fine_aligner.hpp>
#include <src_jf_aligner/misc.hpp>
#include <debug.hpp>
#include <src_jf_aligner/jf_aligner_cmdline.hpp>

using align_pb::coarse_aligner;
using align_pb::fine_aligner;
using align_pb::coords_info_type;

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
                  const size_t format, const align_pb::coords_info_type& coords,
                  const std::vector<int>& order, bool zero_skip = true, const size_t num_matches = 0) {
  const size_t nb_lines = coords.size();
  size_t min_output_coverage = 0;
  size_t matchq = 60;
  if(nb_lines == 2) 
    matchq = 30;
  else
    if(nb_lines > 2) 
      matchq = 0;
  if(nb_lines == 0 && zero_skip) return;
  //format is 0 is compact, 1 is full and 2 is paf
  if(format==2){//PAF
    if(num_matches > 0){//determine minimum output coverage;
      std::vector<size_t> coverages;
      for(size_t i = 0; i < nb_lines; ++i) {
        const auto& it = coords[order[i]];
        coverages.push_back(it.pb_cover);
      }
      if(coverages.size() > num_matches){
        std::sort(coverages.begin(),coverages.end());
        min_output_coverage = coverages[coverages.size()-num_matches];
      }
    }
    size_t matches_output = 0;
    for(size_t i = 0; i < nb_lines; ++i) {
      const auto& it = coords[order[i]];
      char* ori="+";
      size_t astart=it.qs;
      size_t aend=it.qe;
      if(astart>aend){
        ori="-";
        astart=it.qe;
        aend=it.qs;
      }
      if((it.pb_cover >= min_output_coverage && matches_output <= num_matches) || num_matches == 0){
        out << pb_name << "\t" << pb_size << "\t" << it.rs << "\t" << it.re << "\t" <<  ori << "\t" << it.name_u->name << "\t" << it.ql << "\t" << astart << "\t" << aend << "\t" << it.pb_cover << "\t" << it.nb_mers << "\t" << matchq << "\n";
        matches_output++;
      }
      //out << "ALL matches:" << matches_output << " min_cov:" << min_output_coverage << " " << pb_name << "\t" << pb_size << "\t" << it.rs << "\t" << it.re << "\t" <<  ori << "\t" << it.name_u->name << "\t" << it.ql << "\t" << astart << "\t" << aend << "\t" << it.pb_cover << "\t" << it.nb_mers << "\t255\n";
    }
  }else{
      if(format==1)
      out << ">" << nb_lines << " " << pb_name << "\n";
      for(size_t i = 0; i < nb_lines; ++i) {
        const auto& it = coords[order[i]];
        if(format==0)
          out << pb_name << " ";
        out << it.rs << " " << it.re << " " << it.qs << " " << it.qe << " "
            << it.nb_mers << " "
            << it.pb_cons << " " << it.sr_cons << " "
            << it.pb_cover << " " << it.sr_cover << " "
            << pb_size << " " << it.ql
            << " " << it.stretch << " " << it.offset << " " << it.avg_err;
//    if(!compact)
//      out << " " << pb_name;
        assert(it.name_u == &it.qfrag->fwd || it.name_u == &it.qfrag->bwd);
        out << " " << it.name_u->name;
        auto mit = it.kmers_info.cbegin();
        auto bit = it.bases_info.cbegin();
        for( ; mit != it.kmers_info.cend(); ++mit, ++bit)
          out << " " << *mit << ":" << *bit;
        out << "\n";
      }
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

void print_alignments(read_parser* reads, Multiplexer* details_m, Multiplexer* coords_m,
                      const coarse_aligner* align_data, const fine_aligner* short_align_data,
                      bool skip_zero) {
  parse_sequence                        parser;
  coarse_aligner::thread                aligner(*align_data);
  std::unique_ptr<fine_aligner::thread> short_aligner;
  std::unique_ptr<short_parse_sequence> short_parser;

  mstream          details_io(details_m);
  mstream          coords_io(coords_m);
  std::string      name;
  std::vector<int> sort_array;

  if(short_align_data) {
    short_aligner.reset(new fine_aligner::thread(*short_align_data));
    short_parser.reset(new short_parse_sequence);
  }

  while(true) {
    read_parser::job job(*reads);
    if(job.is_empty()) break;

    for(size_t i = 0; i < job->nb_filled; ++i) { // Process each read
      auto name_end = job->data[i].header.find_first_of(" \t\n\v\f\r");
      name = job->data[i].header.substr(0, name_end);
      const size_t pb_size = job->data[i].seq.size();
      parser.reset(job->data[i].seq);
      aligner.align_sequence_max(parser, pb_size);

      const coords_info_type* coords = &aligner.coords();

      if(short_aligner) { // Lets recompute the alignment with smaller mers
        short_parser->reset(job->data[i].seq);
        short_aligner->align_sequence(*short_parser, pb_size, *coords);
        coords = &short_aligner->coords();
      }

      //      std::cerr << name << ' ' << coords->size() << '\n';
      const int n = coords->size();
      if((int)sort_array.size() < n)
        sort_array.resize(n);
      for(int i = 0; i < n; ++i)
        sort_array[i] = i;
      std::sort(sort_array.begin(), sort_array.begin() + n, [coords] (int i, int j) { return (*coords)[i] < (*coords)[j]; });
      print_coords(*coords_io, name, pb_size, args.format_arg, *coords, sort_array, skip_zero, args.num_matches_arg);
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
  if(args.unitigs_lengths_given || args.unitigs_sequences_given) {
    if(!args.k_mer_given)
      jf_aligner_cmdline::error() << "The mer length used for generating the k-unitigs (-k, --k-mer) is required if the unitig lengths (-l, --unitig-lengths or -u, --unitigs-sequences) is passed.";
    unitigs_lengths.reset(new std::vector<int>);
    if(args.unitigs_lengths_given) {
      std::ifstream is(args.unitigs_lengths_arg);
      if(!is.good())
        jf_aligner_cmdline::error() << "Failed to open unitig lengths map file '" << args.unitigs_lengths_arg << "'";
      read_unitigs_lengths(is, *unitigs_lengths);
    } else {
      std::ifstream is(args.unitigs_sequences_arg);
      if(!is.good())
        jf_aligner_cmdline::error() << "Failed to open unitig sequence file '" << args.unitigs_sequences_arg << "'";
      read_unitigs_sequences(is, *unitigs_lengths);
    }
  }

  // Read the super reads
  if(args.fine_mer_given)
    short_mer_type::k(args.fine_mer_arg);
  auto psa = superread_parse(args.reference_arg.cbegin(), args.reference_arg.cend(),
                             std::min(short_mer_type::k(), args.psa_min_arg), mer_dna::k());

  // Prepare I/O
  stream_manager streams(args.query_arg.cbegin(), args.query_arg.cend());
  read_parser    reads(4 * args.threads_arg, 100, 1, streams);

  // Create aligners
  coarse_aligner align_data(psa, mer_dna::k(),
                            args.stretch_factor_arg, args.stretch_constant_arg, args.stretch_cap_arg, args.window_size_arg,
                            args.forward_flag, args.max_match_flag,
                            args.max_count_arg ? args.max_count_arg : std::numeric_limits<int>::max(),
                            args.mers_matching_arg / 100.0, args.bases_matching_arg / 100.0);
  if(unitigs_lengths) align_data.unitigs_lengths(unitigs_lengths.get(), args.k_mer_arg);
  std::unique_ptr<fine_aligner> short_align_data;
  if(args.fine_mer_given)
    short_align_data.reset(new fine_aligner(psa, args.fine_mer_arg, unitigs_lengths.get(), args.k_mer_arg));

  // Output header if necessary
  if(!args.no_header_flag && args.format_arg <= 1)
    print_coords_header(coords.multiplexer(), args.format_arg);

  // Output alignements
  std::vector<std::thread> threads;
  for(unsigned int i = 0; i < args.threads_arg; ++i)
    threads.push_back(std::thread(print_alignments, &reads, details.multiplexer(), coords.multiplexer(),
                                  &align_data, short_align_data.get(), !args.zero_match_flag));
  for(auto& th : threads)
    th.join();

  return 0;
}
