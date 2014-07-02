#include <memory>
#include <ios>
#include <thread>

#include <src_jf_aligner/create_mega_reads_cmdline.hpp>
#include <src_jf_aligner/superread_parser.hpp>
#include <src_jf_aligner/pb_aligner.hpp>
#include <src_jf_aligner/union_find.hpp>

create_mega_reads_cmdline args;

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

struct node_info {
  bool            start_node;
  bool            end_node;
  double          imp_s, imp_e;
  union_find::set component;

  node_info() = default;

  void reset(const align_pb::coords_info& info) {
    start_node = true;
    end_node   = true;
    imp_s = info.imp_s();
    imp_e = info.imp_e();
    component.reset();
  }
};

struct overlap_graph {
  double overlap_play;
  unsigned int k_len;
  const std::vector<int>& unitigs_lengths;
  double nb_errors;

  overlap_graph(double play, unsigned int len, const std::vector<int>& lengths, double errors) :
    overlap_play(play), k_len(len), unitigs_lengths(lengths), nb_errors(errors)
  { }

  void traverse(const std::vector<int>& sort_array, const align_pb::coords_info_type& coords,
                std::vector<node_info>& nodes) {
    for(size_t it_i = 0; it_i != sort_array.size(); ++it_i) {
      auto&       node_i   = nodes[it_i];
      const auto& coords_i = coords[it_i];
      if(node_i.imp_e >= coords_i.rl)
        continue; // Hanging off 3' end of PacBio read: no overlap
      for(size_t it_j = it_i + 1; it_j != sort_array.size(); ++it_j) {
        auto&       node_j   = nodes[it_j];
        const auto& coords_j = coords[it_j];
        if(node_j.imp_s <= 1)
          continue; // Hanging off 5' end of PacBio read: no overlap
        const double position_len = node_i.imp_e - node_j.imp_s + 1;
        if(position_len * overlap_play < k_len)
          break; // maximum implied overlap is less than a k-mer length
        const int nb_u_overlap = coords_i.unitigs.overlap(coords_j.unitigs);
        if(!nb_u_overlap)
          continue; // No overlap according to unitig names
        int u_overlap_len = 0;
        for(int i = 0; i < nb_u_overlap; ++i)
          u_overlap_len += unitigs_lengths[coords_j.unitigs.unitig_id(i)];
        u_overlap_len -= (nb_u_overlap - 1) * k_len;
        double error = nb_errors * (coords_i.avg_err + coords_j.avg_err);
        if(u_overlap_len > overlap_play * position_len + error || position_len > overlap_play * u_overlap_len + error)
          continue; // Overlap lengths (position, unitigs) do not agree

        // We have an overlap between nodes i and j
        node_i.end_node    = false;
        node_j.start_node  = false;
        node_i.component  |= node_j.component;
        
      }
    }
  }
};


void create_mega_reads(read_parser* reads, Multiplexer* output_m, const align_pb* align_data,
                       overlap_graph* graph_walker) {
  mer_dna                tmp_m;
  parse_sequence         parser;
  align_pb::thread       aligner(*align_data);
  std::vector<int>       sort_array;
  std::vector<node_info> nodes;

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

      const auto& coords = aligner.coords();
      const int n = coords.size();
      if((int)sort_array.size() < n)
        sort_array.resize(n);
      if((int)nodes.size() < n)
        nodes.resize(n);
      for(int i = 0; i < n; ++i) {
        sort_array[i] = i;
        nodes[i].reset(coords[i]);
      }
      std::sort(sort_array.begin(), sort_array.end() + n,
                [&nodes] (int i, int j) { return nodes[i].imp_s < nodes[j].imp_s || (nodes[i].imp_s == nodes[j].imp_s &&
                                                                                     nodes[i].imp_e < nodes[j].imp_e); });
      graph_walker->traverse(sort_array, coords, nodes);
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
                      true /* forward */, args.max_match_flag, args.max_count_arg,
                      args.mers_matching_arg / 100.0, args.bases_matching_arg / 100.0);
  align_data.unitigs_lengths(&unitigs_lengths, args.k_mer_arg);

  // Output candidate mega_reads
  overlap_graph graph_walker(args.overlap_play_arg, args.k_mer_arg, unitigs_lengths, args.errors_arg);
  std::vector<std::thread> threads;
  for(unsigned int i = 0; i < args.threads_arg; ++i)
    threads.push_back(std::thread(create_mega_reads, &reads, output.multiplexer(), &align_data,
                                  &graph_walker));
  for(unsigned int i = 0; i < args.threads_arg; ++i)
    threads[i].join();

  return 0;
}
