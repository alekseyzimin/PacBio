#include <cerrno>
#include <climits>
#include <ios>
#include <thread>

#include <src_jf_aligner/longest_path_overlap_graph2_cmdline.hpp>
#include <src_jf_aligner/output_file.hpp>
#include <src_jf_aligner/pb_aligner.hpp>
#include <src_jf_aligner/overlap_graph.hpp>
#include <jflib/multiplexed_parser.hpp>

typedef longest_path_overlap_graph2_cmdline cmdline_args;
cmdline_args args;

struct coords_lines {
  std::string              header;
  std::vector<std::string> lines;
};

class coords_parser : public multiplexed_parser<coords_lines> {
  std::ifstream file_;
public:
  coords_parser(int nb_threads, const char* file) : multiplexed_parser(nb_threads), file_(file) { }

  void parser_loop() {
    std::string line;
    while(file_.good()) {
      elt e(elt_init());
      size_type& i = e->nb_filled;
      for(i = 0; i < group_size(); ++i) {
        auto& elt = e->elements[i];
        if(!std::getline(file_, line)) break;
        if(line[0] != '>') {
          std::cerr << "Invalid input file. Line expected to match /^>/ but got: " << line << std::endl;
          file_.close();
          break;
        }
        char* endptr = 0;
        errno = 0;
        const long nb_lines = std::strtol(line.c_str() + 1, &endptr, 10);
        if(nb_lines == 0 || errno == ERANGE) {
          std::cerr << "Invalid input file. Expected number of lines but got: " << (line.c_str() + 1) << std::endl;
          file_.close();
          break;
        }
        elt.header = endptr + 1;
        elt.lines.resize(nb_lines);
        for(long j = 0; j < nb_lines; ++j) {
          if(!std::getline(file_, elt.lines[j])) {
            std::cerr << "Invalid input file. File truncated" << std::endl;
            file_.close();
            break;
          }
        }
      }
    }
  }
};

void fill_coords(std::vector<std::string>& lines, coords_info_type& coords) {
  coords.resize(lines.size());
  for(size_t i = 0; i < lines.size(); ++i) {
    std::istringstream is(lines[i]);
    auto&              c = coords[i];
    is >> c.rs >> c.re >> c.qs >> c.qe >> c.nb_mers >> c.pb_cons >> c.sr_cons
       >> c.pb_cover >> c.sr_cover >> c.rl >> c.ql
       >> c.stretch >> c.offset >> c.avg_err
       >> c.qname;
    c.unitigs = qname;
    c.kmers_info.clear();
    c.bases_info.clear();
    char sep;
    int mers, bases;
    while(is >> mers >> sep >> bases) {
      c.kmers_info.push_back(mers);
      c.base_info.push_back(mers);
    }
  }
}

void create_mega_reads(coords_parser* parser, Multiplexer* output_m, overlap_graph* graph_walker) {
  coords_info_type coords;

  for(coords_parser::stream coords_stream(*parser); coords_stream; ++coords_stream) {
    
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
