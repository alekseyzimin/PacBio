/******************************************
Copyright University of Maryland 2015
******************************************/
#include <vector>
#include <fstream>

#include <src_jf_aligner/coords_parsing.hpp>
#include <src_jf_aligner/merge_coords_cmdline.hpp>
#include "zstr.hpp"


int main(int argc, char *argv[]) {
  merge_coords_cmdline        args(argc, argv);

  // Open output
  std::ostream os(std::cout.rdbuf());
  std::ofstream output_file;
  if(args.output_given) {
    output_file.open(args.output_arg);
    if(!output_file.good())
      merge_coords_cmdline::error() << "Error opening output file '" << args.output_arg << '\'';
    os.rdbuf(output_file.rdbuf());
  }

  // Special cases: no merging needed
  if(args.coords_arg.size() == 0)
    return 0;
  if(args.coords_arg.size() == 1) {
    std::ifstream is(args.coords_arg[0]);
    if(!is.good())
      merge_coords_cmdline::error() << "Error opening coords file '" << args.coords_arg[0] << '\'';
    os << is.rdbuf();
    return 0;
  }

  // General case: merge inputs
  frag_lists                  frags(1);
  std::vector<zstr::ifstream*> inputs;

  for(auto path : args.coords_arg) {
    auto is = new zstr::ifstream(path);
    if(!is->good())
      merge_coords_cmdline::error() << "Error opening coords file '" << path << '\'';
    inputs.push_back(is);
  }

  std::string              pb_name;
  std::string              line;
  std::vector<std::string> lines;
  while(inputs.front()->peek() == '>') {
    pb_name.clear();
    lines.clear();
    for(auto is : inputs) {
      if(is->peek() != '>')
        merge_coords_cmdline::error() << "Reached end of file prematurely";
      std::getline(*is, line);
      char* next;
      const int nb_aligns = std::strtol(line.c_str() + 1, &next, 10);
      if(!*next)
        merge_coords_cmdline::error() << "Invalid format: query sequence name missing";
      if(!pb_name.empty()) {
        if(pb_name != next + 1)
          merge_coords_cmdline::error() << "Invalid order of query sequence: expected '"
                                        << pb_name << "' and got '" << next;
      } else {
        pb_name = next + 1;
      }
      for(int i = 0; i < nb_aligns; ++i) {
        std::getline(*is, line);
        lines.push_back(line);
      }
    }
    os << '>' << lines.size() << ' ' << pb_name << '\n';
    for(const auto& l : lines)
      os << l << '\n';
  }

  if(!std::all_of(inputs.cbegin(), inputs.cend(), [](std::istream* is) { return is->peek() == EOF; }))
        merge_coords_cmdline::error() << "Reached end of file prematurely";

  output_file.close();

  return 0;
}
