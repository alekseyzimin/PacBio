#include <src_jf_aligner/coords_parsing.hpp>
#include <src_jf_aligner/frag_info.hpp>

//std::istream& operator>>(std::istream& is, align_pb::coords_info& c) {
void parse_coords(int thid, std::istream& is, align_pb::coords_info& c, frag_lists& frags) {
  std::string qname;
  is >> c.rs >> c.re >> c.qs >> c.qe >> c.nb_mers >> c.pb_cons >> c.sr_cons
     >> c.pb_cover >> c.sr_cover >> c.rl >> c.ql
     >> c.stretch >> c.offset >> c.avg_err
     >> qname;
  c.qfrag   = frags.push_back(thid, c.ql, qname.c_str());
  c.name_u  = &c.qfrag->fwd;
  c.kmers_info.clear();
  c.bases_info.clear();
  char sep;
  int mers, bases;
  while(is >> mers >> sep >> bases) {
    c.kmers_info.push_back(mers);
    c.bases_info.push_back(bases);
  }
}

inline static void skip_line(std::istream& is) {
  is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

void coords_parser::parser_loop() {
  std::string line;
  // Skip header
  while(file_.peek() != '>' && file_.peek() != EOF)
    skip_line(file_);
  while(file_.good()) {
    elt e(elt_init());
    size_type& i = e->nb_filled;
    for(i = 0; i < group_size(); ++i) {
      auto& elt = e->elements[i];
      if(!std::getline(file_, line)) break;
      if(__builtin_expect(line[0] != '>', 0)) {
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
        if(__builtin_expect(!std::getline(file_, elt.lines[j]), 0)) {
          std::cerr << "Invalid input file. File truncated" << std::endl;
          file_.close();
          break;
        }
      }
    }
  }
}
