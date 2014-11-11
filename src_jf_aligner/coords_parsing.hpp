#ifndef __COORDS_PARSING_H__
#define __COORDS_PARSING_H__

#include <istream>
#include <vector>
#include <string>
#include <jflib/multiplexed_parser.hpp>
#include <src_jf_aligner/pb_aligner.hpp>

// parse one coords_info record in compact format
//std::istream& operator>>(std::istream& is, align_pb::coords_info& c);
void parse_coords(int thid, std::istream& is, align_pb::coords_info& c, frag_lists& frags);

struct coords_lines {
  std::string              header;
  std::vector<std::string> lines;
};

// Multiplexed parser for compact coords file
class coords_parser : public multiplexed_parser<coords_lines> {
  std::ifstream file_;
public:
  coords_parser(int nb_threads, const char* file) : multiplexed_parser(nb_threads), file_(file) { }

  void parser_loop();
};



#endif /* __COORDS_PARSING_H__ */
