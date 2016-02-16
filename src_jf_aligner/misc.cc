/******************************************
Copyright University of Maryland 2015
******************************************/
#include <src_jf_aligner/misc.hpp>
#include <limits>

inline static std::istream& skip_header(std::istream& is) {
  return is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

void read_unitigs_lengths(std::istream& is, std::vector<int>& lengths) {
  std::string  unitig;
  unsigned int len;
  is >> unitig >> len;
  while(is.good()) {
    lengths.push_back(len);
    is >> unitig >> len;
  }
}

void read_unitigs_sequences(std::istream& is, std::vector<int>& lengths) {
  std::string sequence;
  while(skip_header(is)) {
    std::getline(is, sequence);
    lengths.push_back(sequence.size());
  }
}


void read_unitigs_sequences(std::istream& is, std::vector<int>& lengths,
                            std::vector<std::string>& sequences) {
  while(skip_header(is)) {
    sequences.push_back("");
    std::getline(is, sequences.back());
    lengths.push_back(sequences.back().size());
  }
}
