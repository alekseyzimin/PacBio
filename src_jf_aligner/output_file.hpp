#ifndef _OUTPUT_FILE_H_
#define _OUTPUT_FILE_H_

#include <fstream>
#include <multiplexer.hpp>
#include <jellyfish/err.hpp>
namespace err = jellyfish::err;

// Manage an output file.
class output_file {
  std::unique_ptr<std::ostream> file_;
  std::unique_ptr<Multiplexer>  multiplexer_;

public:
  output_file() = default;
  output_file(const char* path, int threads) { open(path, threads); }
  output_file(std::ostream& out, int threads) { set(out, threads); }
  void open(const char* path, int threads) {
    file_.reset(new std::ofstream(path));
    if(!file_->good())
      throw std::runtime_error(err::msg() << "Failed to open file '" << path << "'");
    multiplexer_.reset(new Multiplexer(*file_.get(), 1024 * 1024, 10 * 1024 * 1024));
  }
  void set(std::ostream& out, int threads) {
    file_.reset(&out);
    multiplexer_.reset(new Multiplexer(*file_.get(), 1024 * 1024, 10 * 1024 * 1024));
  }
  std::ostream& file() { return *file_; }
  Multiplexer* multiplexer() { return multiplexer_.get(); }
};


#endif /* _OUTPUT_FILE_H_ */
