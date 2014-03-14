#include <fstream>
#include <thread>
#include <vector>

#include <gtest/gtest.h>
#include <multiplexer.hpp>

namespace {
struct remove_file {
  const char* path;
  bool do_unlink;
  remove_file(const char* p, bool unlink = true) : path(p), do_unlink(unlink) { }
  ~remove_file() { if(do_unlink) unlink(path); }
};

TEST(Multiplexer, OneThread) {
  remove_file tmp("/tmp/output");
  tmp.do_unlink = false;

  {
    std::ofstream os(tmp.path);
    Multiplexer m(os);
    Multiplexer::ostream mos(m);
    mos << "Hello\n" << 5 << "\n";
    mos << 6 << "\n" << "toto" << "\n";
    //    mos.end_record();
    mos.do_flush();
    mos << "Test " << 3.14 << "\n";
  }

  std::ifstream is(tmp.path);
  std::string line;
  EXPECT_TRUE(std::getline(is, line));
  EXPECT_EQ("Hello", line);
  EXPECT_TRUE(std::getline(is, line));
  EXPECT_EQ("5", line);
  EXPECT_TRUE(std::getline(is, line));
  EXPECT_EQ("6", line);
  EXPECT_TRUE(std::getline(is, line));
  EXPECT_EQ("toto", line);
  EXPECT_TRUE(std::getline(is, line));
  EXPECT_EQ("Test 3.14", line);
  EXPECT_FALSE(std::getline(is, line));
} // Multiplexer.OneThread

void output_lines(int id, Multiplexer* m) {
  Multiplexer::ostream os(m);
  std::cout << id << "\n";
  for(int i = 0; i < 1000; ++i) {
    if(i % 10 == 9)
      os.end_record();
  }
}

TEST(Multiplexer, ManyThreads) {
  static const int nb_threads = 10;
  remove_file tmp("/tmp/manythreads_output");
  tmp.do_unlink = false;

  {
    std::ofstream os(tmp.path);
    Multiplexer m(os);
    std::vector<std::thread> threads;

    for(int i = 0; i < nb_threads; ++i)
      threads.push_back(std::thread(output_lines, i, &m));
    for(auto& th: threads) th.join();
  }

  //  std::set<std::pair<int, int> > res;
  std::ifstream is(tmp.path);
  std::string line;
  while(std::getline(is, line)) {
    int thid = std::stoi(line);
    EXPECT_TRUE(thid >= 0 && thid <= nb_threads);
    size_t lstart = line.find_first_not_of("0123456789");
    EXPECT_EQ(" line ", line.substr(lstart, 6));
  }
}
} // namespace
