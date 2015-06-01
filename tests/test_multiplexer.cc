#include <fstream>
#include <thread>
#include <vector>

#include <gtest/gtest.h>
#include <multiplexer.hpp>
#include <tests/misc.hpp>

namespace {
using misc::remove_file;

TEST(Multiplexer, OneThread) {
  remove_file tmp;

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
  EXPECT_TRUE((bool)std::getline(is, line));
  EXPECT_EQ("Hello", line);
  EXPECT_TRUE((bool)std::getline(is, line));
  EXPECT_EQ("5", line);
  EXPECT_TRUE((bool)std::getline(is, line));
  EXPECT_EQ("6", line);
  EXPECT_TRUE((bool)std::getline(is, line));
  EXPECT_EQ("toto", line);
  EXPECT_TRUE((bool)std::getline(is, line));
  EXPECT_EQ("Test 3.14", line);
  EXPECT_FALSE((bool)std::getline(is, line));
} // Multiplexer.OneThread

void output_lines(int id, Multiplexer* m, const int nb_lines) {
  Multiplexer::ostream os(m);

  for(int i = 0; i < nb_lines; ++i) {
    os << id << " line " << i << "\n";
    if(i % 10 == 9)
      os.end_record();
  }
}

TEST(Multiplexer, ManyThreads) {
  static const int nb_threads = 10;
  static const int lines_per_thread = 1000;
  remove_file tmp;

  {
    std::ofstream os(tmp.path);
    Multiplexer m(os);
    std::vector<std::thread> threads;

    for(int i = 0; i < nb_threads; ++i)
      threads.push_back(std::thread(output_lines, i, &m, lines_per_thread));
    for(auto& th: threads) th.join();
  }

  std::set<std::pair<int, int> > res;
  std::ifstream is(tmp.path);
  std::string line;
  while(std::getline(is, line)) {
    int thid = std::stoi(line);
    EXPECT_TRUE(thid >= 0 && thid < nb_threads);
    size_t lstart = line.find_first_not_of("0123456789");
    EXPECT_EQ(" line ", line.substr(lstart, 6));
    int line_i = std::stoi(line.substr(lstart + 6));
    EXPECT_TRUE(line_i >= 0 && line_i < lines_per_thread);
    //    std::cerr << line << " -> " << thid << " " << line_i << "\n";
    auto insert_result = res.insert(std::make_pair(thid, line_i));
    EXPECT_TRUE(insert_result.second);
  }
  ASSERT_EQ((size_t)(nb_threads * lines_per_thread), res.size());
  auto it = res.cbegin();
  for(int i = 0; i < nb_threads; ++i) {
    for(int j = 0; j < lines_per_thread; ++j, ++it) {
      ASSERT_NE(res.cend(), it);
      EXPECT_EQ(i, it->first);
      EXPECT_EQ(j, it->second);
    }
  }
  ASSERT_EQ(res.cend(), it);
}
} // namespace
