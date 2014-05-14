#include <random>
#include <iostream>
#include <fstream>
#include <gtest/gtest.h>
#include <src_jf_aligner/least_square_2d.hpp>



namespace {
struct least_square_2d_trivial {
  long n;
  double sum_x, sum_y, sum_xy, sum_xx;
  least_square_2d_trivial() :
    n(0), sum_x(0), sum_y(0), sum_xy(0), sum_xx(0)
  { }

  void add(double x, double y) {
    ++n;
    sum_x  += x;
    sum_y  += y;
    sum_xy += x * y;
    sum_xx += x * x;
  }

  double var_x() const { return n * sum_xx - sum_x * sum_x; }
  double covar_xy() const { return n * sum_xy - sum_x * sum_y; }
  double nb() const { return sum_xx * sum_y - sum_x * sum_xy; }
  double a() const { return covar_xy() / var_x(); }
  double b() const { return nb() / var_x(); }
};

TEST(LeastSquare, Trivial) {
  static const int nb_points = 1000;

  std::default_random_engine gen;
  std::normal_distribution<double> normal;
  std::uniform_real_distribution<double> uniform(-10.0, 10.0);
  const double a = uniform(gen);
  const double b = uniform(gen);

  std::vector<std::pair<double, double> > points;
  least_square_2d         ls;
  least_square_2d_trivial lst;

  for(int i = 0; i < nb_points; ++i) {
    const double x = uniform(gen);
    const double y = a * x + b + normal(gen);
    points.push_back(std::make_pair(x, y));
    ls.add(x, y);
    lst.add(x, y);
  }

  EXPECT_NEAR(a, lst.a(), 1e-1);
  EXPECT_NEAR(b, lst.b(), 1e-1);
  EXPECT_NEAR(lst.a(), ls.a(), 1e-4);
  EXPECT_NEAR(lst.b(), ls.b(), 1e-4);

  double r = 0;
  const double a_ = ls.a();
  const double b_ = ls.b();
  for(auto v : points)
    r += a_ * v.first + b_ - v.second;
  EXPECT_NEAR(0, r, 1e-4);
} // LeastSquare.Trivial

} // namespace
