#ifndef __LEAST_SQUARE_2D_H__
#define __LEAST_SQUARE_2D_H__

#include <cstdlib>


/* Compute the parameters a and b for the least square optimization:
 * ||a*X + b -Y||^2. They are computed as:
 *
 * a = cov(X,Y) / var(X)
 * b = (cov(XY,X) - cov(XX,Y)) / var(X)
 *
 * where cov(,) is the covariance and var is the variance. It is
 * computed with a stable online algorithm. If the points are (x1,y1),
 * (x2,y2), ..., then:
 *
 * expectancy of {x1,x2,...,xk}
 *  E(X,0) = 0
 *  E(X,k) = E(X,k-1) + (xk - E(X,k-1)) / k
 *
 * covariance of {(x1,y1),(x2,y2),...,(xk,yk)}
 *  COV(X,Y,0) = 0
 *  COV(X,Y,k) = COV(X,Y,k-1) + (xk - E(X,k)) * (yk - E(Y,k-1))
 *             = COV(X,Y,k-1) + (xk - E(X,k-1)) * (yk - E(Y,k))
 *
 * variance of {x1,x2,....,xk}
 *  VAR(X,0) = 0
 *  VAR(X,k) = VAR(X,k-1) + (xk - E(X,k-1)) * (yk - E(X,k))
 *
 * the numerator of b
 *  NB(X,Y,0) = 0
 *  NB(X,Y,k) = NB(X,Y,k-1) + (xk*yk - E(XY,k-1)) * (xk * xk - E(XX,k-1)) * (yk - E(Y,k))
 */
struct least_square_2d {
  double EX, EY, EXX, EXY;
  double VX, CXY, NB;
  long n;
  least_square_2d() :
    EX(0), EY(0), EXX(0), EXY(0),
    VX(0), CXY(0), NB(0),
    n(0)
  { }

  void add(const double x, const double y) {
    ++n;

    const double deltaX = x - EX;
    EX += deltaX / n;
    const double ndeltaX = x - EX;
    VX += deltaX * ndeltaX;

    const double deltaY = y - EY;
    EY += deltaY / n;
    const double ndeltaY = y - EY;

    const double deltaXX = x * x - EXX;
    EXX += deltaXX / n;

    const double deltaXY = x * y - EXY;
    EXY += deltaXY / n;

    CXY += deltaX * ndeltaY;
    NB += deltaXY * ndeltaX - deltaXX * ndeltaY;
  }
  void add(std::pair<double, double>& p) {
    add(p.first, p.second);
  }

  // template<typename Iterator>
  // void add(Iterator start, const Iterator end) {
  //   for( ; start != end; ++start)
  //     add(start->first, start->second);
  // }

  double a() const { return CXY / VX; }
  double b() const { return NB / VX; }

  // template<typename Iterator>
  // double error(Iterator start, const Iterator end) {
  //   double       e  = 0;
  //   const double a_ = a();
  //   const double b_ = b();
  //   long         n_ = 0;
  //   for( ; start != end; ++start, ++n_)
  //     e += abs(a_ * start->first + b_ - start->second);
  //   return e / n_;
  // }
};


#endif /* __LEAST_SQUARE_2D_H__ */
