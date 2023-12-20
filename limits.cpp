#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

template <typename FunType>
std::pair<double, bool> Limit(FunType f, double x0, double h = 0, double eps = 1e-8, double nmax = 20) {

  if (nmax<3||nmax>30||eps<=0) throw std::domain_error("Invalid parameters");
  if (std::abs(h) < eps) {
    if (x0<=1) h=0.001;
    else h=0.001*std::abs(x0);
  }
  std::vector<double>y(nmax,0);
  int t;
  double yl=std::numeric_limits<double>::infinity();

  for (int i=0;i<nmax;i++) {

    y[i]=f(x0+h);
    t=2;
    for (int k=i-1;k>=0;k--){
      y[k]=(t*y[k+1]-y[k])/(t-1);
      t*=2;
    }

    if (std::fabs(y[0]-yl)<eps) return std::pair<double, bool>(y[0],true);
    yl=y[0];
    h=h/2;
  }
  
  return std::pair<double, bool>(y[0], false);
}

int main() {
  // AT1 - Limit test: Limes (tan(x)-sin(x))/(x*x*x) x->0:
  auto limes = Limit(
      [](double x) { return (std::tan(x) - std::sin(x)) / (x * x * x); }, 0);
  std::cout << limes.first << " " << limes.second;
  return 0;
}