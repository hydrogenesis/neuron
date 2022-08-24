#pragma once
#include <vector>
#include <cmath>

using std::vector;

class HH {
private:
public:
  long long i;
  double I;
  double dt;
  long long loop;
  double v;
  double mi;
  double hi;
  double ni;
  double gNa;
  double eNa;
  double gK;
  double eK;
  double gL;
  double eL;
  double g;
  vector<double> t;
  vector<double> V;
  vector<double> m;
  vector<double> h;
  vector<double> n;
  vector<HH*> prev;
  vector<HH*> next;

  HH(double I, double tspan, double dt, double v, double mi, double hi, double ni,
      double gNa, double eNa, double gK, double eK, double gL, double eL, double g);
  ~HH();

  double DynamicI(double t, double I) const;
  void Append(HH* next_hh);
  inline void Advance() { ++i; }
  void Advance_Euler();
  inline double AlphaM(double V) const {
    return (2.5-0.1*(V+65)) / (std::exp(2.5-0.1*(V+65)) - 1);
  }
  inline double BetaM(double V) const {
    return 4.0*std::exp(-(V+65)/18);
  }
  inline double AlphaH(double V) const {
    return 0.07*std::exp(-(V+65)/20);
  }
  inline double BetaH(double V) const {
    return 1.0/(std::exp(3.0-0.1*(V+65))+1);
  }
  inline double AlphaN(double V) const {
    return (0.1-0.01*(V+65)) / (std::exp(1-0.1*(V+65)) -1);
  }
  inline double BetaN(double V) const {
    return 0.125*std::exp(-(V+65)/80);
  }
};  // class HH
