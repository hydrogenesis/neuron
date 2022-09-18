#pragma once
#include <vector>
#include <cmath>
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

using std::vector;

class HH {
private:
public:
  long long i;
  long long index;
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
  double C;
  double t;
  double V;
  double m;
  double h;
  double n;
  double V1;
  double m1;
  double h1;
  double n1;
  long long i_backup;
  double dt_backup;
  double t_backup;
  double V_backup;
  double m_backup;
  double h_backup;
  double n_backup;
  double x, y, z;
  int type;
  double radius;
  double area;
  double length;
  vector<HH*> prev;
  vector<HH*> next;

  CUDA_HOSTDEV HH();
  HH(double I, double tspan, double dt, double v, double mi, double hi, double ni,
      double gNa, double eNa, double gK, double eK, double gL, double eL, double g);
  CUDA_HOSTDEV ~HH();

  void Backup();
  void Restore();
  double DynamicI(double t, double I) const;
  void Append(HH* next_hh);
  void SetParams(double x, double y, double z, int type, double length, double radius);
  void Advance();
  inline void Record(std::vector<double>* t_rec, std::vector<double>* V_rec) { t_rec->push_back(t); V_rec->push_back(V); }
  void Advance_Euler();
  void Advance_Crank_Nicolson(int iter_num);
  inline double Vtrap(double x, double y) const {
    if (std::abs(x/y) < 1e-6)
      return y*(1 - x/y/2.0);
    else
      return x/(std::exp(x/y) - 1);
  }
  inline double AlphaM(double V) const {
    return .1 * Vtrap(-(V+40),10);
    // double x = 2.5-0.1*(V+65);
    // if (std::abs(x) > 1e-6) {
    //   return (2.5-0.1*(V+65)) / (std::exp(x) - 1);
    // } else {
    //   return (2.5-0.1*(V+65)) / (0.5*x - 1);
    // }
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
    return .01*Vtrap(-(V+55),10);
    // double x = 1-0.1*(V+65);
    // if (std::abs(x) > 1e-6) {
    //   return (0.1-0.01*(V+65)) / (std::exp(x)-1);
    // } else {
    //   return (0.1-0.01*(V+65)) / (0.5*x-1);
    // }
  }
  inline double BetaN(double V) const {
    return 0.125*std::exp(-(V+65)/80);
  }
  CUDA_HOSTDEV void PrintDebugInfo() const;
};  // class HH
