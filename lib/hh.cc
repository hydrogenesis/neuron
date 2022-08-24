#include "hh.h"
#include <cmath>
#include <stdio.h>

HH::HH(double I, double tspan, double dt, double v, double mi, double hi, double ni,
      double gNa, double eNa, double gK, double eK, double gL, double eL, double g) {
  this->i = 0;
  this->I = I;
  this->dt = dt;
  this->loop = std::ceil(tspan / dt);
  this->v = v;
  this->mi = mi;
  this->hi = hi;
  this->ni = ni;
  this->gNa = gNa;
  this->eNa = eNa;
  this->gK = gK;
  this->eK = eK;
  this->gL = gL;
  this->eL = eL;
  this->g = g;
  t.resize(loop);
  V.resize(loop);
  m.resize(loop);
  h.resize(loop);
  n.resize(loop);

  for (long long k = 0; k < loop; k++) {
    t[k] = k * dt;
    V[k] = 0;
    m[k] = 0;
    h[k] = 0;
    n[k] = 0;
  }
  V[0] = v;
  m[0] = mi;
  h[0] = hi;
  n[0] = ni;
}

HH::~HH() {

}

double HH::DynamicI(double t, double I) const {
  return I;
}

void HH::Append(HH* next_hh) {
  next.push_back(next_hh);
  next_hh->prev.push_back(this);
}

void HH::Advance_Euler() {
  // V[i] = v;
  // m[i] = mi;
  // h[i] = hi;
  // n[i] = ni;
  V[i+1] = V[i] + dt*(
      gNa*m[i]*m[i]*m[i]*h[i]*(eNa-(V[i]+65))
      + gK*n[i]*n[i]*n[i]*n[i]*(eK-(V[i]+65))
      + gL*(eL-(V[i]+65))
      + DynamicI(i*dt, I)
      );
  m[i+1] = m[i] + dt*(AlphaM(V[i])*(1-m[i])-BetaM(V[i])*m[i]);
  h[i+1] = h[i] + dt*(AlphaH(V[i])*(1-h[i])-BetaH(V[i])*h[i]);
  n[i+1] = n[i] + dt*(AlphaN(V[i])*(1-n[i])-BetaN(V[i])*n[i]);
  for (auto hh : next) {
    V[i+1] += dt * (hh->V[hh->i] - V[i]) * hh->g;
  }
  for (auto hh : prev) {
    V[i+1] -= dt * (V[i] - hh->V[hh->i]) * hh->g;
  }
  // printf("%.06f, %.06f, %.06f, %.06f\n", V[i+1], m[i+1], h[i+1], n[i+1]);
}
