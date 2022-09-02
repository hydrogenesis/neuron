#include "hh.h"
#include <cmath>
#include <stdio.h>

// HH::HH():HH(0.0, 25, 0.025, -65, 0.5, 0.06, 0.5,
HH::HH():HH(0.0, 25, 0.025, -65, 0.05293248525724958, 0.5961207535084603, 0.3176769140606974,
     120.0, 115.0, 36.0, -12.0, 0.3, 10.6, 3) {
}

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
  t = 0.0;
  V = v;
  m = mi;
  h = hi;
  n = ni;
  t1 = 0.0;
  V1 = 0.0;
  m1 = 0.0;
  h1 = 0.0;
  n1 = 0.0;
}

HH::~HH() {

}

double HH::DynamicI(double t, double I) const {
  // return I*std::sin(t);
  if (5 <= t && t < 6) return I;
  else
    return 0;
  if ((int)(t / 100) % 2 == 0) return I;
  else
    return 0.0;
  // if (t < 50.0) return I;
  // else
  //   return 0.0;
}

void HH::Append(HH* next_hh) {
  next.push_back(next_hh);
  next_hh->prev.push_back(this);
}

void HH::Advance() {
  ++i;
  V = V1;
  m = m1;
  h = h1;
  n = n1;
  t = dt * i;
}

void HH::Advance_Euler() {
  V1 = V + dt*(
      gNa*m*m*m*h*(eNa-(V+65))
      + gK*n*n*n*n*(eK-(V+65))
      + gL*(eL-(V+65))
      + DynamicI(i*dt, I)
      );
  m1 = m + dt*(AlphaM(V)*(1-m)-BetaM(V)*m);
  h1 = h + dt*(AlphaH(V)*(1-h)-BetaH(V)*h);
  n1 = n + dt*(AlphaN(V)*(1-n)-BetaN(V)*n);
  for (auto hh : next) {
    V1 += dt * (hh->V - V) * hh->g;
  }
  for (auto hh : prev) {
    V1 -= dt * (V - hh->V) * hh->g;
  }
  if (m1 > 1) m1 = 0.999;
  if (m1 <= 0) m1 = 0.001;
  if (n1 > 1) n1 = 0.999;
  if (n1 <= 0) n1 = 0.001;
  if (h1 > 1) h1 = 0.999;
  if (h1 <= 0) h1 = 0.001;
  if (V1 >= 100) V1 = -65;
  if (V1 <= -200) V1 = -65;
  // if (V1 >= 200 || V1 <= -200) {
  //   printf("V : %.06f, %.06f, %.06f, %.06f\n", V, m, h, n);
  //   printf("V1: %.06f, %.06f, %.06f, %.06f\n", V1, m1, h1, n1);
  //   printf("ab: %.06f, %.06f, %.06f, %.06f, %.06f, %.06f\n",
  //       gNa*(eNa-(V+65)), gK*(eK-(V+65)), gL*(eL-(V+65)),
  //       BetaM(V), BetaH(V), BetaN(V)
  //       );
  // }
}

void HH::Advance_Crank_Nicolson(int iter_num) {
  for (int it = 0; it < iter_num; ++it) {
    V1 = V + 0.5*dt*(
        gNa*m1*m1*m1*h1*(eNa-(V1+65))
        + gK*n1*n1*n1*n1*(eK-(V1+65)) 
        + gL*(eL-(V1+65))
        + DynamicI((i+1)*dt, I)
        + gNa*m*m*m*h*(eNa-(V+65))
        + gK*n*n*n*n*(eK-(V+65))
        + gL*(eL-(V+65))
        + DynamicI(i*dt, I)
        );
    m1 = m + 0.5*dt*(AlphaM(V1)*(1-m1)-BetaM(V1)*m1 + AlphaM(V)*(1-m)-BetaM(V)*m);
    h1 = h + 0.5*dt*(AlphaH(V1)*(1-h1)-BetaH(V1)*h1 + AlphaH(V)*(1-h)-BetaH(V)*h);
    n1 = n + 0.5*dt*(AlphaN(V1)*(1-n1)-BetaN(V1)*n1 + AlphaN(V)*(1-n)-BetaN(V)*n);
    for (auto hh : next) {
      V1 += 0.5 * dt * (hh->V1 - V1 + hh->V - V) * hh->g;
    }
    for (auto hh : prev) {
      V1 -= 0.5 * dt * (V1 - hh->V1 + V - hh->V) * hh->g;
    }
  }
}

