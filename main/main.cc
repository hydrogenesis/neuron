#include <stdio.h>
#include <vector>

#include "lib/hh.h"
#include "lib/neuromorphic.h"

int main(int argc, char** argv) {

  std::vector<HH*> neuron;
  if (argc == 1) {
    // default setup: an axon
    int branch = 50;
    for (int i = 0; i < branch; ++i) {
      double I = 0.0;
      if (i == 0) I = 50.0;
      auto membrane = new HH(I, 200, 0.001, -65, 0.5, 0.06, 0.5,
        120.0, 115.0, 36.0, -12.0, 0.3, 10.6, 3.0);
      neuron.push_back(membrane);
      if (i > 0) {
        neuron[i - 1]->Append(neuron[i]);
      }
    }
  } else {
    // swc setup
    Neuromorphic::LoadSwc(argv[1], &neuron);
    // neuron[0]->I = 200.0;
    neuron[136]->I = 40.0;
    neuron[140]->I = 40.0;
  }

  for (long long i = 0; i < neuron[0]->loop; ++i) {
    if (i % 1000 == 0) {
      printf("loop: %lld\n", i);
    }
    for (int i = 0; i < (int)neuron.size(); ++i) {
      neuron[i]->Advance_Euler();
    }
    for (int i = 0; i < (int)neuron.size(); ++i) {
      neuron[i]->Advance();
    }
  }

  auto f = fopen("hh.txt", "w");
  for (long long i = 0; i < neuron[0]->loop; ++i) {
    fprintf(f, "%.6f\t%.6f\t%.6f\t%.6f\n", neuron[0]->t[i], neuron[0]->V[i], neuron[neuron.size() - 1]->t[i], neuron[neuron.size() - 1]->V[i]);
  }
  fclose(f);
  
  for (int i = 0; i < (int)neuron.size(); ++i) {
    delete neuron[i];
  }
  return 0;
}
