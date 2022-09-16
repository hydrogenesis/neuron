#pragma once
#include <vector>
#include <string>
#include "hh.h"

class Neuromorphic {
public:
  static void StringSplit(const std::string& s, char delim, std::vector<std::string>* result);

  static long long GetNumExpandedNeurons(const char* filename, int duplicate, double dx);
  static void LoadSwc(const char* filename, std::vector<HH*>* neuron, int duplicate, double dx);
  static void RandomSynapse(std::vector<HH*>* neuron, double factor, double I);
};  // class Neuromorphic
