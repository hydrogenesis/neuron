#pragma once
#include <vector>
#include <string>
#include "hh.h"

class Neuromorphic {
public:
  static void StringSplit(const std::string& s, char delim, std::vector<std::string>* result);

  static void LoadSwc(const char* filename, std::vector<HH*>* neuron, int duplicate);
  static void RandomSynapse(std::vector<HH*>* neuron, double factor, double I);
};  // class Neuromorphic
