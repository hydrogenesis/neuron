#include "neuromorphic.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// trim from start (copying)
static inline std::string ltrim_copy(std::string s) {
    ltrim(s);
    return s;
}

// trim from end (copying)
static inline std::string rtrim_copy(std::string s) {
    rtrim(s);
    return s;
}

// trim from both ends (copying)
static inline std::string trim_copy(std::string s) {
    trim(s);
    return s;
}

void Neuromorphic::StringSplit(const std::string& s, char delim, std::vector<std::string>* result) {
  std::stringstream ss(trim_copy(s));
  std::string item;

  while (getline(ss, item, delim)) {
    result->push_back (item);
  }
}

void Neuromorphic::LoadSwc(const char* filename, std::vector<HH*>* neuron, int duplicate) {
  std::ifstream ifs(filename);
  std::vector<std::vector<std::string>> content;
  for (std::string line; std::getline(ifs, line);) {
    if (line.size() == 0) continue;
    if (line[0] == '#') continue;
    std::vector<std::string> ss;
    StringSplit(line, ' ', &ss);
    // std::cout << "ss size: " << ss.size() << std::endl;
    if (ss.size() != 7) continue;
    content.push_back(ss);
  }
  ifs.close();
  std::cout << "neuron size: " << content.size() << " compartments" << std::endl;
  long long neuron_size = content.size() * duplicate;

  neuron->resize(neuron_size);
  for (long long i = 0; i < neuron_size; ++i) {
    (*neuron)[i] = new HH();
  }
  std::cout << "duplicated neuron size: " << duplicate << " duplicates, " << neuron->size() << " compartments" << std::endl;

  for (auto it = content.begin(); it != content.end(); ++it) {
    // read one line
    int index = atoi((*it)[0].c_str());
    int parent = atoi((*it)[6].c_str());
    if (parent == -1) continue;
    for (int i = 0; i < duplicate; ++i) {
      (*neuron)[i * content.size() + parent - 1]->Append((*neuron)[i * content.size() + index - 1]);
    }
  }
  // connect axon end to next soma
  for (int i = 1; i < duplicate; ++i) {
    (*neuron)[i*content.size() + content.size() - 1]->Append((*neuron)[i * content.size()]);
  }
}

