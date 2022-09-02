#include "neuromorphic.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

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
    int type = atoi((*it)[1].c_str());
    float x = atof((*it)[2].c_str());
    float y = atof((*it)[3].c_str());
    float z = atof((*it)[4].c_str());
    float radius = atof((*it)[5].c_str());
    int parent = atoi((*it)[6].c_str());
    for (int i = 0; i < duplicate; ++i) {
      HH* compartment = (*neuron)[i * content.size() + index - 1];
      compartment->x = x;
      compartment->y = y;
      compartment->z = z;
      compartment->type = type;
      compartment->radius = radius;
      if (parent == -1) {
        compartment->area = 4 * 3.1415927 * radius * radius;
      } else {
        HH* p = (*neuron)[i * content.size() + parent - 1];  // parent
        double length = std::sqrt((p->x - x)*(p->x - x) + (p->y - y)*(p->y - y) + (p->z - z)*(p->z - z));
        compartment->area = length * 2 * 3.1415927 * radius;
        p->Append(compartment);
      }
      // std::cout << "area: " << compartment->area << std::endl;
      compartment->gNa = 0.12 * compartment->area;
      compartment->gK = 0.036 * compartment->area;
      compartment->gL = 0.0003 * compartment->area;
    }
  }
  // connect axon end to next soma
  for (int i = 1; i < duplicate; ++i) {
    (*neuron)[(i-1)*content.size() + content.size() - 1]->Append((*neuron)[i * content.size()]);
  }
}

void Neuromorphic::RandomSynapse(std::vector<HH*>* neuron, double factor, double I) {
  long long neuron_num = neuron->size();
  long long seed = time(NULL);
  // seed = 1661933470;
  printf("seed = %lld\n", seed);
  srand(seed);
  for (long long i = 0; i < neuron_num; ++i) {
    int type = (*neuron)[i]->type;
    if (type != 1 && type != 2) continue;
    double r = (rand() % 1000) / 1000.0;
    if (r < factor) {
      (*neuron)[i]->I = I;
      // printf("r=%f, neuron %lld with type %d has synapse of %f.\n", r, i, type, I);
    }
  }
}
