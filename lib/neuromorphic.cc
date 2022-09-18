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

void LoadFile(const char* filename, std::vector<std::vector<std::string>>* content) {
  std::ifstream ifs(filename);
  for (std::string line; std::getline(ifs, line);) {
    if (line.size() == 0) continue;
    if (line[0] == '#') continue;
    std::vector<std::string> ss;
    Neuromorphic::StringSplit(line, ' ', &ss);
    // std::cout << "ss size: " << ss.size() << std::endl;
    if (ss.size() != 7) continue;
    content->push_back(ss);
  }
  ifs.close();
}

struct Point {
  float x;
  float y;
  float z;
};
long long Neuromorphic::GetNumExpandedNeurons(const char* filename, int duplicate, double dx) {
  std::vector<std::vector<std::string>> content;
  LoadFile(filename, &content);
  std::vector<Point> mapping(content.size() * duplicate);
  long long total_num = 0;
  for (auto it = content.begin(); it != content.end(); ++it) {
    // read one line
    int index = atoi((*it)[0].c_str());
    int type = atoi((*it)[1].c_str());
    float x = atof((*it)[2].c_str());
    float y = atof((*it)[3].c_str());
    float z = atof((*it)[4].c_str());
    // float radius = atof((*it)[5].c_str());
    int parent = atoi((*it)[6].c_str());
    for (int i = 0; i < duplicate; ++i) {
      int the_index = i * content.size() + index - 1;
      mapping[the_index].x = x;
      mapping[the_index].y = y;
      mapping[the_index].z = z;
      if (parent == -1 || type == 1) {  // soma as sphere
        total_num++;
      } else {  // dendrite or axon as cylinder
        int parent_index = i * content.size() + parent - 1;  // parent
        auto p = mapping[parent_index];
        double length = std::sqrt((p.x - x)*(p.x - x) + (p.y - y)*(p.y - y) + (p.z - z)*(p.z - z));
        int num_compartments = std::ceil(length / dx);
        total_num += num_compartments;
      }
    }
  }
  return total_num;
}

void InsertCylinder(int the_index, int parent_index,
    double x, double y, double z, int type, double radius,
    std::vector<long long>& heads,
    std::vector<long long>& tails,
    long long* current_neuron,
    std::vector<HH*>* neuron,
    double dx) {
  HH* p = (*neuron)[heads[parent_index]];  // parent
  double length = std::sqrt((p->x - x)*(p->x - x) + (p->y - y)*(p->y - y) + (p->z - z)*(p->z - z));
  int n = int(std::ceil(length / dx));
  // first compartment
  HH* compartment = (*neuron)[*current_neuron];
  heads[the_index] = *current_neuron;
  double compartment_length = dx;
  if (length < dx) {
    // if the head is the tail
    compartment_length = length;
    compartment_length = dx;
    tails[the_index] = *current_neuron;
  }
  (*current_neuron)++;
  compartment->SetParams(x, y, z, type, compartment_length, radius);
  if (parent_index != -1) {
    (*neuron)[tails[parent_index]]->Append(compartment);
  }
  for (int j = 1; j < n; ++j) {
    // following compartments
    compartment_length = dx;
    if (length < (j + 1) * dx) {
      compartment_length = length - j * dx;
      compartment_length = dx;
    }
    compartment = (*neuron)[*current_neuron];
    // locate tail
    if (j == n - 1) tails[the_index] = *current_neuron;
    compartment->SetParams(x, y, z, type, compartment_length, radius);
    (*neuron)[*current_neuron - 1]->Append(compartment);
    (*current_neuron)++;
  }
}

void Neuromorphic::LoadSwc(const char* filename, std::vector<HH*>* neuron, int duplicate, double dx) {
  std::vector<std::vector<std::string>> content;
  LoadFile(filename, &content);
  int original_neuron_size = content.size() * duplicate;
  std::vector<long long> heads(original_neuron_size);
  std::vector<long long> tails(original_neuron_size);
  long long current_neuron = 0;
  long long neuron_size = Neuromorphic::GetNumExpandedNeurons(filename, duplicate, dx);
  std::cout << "neuron size: " << content.size() << " compartments, expanded: " << neuron_size << std::endl;

  neuron->resize(neuron_size);
  for (long long i = 0; i < neuron_size; ++i) {
    (*neuron)[i] = new HH();
    (*neuron)[i]->index = i;
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
      int the_index = i * content.size() + index - 1;
      if (parent == -1) {  // first compartment as sphere
        if (type == 1) {
          // first compartment is sphere
          HH* compartment = (*neuron)[current_neuron];
          heads[the_index] = current_neuron;
          tails[the_index] = current_neuron;
          current_neuron++;
          compartment->SetParams(x, y, z, type, radius, radius);
        } else {
          // first compartment is cylinder
          std::cout << "ERROR!" << std::endl;
          // InsertCylinder(the_index, -1,
          //     x, y, z, type, radius,
          //     heads,
          //     tails,
          //     &current_neuron,
          //     neuron,
          //     dx);
        }
      } else {
        if (type == 1) {
          // soma as sphere
          HH* compartment = (*neuron)[current_neuron];
          heads[the_index] = current_neuron;
          tails[the_index] = current_neuron;
          int parent_index = i * content.size() + parent - 1;
          (*neuron)[tails[parent_index]]->Append(compartment);
          compartment->SetParams(x, y, z, type, radius, radius);
          current_neuron++;
        } else {  // dendrite or axon as cylinder
          int parent_index = i * content.size() + parent - 1;
          InsertCylinder(the_index, parent_index,
              x, y, z, type, radius,
              heads,
              tails,
              &current_neuron,
              neuron,
              dx);
        }
      }
      // std::cout << "area: " << compartment->area << std::endl;
      // compartment->g = compartment->area / 2.0;
    }
  }
  // connect axon end to next soma
  // for (int i = 1; i < duplicate; ++i) {
  //   (*neuron)[(i-1)*content.size() + content.size() - 1]->Append((*neuron)[i * content.size()]);
  // }
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
