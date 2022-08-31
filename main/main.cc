#include <stdio.h>
#include <vector>
#include <iostream>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <unistd.h>

#include "lib/hh.h"
#include "lib/neuromorphic.h"
#include "boost/lockfree/queue.hpp"

#define CHUNK_SIZE 1000
#define DUPLICATE 1

void worker(int id, int total_threads, std::vector<HH*>* neuron, std::mutex* m, std::condition_variable* cv, int* signal, std::atomic<int>* p) {
  long long total_neurons = neuron->size();
  long long chunk_size = total_neurons / total_threads;
  if (chunk_size * total_threads < total_neurons) chunk_size++;
  while (1) {
    int s = *signal;
    if (s == 0) {
      for (long long i = id * chunk_size; i < (id + 1) * chunk_size && i < total_neurons; ++i) {
        (*neuron)[i]->Advance_Euler();
      }
    } else if (s == 1) {
      break;
    } else {
      // std::this_thread::sleep_for(std::chrono::milliseconds(1));
      usleep(1);
      continue;
    }
    {
      std::unique_lock<std::mutex> lk(*m);
      // std::cout << id << " " << s << " " << *p << std::endl;
      (*p)++;
      cv->wait(lk);
    }
  }
}

void worker_queue(int id, int total_threads, std::vector<HH*>* neuron, boost::lockfree::queue<long long>* queue, int* signal, std::atomic<int>* p) {
  long long total_neurons = neuron->size();
  long long chunk_id = 0;
  long long chunk_size = CHUNK_SIZE;
  while (1) {
    if (*signal == 3) break;
    bool ret = queue->pop(chunk_id);
    if (!ret) {
      // std::this_thread::sleep_for(std::chrono::milliseconds(1));
      usleep(1);
      continue;
    }

    if (*signal == 0) {
      for (long long i = chunk_id * chunk_size; i < (chunk_id + 1) * chunk_size && i < total_neurons; ++i) {
        (*neuron)[i]->Advance_Euler();
      }
    } else if (*signal == 1) {
      for (long long i = chunk_id * chunk_size; i < (chunk_id + 1) * chunk_size && i < total_neurons; ++i) {
        (*neuron)[i]->Advance_Crank_Nicolson(3);
      }
    } else if (*signal == 2) {
      for (long long i = chunk_id * chunk_size; i < (chunk_id + 1) * chunk_size && i < total_neurons; ++i) {
        (*neuron)[i]->Advance();
      }
    }
    (*p)++;
  }
}

int main(int argc, char** argv) {
  std::vector<HH*> neuron;
  int num_threads = 0;
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
    num_threads = atoi(argv[1]);
    // swc setup
    Neuromorphic::LoadSwc(argv[2], &neuron, DUPLICATE);
    // neuron[0]->I = 100.0;
    // neuron[136]->I = 100.0;
    // neuron[200]->I = 100.0;
    // neuron[12221]->I = 100.0;
    Neuromorphic::RandomSynapse(&neuron, 0.001, 100.0);
  }
  long long loop = neuron[0]->loop;
  vector<double> t1, t2, V1, V2;
  t1.reserve(loop);
  t2.reserve(loop);
  V1.reserve(loop);
  V2.reserve(loop);
  neuron[0]->Record(&t1, &V1);
  neuron[neuron.size() - 1]->Record(&t2, &V2);
  // threading vars
  std::mutex m;
  std::condition_variable cv;
  int signal = 2;
  std::atomic<int> p(0);
  vector<std::thread*> threads(num_threads);

  long long total_neurons = neuron.size();
  long long chunk_num = total_neurons / CHUNK_SIZE;
  if (CHUNK_SIZE * chunk_num < total_neurons) chunk_num++;
  boost::lockfree::queue<long long> queue(chunk_num);
  signal = 0;

  for (int i = 0; i < num_threads; ++i) {
    // threads[i] = new std::thread(worker, i, num_threads, &neuron, &m, &cv, &signal, &p);
    threads[i] = new std::thread(worker_queue, i, num_threads, &neuron, &queue, &signal, &p);
  }
  p = 0;
  signal = 0;
  for (long long i = 0; i < loop; ++i) {
    if (i % 1000 == 0) {
      printf("loop: %lld\n", i);
    }
    if (num_threads == 0) {
      for (int i = 0; i < (int)neuron.size(); ++i) {
        neuron[i]->Advance_Euler();
      }
      for (int j = 0; j < (int)neuron.size(); ++j) {
        neuron[j]->Advance();
      }
    } else {
      signal = 0;
      p = 0;
      for (long long chunk_id = 0; chunk_id < chunk_num; ++chunk_id) {
        queue.push(chunk_id);
      }
      while (p != chunk_num) {
        usleep(1);
      }
      signal = 2;
      p = 0;
      for (long long chunk_id = 0; chunk_id < chunk_num; ++chunk_id) {
        queue.push(chunk_id);
      }
      while (p != chunk_num) {
        usleep(1);
      }
      // join all threads
      // while (p != num_threads) {
      //   usleep(1);
      // }
    }
    neuron[0]->Record(&t1, &V1);
    neuron[neuron.size() - 1]->Record(&t2, &V2);
    // neuron[12222]->Record(&t2, &V2);
    // if (num_threads > 0) {
    //   p = 0;
    //   signal = 0;
    //   cv.notify_all();
    // }
  }
  if (num_threads > 0) {
    signal = 3;
    for (int j = 0; j < num_threads; ++j) {
      threads[j]->join();
    }
    // while (p != num_threads) {
    //   std::this_thread::sleep_for(std::chrono::microseconds(1));
    //   usleep(1);
    // }
    // signal = 1;
    // p = 0;
    // cv.notify_all();
    // for (int j = 0; j < num_threads; ++j) {
    //   threads[j]->join();
    // }
  }

  auto f = fopen("hh.txt", "w");
  for (long long i = 0; i < loop; ++i) {
    fprintf(f, "%.6f\t%.6f\t%.6f\t%.6f\n", t1[i], V1[i], t2[i], V2[i]);
  }
  fclose(f);
  
  for (int i = 0; i < (int)neuron.size(); ++i) {
    delete neuron[i];
  }
  return 0;
}
