#include "scheme.h"
#include "protocol.h"
#include "sha256.h"
#include <iostream>
#include <chrono>
#include <numeric>


template <typename F>
auto timed(F f) {
  auto start = std::chrono::high_resolution_clock::now();
  f();
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  return elapsed.count();
}



constexpr std::size_t n_repetitions = 1;


emp::NetIO* genChannel;
emp::NetIO* evalChannel;

double experiment(const Circuit& c) {
  std::vector<double> results;
  for (std::size_t i = 0; i < n_repetitions; ++i) {
    results.push_back(timed([&] {
      semihonest(c, *genChannel, *evalChannel);
    }));
  }

  return std::accumulate(results.begin(), results.end(), 0.0)/results.size();
}


int main(int argc, char** argv) {
  const auto sha_netlist = Circuit {
    sha,
    sha.desc.nInp,
    sha.desc.nOut,
    sha.desc.nRow,
  };

  std::thread th { [&] {
    genChannel = new emp::NetIO(nullptr, 44444, true);
  }};
  sleep(1);
  evalChannel = new emp::NetIO("127.0.0.1", 44444, true);
  th.join();

  std::vector<Circuit> cs;
  for (std::size_t j = 0; j < (1 << 4); ++j) {
    cs.push_back(sha_netlist);
  }
  if (cs.size() == 1) {
    std::cout << experiment(cs[0]) << '\n';
  } else {
    std::cout << experiment(conditional(cs)) << '\n';
  }
}
