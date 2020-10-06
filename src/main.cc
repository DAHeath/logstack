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

  /* std::cout << (genChannel->counter + evalChannel->counter) << '\n'; */
  /* genChannel->counter = 0; */
  /* evalChannel->counter = 0; */

  return std::accumulate(results.begin(), results.end(), 0.0)/results.size();
}


int main(int argc, char** argv) {
  const auto sha_netlist = Circuit {
    sha,
    sha.desc.nInp,
    sha.desc.nOut,
    sha.desc.nRow,
  };
  /* std::vector<Circuit> cs = { */
  /*   sha_netlist, sha_netlist, */
  /*   sha_netlist, sha_netlist, */
  /*   sha_netlist, sha_netlist, */
  /*   sha_netlist, sha_netlist, */
  /* }; */
  /* std::cout << experiment(conditional(cs)) << '\n'; */


  /* std::vector<Circuit> cs = { sha_netlist }; */
  /* for (std::size_t i = 0; i <= 10; ++i) { */
  /*   if (cs.size() == 1) { */
  /*     std::cout << experiment(cs[0]) << '\n'; */
  /*   } else { */
  /*     std::cout << experiment(conditional(cs)) << '\n'; */
  /*   } */

  /*   for (std::size_t j = 0; j < (1 << i); ++j) { */
  /*     cs.push_back(sha_netlist); */
  /*   } */
  /* } */


  std::thread th { [&] {
    genChannel = new emp::NetIO(nullptr, 55556, true);
  }};
  sleep(1);
  evalChannel = new emp::NetIO("127.0.0.1", 55555, true);
  th.join();


  /* int n = atoi(argv[1]); */


  /* std::vector<Circuit> cs; */
  /* for (std::size_t i = 0; i < n; ++i) { */
  /*   cs.push_back(sha_netlist); */
  /* } */
  /* if (cs.size() == 1) { */
  /*   /1* std::cout << experiment(cs[0]) << '\n'; *1/ */
  /*   experiment(cs[0]); */
  /* } else { */
  /*   /1* std::cout << experiment(conditional(cs)) << '\n'; *1/ */
  /*   experiment(conditional(cs)); */
  /* } */



  /* /1* for (std::size_t i = 0; i < 128; ++i) { *1/ */
  /* /1*   cs.push_back(sha_netlist); *1/ */
  /* /1* } *1/ */
  /* /1* std::cout << experiment(conditional(cs)) << '\n'; *1/ */

  for (std::size_t i = 0; i <= 14; ++i) {
    std::vector<Circuit> cs;
    for (std::size_t j = 0; j < (1 << i); ++j) {
      cs.push_back(sha_netlist);
    }
    if (cs.size() == 1) {
      std::cout << experiment(cs[0]) << '\n';
      /* experiment(cs[0]); */
    } else {
      std::cout << experiment(conditional(cs)) << '\n';
      /* experiment(conditional(cs)); */
    }
  }


/*   PRG prg; */
/*   PRF f; */


/*   const auto enc = genEncoding(prg, 3); */

/*   const auto b = 2; */

/*   std::vector<Label> seeds(2*b - 2); */
/*   for (std::size_t i = 0; i < 2*b - 2; ++i) { */
/*     seeds[i] = prg(); */
/*   } */

/*   Material m(40); */
/*   for (auto& r: m) { r = 0; } */
/*   std::span<Label> mat(m); */
/*   const auto bad = gbGadget(b, f, enc.delta, seeds, enc.zeros, mat); */


/*   mat = m; */

/*   Labelling index = enc.zeros; */
/*   index[0] ^= enc.delta; */
/*   index[1] ^= enc.delta; */

/*   const auto actual = evGadget(b, f, index, mat); */


/*   for (std::size_t i = 0; i < seeds.size(); ++i) { */
/*     show(actual[i]); */
/*     show(seeds[i]); */
/*     show(bad[i]); */
/*   } */
}
