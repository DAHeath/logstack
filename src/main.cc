#include "scheme.h"
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



Garbling garble(const PRF& f, const Circuit& c, PRG& prg) {
  auto inpEnc = genEncoding(prg, c.nInp);
  Material m(c.nRow);
  auto outEnc = gb(f, c, inpEnc, m);
  return { m, { inpEnc, outEnc } };
}


constexpr std::size_t n_repetitions = 1;

double experiment(const Circuit& c) {
  std::vector<double> results;
  for (std::size_t i = 0; i < n_repetitions; ++i) {
    results.push_back(timed([&] {
      PRG seed;
      PRF k;
      auto [material, interface] = garble(k, c, seed);


      const auto delta1 = interface.inpEnc.delta;
      const auto delta2 = interface.outEnc.delta;

      Labelling inp(c.nInp);
      for (std::size_t i = 0; i < c.nInp; ++i) {
        inp[i] = interface.inpEnc.zeros[i];
      }

      std::cout << "EV\n";

      const auto out = ev(k, c, inp, material);


      for (std::size_t i = 0; i < c.nOut; ++i) {
        if (out[i] == interface.outEnc.zeros[i]) {
          std::cout << '0';
        } else if (out[i] == (interface.outEnc.zeros[i] ^ delta2)) {
          std::cout << '1';
        } else {
          /* std::cerr << "ERROR!\n"; */
          /* std::cerr << out[i] << '\n'; */
          /* std::cerr << g.outputEncoding.zeros[i] << '\n'; */
          /* std::cerr << (g.outputEncoding.zeros[i] ^ delta2) << '\n'; */
          /* std::exit(1); */
        }
      }
      std::cout << '\n';
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
  std::vector<Circuit> cs = {
    sha_netlist, sha_netlist,
    sha_netlist, sha_netlist,
    sha_netlist, sha_netlist,
    sha_netlist, sha_netlist,
  };
  std::cout << experiment(conditional(cs)) << '\n';


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


  /* std::vector<Circuit> cs; */
  /* for (std::size_t i = 1; i <= 16; ++i) { */
  /*   cs.push_back(sha_netlist); */
  /*   if (cs.size() == 1) { */
  /*     std::cout << experiment(cs[0]) << '\n'; */
  /*   } else { */
  /*     std::cout << experiment(conditional(cs)) << '\n'; */
  /*   } */
  /* } */


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
