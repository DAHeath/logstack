#include "scheme.h"
#include "compiler.h"
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

void showGate(const Gate& g) {
  switch (g.type) {
    case GateType::INPUT: std::cout << "INPUT " << g.out; break;
    case GateType::OUTPUT: std::cout << "OUTPUT " << g.inp0; break;
    case GateType::AND: std::cout << "AND " << g.inp0 << " " << g.inp1 << " " << g.out; break;
    case GateType::XOR: std::cout << "XOR " << g.inp0 << " " << g.inp1 << " " << g.out; break;
    case GateType::NOT: std::cout << "NOT " << g.inp0 << " " << g.out; break;
  }
}


void showNetlist(const Netlist& nl) {
  for (const auto& g: nl) {
    showGate(g);
    std::cout << '\n';
  }
}


void showCircuit(const Circuit& c) {
  std::cout << "INPUTS: " << c.nInp << '\n';
  std::cout << "OUTPUTS: " << c.nOut << '\n';
  std::cout << "ROWS: " << c.nRow << '\n';
  std::visit(overloaded {
    [&](const Netlist& n) {
      showNetlist(n);
    },
    [&](const Conditional& cond) {
    },
    [&](const Sequence& seq) {
    },
  }, c.content);
}


Circuit conditional(const std::vector<Circuit>& cs) {
  const auto b = cs.size();
  const auto transSize = cs[0].nInp + 1;
  const auto demSize = 3*cs[0].nInp + 2;
  const auto muxSize = (b-1) * (cs[0].nOut + 2);

  const auto logb = ilog2(b);


  const auto fullDemSize = demSize * logb + 3*logb*logb;

  std::size_t nRow = 0;
  for (const auto& c: cs) { nRow = std::max(nRow, c.nRow); }
  nRow += transSize + fullDemSize + muxSize;

  return Circuit {
    Conditional { cs },
    cs[0].nInp + logb, // nInp
    cs[0].nOut, // nInp
    nRow, // nRow
  };
}


constexpr std::size_t n_repetitions = 1;

double experiment(const Circuit& c) {
  std::vector<double> results;
  for (std::size_t i = 0; i < n_repetitions; ++i) {
    results.push_back(timed([&] {
      PRG seed;
      PRF k;

      Material material(c.nRow);
      auto g = garble(seed, k, c, material);

      /* std::cout << g.inputEncoding.zeros.size() << '\n'; */

      const auto delta1 = g.inputEncoding.delta;
      const auto delta2 = g.outputEncoding.delta;

      Labelling inp(c.nInp);
      for (std::size_t i = 0; i < c.nInp; ++i) {
        inp[i] = g.inputEncoding.zeros[i];
      }

      ev(k, c, inp, material);
    }));
  }

  return std::accumulate(results.begin(), results.end(), 0.0)/results.size();
}


int main(int argc, char** argv) {
  {
    std::array<U32, 16> inp;
    for (auto& i: inp) { i = U32::input(); }

    const auto out = sha256(inp);

    for (const auto& o: out) { o.output(); }
  }

  const auto sha = Bool::compile();


  /* Circuit c; */
  /* if (n <= 1) { */
  /*   c = sha; */
  /* } else { */
  /*   std::vector<Circuit> cs; */
  /*   for (std::size_t i = 0; i < n; ++i) { */
  /*     cs.push_back(sha); */
  /*   } */
  /*   c = conditional(cs); */
  /* } */

  std::vector<Circuit> cs;
  for (std::size_t i = 1; i <= 64; ++i) {
    cs.push_back(sha);
    experiment(conditional(cs));
    /* std::cout << experiment(conditional(cs)) << '\n'; */
  }





  /* for (std::size_t i = 0; i < c.nOut; ++i) { */
  /*   if (out[i] == g.outputEncoding.zeros[i]) { */
  /*     std::cout << '0'; */
  /*   } else if (out[i] == (g.outputEncoding.zeros[i] ^ delta2)) { */
  /*     std::cout << '1'; */
  /*   } else { */
  /*     std::cerr << "ERROR!\n"; */
  /*     std::cerr << out[i] << '\n'; */
  /*     std::cerr << g.outputEncoding.zeros[i] << '\n'; */
  /*     std::cerr << (g.outputEncoding.zeros[i] ^ delta2) << '\n'; */
  /*     std::exit(1); */
  /*   } */
  /* } */
  /* std::cout << '\n'; */
}
