#include "scheme.h"
#include "compiler.h"
#include "sha256.h"
#include <iostream>


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


int main() {

  {
    std::array<U32, 16> inp;
    for (auto& i: inp) { i = U32::input(); }

    const auto out = sha256(inp);

    for (const auto& o: out) { o.output(); }
  }

  /* { */
  /*   U32 x = U32::input(); */
  /*   U32 y = U32::input(); */

  /*   (x + y).output(); */
  /* } */

  const auto c = Bool::compile();
  /* showCircuit(c); */



  /* Circuit c { */
  /*   Conditional { { andc, andc, xorc, nand } }, */
  /*   4, // nInp */
  /*   1, // nOut */
  /*   40 */
  /* }; */

  for (std::size_t i = 0; i < 100; ++i) {
  PRG seed;
  PRF k;

  Material material(c.nRow);
  auto g = garble(seed, k, c, material);

  const auto delta1 = g.inputEncoding.delta;
  const auto delta2 = g.outputEncoding.delta;

  Labelling inp(c.nInp);
  for (std::size_t i = 0; i < c.nInp; ++i) {
    inp[i] = g.inputEncoding.zeros[i];
  }

  const auto out = ev(k, c, inp, material);

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
  }
}
