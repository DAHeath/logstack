#include "scheme.h"
#include <iostream>


int main() {

  {
  PRG prg;
  PRF f;

  auto e = genEncoding(prg, 2);
  auto e0 = genEncoding(prg, 1);
  auto e1 = genEncoding(prg, 1);

  std::span<Label> zeros { e.zeros };
  zeros = zeros.subspan(1);

  Material material(3*1 + 2);
  std::span<Label> mat { material };
  auto [b0, b1] = gbDem(prg, f, e.delta, e.zeros[0], zeros, e0, e1, mat);

  mat = material;
  /* zeros[0] ^= e.delta; */
  const auto [X, Y] = evDem(f, e.zeros[0] ^ e.delta, zeros, mat);

  const auto good00 = e0.zeros[0];
  const auto good10 = e1.zeros[0];
  const auto good01 = e0.zeros[0] ^ e0.delta;
  const auto good11 = e1.zeros[0] ^ e1.delta;
  const auto bad0 = b0[0];
  const auto bad1 = b1[0];
  std::cout << "X: " << (X[0] == good00) << (X[0] == good01) << (X[0] == bad0) << '\n';
  std::cout << "Y: " << (Y[0] == good10) << (Y[0] == good11) << (Y[0] == bad1) << '\n';
  }

  /* std::cout << "good0: " << good << '\n'; */
  /* std::cout << "good1: " << (good ^ delta) << '\n'; */
  /* std::cout << "bad:   " << bad << '\n'; */
  /* std::cout << "X:     " << X << '\n'; */
  /* std::cout << "Y:     " << Y << '\n'; */


/*   Circuit c { */
/*     Netlist { */
/*       Gate { GateType::INPUT, 0, 0, 0 }, */
/*       Gate { GateType::INPUT, 0, 0, 1 }, */
/*       Gate { GateType::AND, 0, 1, 2 }, */
/*       Gate { GateType::OUTPUT, 2, 0 , 0 }, */
/*     }, */
/*     3, // nWires */
/*     2, // nInp */
/*     1, // nOut */
/*     2, // nRow */
/*   }; */


/*   PRG seed; */
/*   PRF k; */

/*   auto g = garble(seed, k, c); */

/*   const auto delta = g.inputEncoding.delta; */

/*   const Labelling inp = { */
/*     g.inputEncoding.zeros[0] ^ delta, */
/*     g.inputEncoding.zeros[1] ^ delta, */
/*   }; */
/*   std::span<Label> mat(g.material); */
/*   const auto out = ev(k, c, inp, mat); */


/*   std::cout << (out[0] == g.outputEncoding.zeros[0]) << '\n'; */
/*   std::cout << (out[0] == (g.outputEncoding.zeros[0] ^ delta)) << '\n'; */
}
