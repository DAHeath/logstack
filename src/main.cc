#include "scheme.h"
#include <iostream>


int main() {

  /* { */
  /*   PRG prg; */
  /*   PRF f; */

  /*   auto e = genEncoding(prg, 2); */
  /*   auto e0 = genEncoding(prg, 1); */
  /*   auto e1 = genEncoding(prg, 1); */

  /*   std::span<Label> zeros { e.zeros }; */
  /*   zeros = zeros.subspan(1); */

  /*   Material material(3*1 + 2); */
  /*   std::span<Label> mat { material }; */
  /*   auto [b0, b1] = gbDem(prg, f, e.delta, e.zeros[0], zeros, e0, e1, mat); */

  /*   mat = material; */
  /*   /1* zeros[0] ^= e.delta; *1/ */
  /*   const auto [X, Y] = evDem(f, e.zeros[0], zeros, mat); */

  /*   const auto good00 = e0.zeros[0]; */
  /*   const auto good10 = e1.zeros[0]; */
  /*   const auto good01 = e0.zeros[0] ^ e0.delta; */
  /*   const auto good11 = e1.zeros[0] ^ e1.delta; */
  /*   const auto bad0 = b0[0]; */
  /*   const auto bad1 = b1[0]; */
  /*   std::cout << "X: " << (X[0] == good00) << (X[0] == good01) << (X[0] == bad0) << '\n'; */
  /*   std::cout << "Y: " << (Y[0] == good10) << (Y[0] == good11) << (Y[0] == bad1) << '\n'; */
  /* } */

  /* { // test mux */

  /*   PRG prg; */
  /*   PRF f; */

  /*   const auto e = genEncoding(prg, 1); */
  /*   const auto e0 = genEncoding(prg, 1); */
  /*   const auto e1 = genEncoding(prg, 1); */

  /*   const auto b0 = prg(); */
  /*   const auto b1 = prg(); */


  /*   Material material(1 + 2); */

  /*   const auto eout = gbMux(prg, f, e.delta, e.zeros[0], e0, e1, { b0 }, { b1 }, material); */

  /*   /1* const auto out = evMux(f, e.zeros[0], { e0.zeros[0]}, { b1 }, material); *1/ */
  /*   const auto out = evMux(f, e.zeros[0] ^ e.delta, { b0 }, { e1.zeros[0] ^ e1.delta }, material); */

  /*   std::cout << (out[0] == eout.zeros[0]) << '\n'; */
  /*   std::cout << (out[0] == (eout.zeros[0] ^ eout.delta)) << '\n'; */
  /* } */

  {

  Circuit andc {
    Netlist {
      Gate { GateType::INPUT, 0, 0, 0 },
      Gate { GateType::INPUT, 0, 0, 1 },
      Gate { GateType::AND, 0, 1, 2 },
      Gate { GateType::OUTPUT, 2, 0 , 0 },
    },
    2, // nInp
    1, // nOut
    2, // nRow
  };
  Circuit xorc {
    Netlist {
      Gate { GateType::INPUT, 0, 0, 0 },
      Gate { GateType::INPUT, 0, 0, 1 },
      Gate { GateType::XOR, 0, 1, 2 },
      Gate { GateType::OUTPUT, 2, 0 , 0 },
    },
    2, // nInp
    1, // nOut
    0, // nRow
  };


  Circuit c {
    Conditional { { andc, xorc } },
    3, // nInp
    1, // nOut
    20
  };


  PRG seed;
  PRF k;

  Material material(c.nRow);
  auto g = garble(seed, k, c, material);

  const auto delta1 = g.inputEncoding.delta;
  const auto delta2 = g.outputEncoding.delta;

  const Labelling inp = {
    g.inputEncoding.zeros[0],
    g.inputEncoding.zeros[1] ^ delta1,
    g.inputEncoding.zeros[2] ^ delta1,
  };
  const auto out = ev(k, c, inp, material);

  std::cout << (out[0] == g.outputEncoding.zeros[0]) << '\n';
  std::cout << (out[0] == (g.outputEncoding.zeros[0] ^ delta2)) << '\n';

  }
}
