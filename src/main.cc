#include "scheme.h"
#include <iostream>

int main() {

  Circuit c {
    Netlist {
      Gate { GateType::INPUT, 0, 0, 0 },
      Gate { GateType::INPUT, 0, 0, 1 },
      Gate { GateType::AND, 0, 1, 2 },
      Gate { GateType::OUTPUT, 2, 0 , 0 },
    },
    3, // nWires
    2, // nInp
    1, // nOut
    2, // nRow
  };


  PRG seed;
  PRF k;

  auto g = garble(seed, k, c);

  const auto delta = g.inputEncoding.delta;

  const Labelling inp = {
    g.inputEncoding.zeros[0],
    g.inputEncoding.zeros[1],
  };
  std::span<Label> mat(g.material);
  const auto out = ev(k, c, inp, mat);


  std::cout << (out[0] == g.outputEncoding.zeros[0]) << '\n';
  std::cout << (out[0] == (g.outputEncoding.zeros[0] ^ delta)) << '\n';
}
