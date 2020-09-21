#include "scheme.h"
#include <iostream>


std::tuple<Label, Label, Label> gbDem(
    const PRF& f,
    const Label& delta,
    const Label& S0,
    const Label& A0) {

  const auto s = S0[0];
  const auto a = A0[0];

  const auto hA0 = f(A0);
  const auto hA1 = f(A0 ^ delta);

  const auto row0 = hA0 ^ ((a != s) ? A0 : 0);
  const auto row1 = hA1 ^ ((a != s) ? delta : A0);

  const auto good = a ? row1 : row0;
  const auto bad = good ^ A0;
  const auto mat = row0 ^ row1;

  return { good, bad, mat };
}


std::pair<Label, Label> evDem(
    const PRF& f,
    const Label& S,
    const Label& A,
    const Label& mat) {
  const auto s = S[0];
  const auto a = A[0];

  const auto hA = f(A);

  const auto r0 = hA ^ (a ? mat : 0) ^ ((a != s) ? A : 0);
  const auto r1 = hA ^ (a ? mat : 0) ^ ((a != s) ? 0 : A);

  return { r0, r1 };
}


int main() {

  PRG prg;
  const auto A0 = prg();
  const auto S0 = prg();
  auto delta = prg();
  delta[0] = 1;

  PRF f;

  const auto [good, bad, mat] = gbDem(f, delta, S0, A0);

  const auto [X, Y] = evDem(f, S0 ^ delta, A0, mat);

  const auto good0 = good;
  const auto good1 = good ^ delta;

  std::cout << "X: " << (X == good0) << (X == good1) << (X == bad) << '\n';
  std::cout << "Y: " << (Y == good0) << (Y == good1) << (Y == bad) << '\n';

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
