#include "scheme.h"


template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;


Material operator^(const Material& m0, const Material& m1) {
  const auto n = m0.size();
  assert(m1.size() == n);

  Material out(n);

  for (std::size_t i = 0; i < n; ++i) {
    out[i] = m0[i] ^ m1[i];
  }
  return out;
}


Material& operator^=(Material& m0, const Material& m1) {
  const auto n = m1.size();
  assert(m0.size() >= n);

  for (std::size_t i = 0; i < n; ++i) {
    m0[i] ^= m1[i];
  }
  return m0;
}


Material& operator^=(Material& m0, std::span<Label> m1) {
  const auto n = m1.size();
  assert(m0.size() >= n);

  for (std::size_t i = 0; i < n; ++i) {
    m0[i] ^= m1[i];
  }
  return m0;
}


template <typename T>
std::vector<T> operator|(std::vector<T> m0, std::vector<T> m1) {
  const auto n0 = m0.size();
  const auto n1 = m1.size();
  std::vector<T> out(n0 + n1);

  for (std::size_t i = 0; i < n0; ++i) { out[i] = std::move(m0[i]); }
  for (std::size_t i = 0; i < n1; ++i) { out[i + n0] = std::move(m1[i]); }
  return out;
}


std::tuple<Label, Label, Label> gbDem1(
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


std::tuple<Material, Labelling, Labelling> gbDem(
    PRG& seed,
    const PRF& f,
    const Label& delta,
    const Label& S0,
    std::span<const Label> zeros,
    const Encoding& e0,
    const Encoding& e1) {

  const auto n = zeros.size();

  Material material(3*n + 2);
  std::span<Label> mat(material);

  const auto hS0 = f(S0);
  const auto hS1 = f(S0 ^ delta);


  auto dd0 = delta ^ e0.delta;
  auto dd1 = delta ^ e1.delta;
  // some extra random bits to mask the definitely zero lsb of delta ^ delta0/delta1
  const auto mask = seed();
  dd0[0] = mask[0];
  dd1[0] = mask[1];

  mat[0] = hS0 ^ dd0;
  mat[1] = hS1 ^ dd1;

  auto bad_diff0 = hS1 ^ mat[0];
  auto bad_diff1 = hS0 ^ mat[1];
  // we know that the lsb of delta ^ delta0/delta1 should be 0
  bad_diff0[0] = 0;
  bad_diff1[0] = 0;

  mat = mat.subspan(2);

  Labelling bad0(n);
  Labelling bad1(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto [g, b, m] = gbDem1(f, delta, S0, zeros[i]);
    mat[0] = m;
    mat = mat.subspan(1);
    const auto good0 = g ^ (g[0] ? (delta ^ e0.delta) : 0);
    const auto good1 = g ^ (g[0] ? (delta ^ e1.delta) : 0);
    bad0[i] = b ^ (b[0] ? bad_diff0 : 0);
    bad1[i] = b ^ (b[0] ? bad_diff1 : 0);
    mat[0] = hS0 ^ good0 ^ e0.zeros[i];
    mat[1] = hS1 ^ good1 ^ e1.zeros[i];
    bad0[i] ^= hS1 ^ mat[0];
    bad1[i] ^= hS0 ^ mat[1];
    mat = mat.subspan(2);
  }

  return { material, bad0, bad1 };
}


std::pair<Label, Label> evDem1(
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


std::pair<Labelling, Labelling> evDem(
    const PRF& f,
    const Label& S,
    std::span<const Label> inp,
    std::span<Label>& mat) {

  const auto n = inp.size();

  const auto hS = f(S);
  auto diff0 = hS ^ mat[0];
  auto diff1 = hS ^ mat[1];
  diff0[0] = 0;
  diff1[0] = 0;

  mat = mat.subspan(2);

  Labelling inp0(n);
  Labelling inp1(n);
  for (std::size_t i = 0; i < n; ++i) {
    auto [i0, i1] = evDem1(f, S, inp[i], mat[0]);
    inp0[i] = i0 ^ (i0[0] ? diff0 : 0);
    inp1[i] = i1 ^ (i1[0] ? diff1 : 0);
    mat = mat.subspan(1);

    inp0[i] ^= hS ^ mat[0];
    inp1[i] ^= hS ^ mat[1];
    mat = mat.subspan(2);
  }

  return { inp0, inp1 };
}




void gbGate(
    const PRF& f,
    const Gate& g,
    const Label& delta,
    NetlistCtxt& ctxt) {

  const auto& A0 = ctxt.w[g.inp0];
  const auto& B0 = ctxt.w[g.inp1];
  auto& C0 = ctxt.w[g.out];
  switch (g.type) {
    case GateType::INPUT:
      // read in the next label from the input labelling
      C0 = ctxt.inp[0];
      ctxt.inp = ctxt.inp.subspan(1);
      break;

    case GateType::OUTPUT:
      // write out to the end of the output labelling
      ctxt.out[0] = A0;
      ctxt.out = ctxt.out.subspan(1);
      break;

    case GateType::AND: {
      const auto nonce0 = Label { ctxt.nonce };
      const auto nonce1 = Label { ctxt.nonce + 1 };
      const auto hA0 = f(A0 ^ nonce0);
      const auto hA1 = f(A0 ^ delta ^ nonce0);
      const auto hB0 = f(B0 ^ nonce1);
      const auto hB1 = f(B0 ^ delta ^ nonce1);
      const auto A0D = A0 & delta;
      const auto B0D = B0 & delta;

      const auto X = A0[0] ? (hA1 ^ B0D) : hA0;
      const auto Y = B0[0] ? (hB1 ^ A0D) : hB0;

      ctxt.material[0] = hA0 ^ hA1 ^ B0D;
      ctxt.material[1] = hB0 ^ hB1 ^ A0D;

      ctxt.nonce += 2;
      ctxt.material = ctxt.material.subspan(2);

      C0 = A0&B0 ^ X ^ Y;
      break;
    }

    case GateType::XOR:
      C0 = A0 ^ B0;
      break;
  }
}


void evGate(
    const PRF& f,
    const Gate& g,
    NetlistCtxt& ctxt) {
  const auto& A = ctxt.w[g.inp0];
  const auto& B = ctxt.w[g.inp1];
  auto& C = ctxt.w[g.out];
  switch (g.type) {
    case GateType::INPUT:
      // read in the next label from the input labelling
      C = ctxt.inp[0];
      ctxt.inp = ctxt.inp.subspan(1);
      break;

    case GateType::OUTPUT:
      // write out to the end of the output labelling
      ctxt.out[0] = A;
      ctxt.out = ctxt.out.subspan(1);
      break;

    case GateType::AND: {
      const auto nonce0 = Label { ctxt.nonce };
      const auto nonce1 = Label { ctxt.nonce + 1 };
      const auto hA = f(A ^ nonce0);
      const auto hB = f(B ^ nonce1);
      const auto X = A[0] ? hA ^ ctxt.material[0] : hA;
      const auto Y = B[0] ? hB ^ ctxt.material[1] : hB;
      ctxt.material = ctxt.material.subspan(2);
      ctxt.nonce += 2;
      C = (A&B) ^ X ^ Y;
      break;
    }

    case GateType::XOR:
      C = A ^ B;
      break;
  }
}


CondGarbling_ gbCond_(const PRF& f, std::span<const Circuit> cs, const Label& seed) {
  const auto n = cs.size();
  PRG prg(seed);

  if (n == 1) {
    const auto g = garble(prg, f, cs[0]);
    return { std::move(g.inputEncoding), std::move(g.material) };
  } else {
    // Generate a fresh encoding.
    // All branches should have same number of inputs.
    auto e = genEncoding(prg, cs[0].nInp);
    const auto s0 = e.zeros[0];
    const auto s1 = e.zeros[0] ^ e.delta;

    // split the vector of circuits
    const std::span<const Circuit> cs0 = cs.subspan(0, n/2);
    const std::span<const Circuit> cs1 = cs.subspan(n/2);

    const auto g0 = gbCond_(f, cs0, s1);
    const auto g1 = gbCond_(f, cs1, s0);

    std::span<Label> zeros { e.zeros };
    zeros = zeros.subspan(1);
    auto [mdem, e0, e1] = gbDem(prg, f, e.delta, s0, zeros, g0.inputEncoding, g1.inputEncoding);

    auto stack = g0.material ^ g1.material;
    return { e, std::move(mdem) | std::move(stack) };
  }
}


Labelling evCond(
    const PRF& f,
    std::span<const Circuit> cs,
    const Labelling& inp,
    std::span<Label>& mat) {

  const auto n = cs.size();

  if (n == 1) {
    return ev(f, cs[0], inp, mat);
  } else {
    // split input into the branch condition and the remaining wires
    const auto s = inp[0];
    std::span<const Label> inp_(inp);
    inp_ = inp_.subspan(1);


    // split the vector of circuits
    const std::span<const Circuit> cs0 = cs.subspan(0, n/2);
    const std::span<const Circuit> cs1 = cs.subspan(n/2);

    // garble both subtrees

    auto g0 = gbCond_(f, cs0, s);
    auto g1 = gbCond_(f, cs1, s);

    const auto [inp0, inp1] = evDem(f, s, inp_, mat);

    // unstack!
    g0.material ^= mat;
    g1.material ^= mat;

    std::span mat0 { g1.material };
    std::span mat1 { g0.material };

    return evCond(f, cs0, inp0, mat0) ^ evCond(f, cs1, inp1, mat1);
  }
}


CondGarbling gbCond(
    const PRF& f,
    std::span<const Circuit> cs,
    const Label& seed) {
  const auto n = cs.size();
  PRG prg(seed);

  if (n == 1) {
    const auto g = garble(prg, f, cs[0]);
    return {
      std::move(g.material),
      std::move(g.inputEncoding),
      { std::move(g.outputEncoding) },
      { }
    };
  } else {
    auto e = genEncoding(prg, cs[0].nInp);
    const auto s0 = e.zeros[0];
    const auto s1 = e.zeros[0] ^ e.delta;

    // split the vector of circuits
    const std::span<const Circuit> cs0 = cs.subspan(0, n/2);
    const std::span<const Circuit> cs1 = cs.subspan(n/2);

    auto g0 = gbCond(f, cs0, s1);
    auto g1 = gbCond(f, cs1, s0);

    auto branchmat = g0.material ^ g1.material;

    auto g0_ = gbCond_(f, cs0, s0);
    auto g1_ = gbCond_(f, cs1, s1);

    std::span<Label> zeros { e.zeros };
    zeros = zeros.subspan(1);
    auto [mdem, bad0, bad1] = gbDem(prg, f, e.delta, s0, zeros, g0.inputEncoding, g1.inputEncoding);

    g0_.material ^= branchmat;
    g1_.material ^= branchmat;

    std::span mat0 { g1_.material };
    std::span mat1 { g0_.material };
    auto badout0 = evCond(f, cs0, bad0, mat0);
    auto badout1 = evCond(f, cs1, bad1, mat1);

    return {
      std::move(mdem) | std::move(branchmat),
      e,
      std::move(g0.outputEncodings) | std::move(g1.outputEncodings),
      std::vector<Labelling> { std::move(badout0) }
        | std::move(g0.bad)
        | std::vector<Labelling> { std::move(badout1) }
        | std::move(g1.bad)
    };
  }
}



Encoding gb(const PRF& f, const Circuit& c, const Encoding& inputEncoding, std::span<Label>& material) {
  return std::visit(overloaded {
    [&](const Netlist& n) {
      const auto delta = inputEncoding.delta;
      Encoding outputEncoding;
      outputEncoding.delta = delta;
      outputEncoding.zeros.resize(c.nOut);

      NetlistCtxt ctxt;
      ctxt.w = Wiring(c.nWires);
      ctxt.material = material;
      ctxt.inp = std::span<const Label>(inputEncoding.zeros);
      ctxt.out = std::span<Label>(outputEncoding.zeros);
      ctxt.nonce = 0;

      for (const auto& g: n) {
        gbGate(f, g, delta, ctxt);
      }

      material = ctxt.material;

      return outputEncoding;
    },
    [&](const Conditional& cond) {
      Encoding e;

      // TODO gen demux
      // TODO gen mux
      std::vector<Garbling> gs;
      gs.reserve(cond.cs.size());
      for (const auto& c: cond.cs) {
        PRG prg { 0 };
        gs.push_back(garble(prg, f, c)); // TODO properly define seed based on condition
      }
      return e;
    },
    [&](const Sequence& seq) {
      auto encoding = inputEncoding;
      for (const auto& c: seq) {
        encoding = gb(f, c, encoding, material);
      }
      return encoding;
    },
  }, c.content);
}


Garbling garble(PRG& seed, const PRF& f, const Circuit& c) {
  Garbling g;
  g.inputEncoding = genEncoding(seed, c.nInp);
  g.material.resize(c.nRow);
  std::span<Label> mat(g.material);
  g.outputEncoding = gb(f, c, g.inputEncoding, mat);
  return g;
}


Labelling ev(const PRF& f, const Circuit& c, const Labelling& input, std::span<Label>& material) {
  return std::visit(overloaded {
    [&](const Netlist& n) {
      Labelling output(c.nOut);

      NetlistCtxt ctxt;
      ctxt.w = Wiring(c.nWires);
      ctxt.material = material;
      ctxt.inp = std::span<const Label>(input);
      ctxt.out = std::span<Label>(output);
      ctxt.nonce = 0;

      for (const auto& g: n) {
        evGate(f, g, ctxt);
      }

      material = ctxt.material;
      return output;
    },
    [&](const Conditional& cond) {
      evCond(f, std::span { cond.cs }, input, material);
      Labelling out;
      return out;
    },
    [&](const Sequence& seq) {
      auto labelling = input;
      for (const auto& c: seq) {
        labelling = ev(f, c, labelling, material);
      }
      return labelling;
    },
  }, c.content);
}


Encoding genEncoding(PRG& prg, std::size_t n) {
  Encoding out;
  out.delta = prg();
  out.delta[0] = 1;
  out.zeros.resize(n);
  for (auto& z: out.zeros) { z = prg(); }
  return out;
}
