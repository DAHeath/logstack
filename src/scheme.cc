#include "scheme.h"

// TODO remove
#include <iostream>

void show(const Label& l) {
  std::uint64_t* xs = (std::uint64_t*)(&l);
  std::cout << std::hex << xs[0] << xs[1] << '\n';
}


template <typename T>
std::span<Label> operator^=(std::span<Label> m0, const T& m1) {
  const auto n = m1.size();
  assert(m0.size() >= n);
  for (std::size_t i = 0; i < n; ++i) { m0[i] ^= m1[i]; }
  return m0;
}


std::pair<Labelling, Labelling> gbDem(
    PRG& seed,
    const PRF& f,
    const Label& delta,
    const Label& S0,
    std::span<const Label> zeros,
    const Encoding& e0,
    const Encoding& e1,
    std::span<Label> mat) {

  const auto n = zeros.size();

  const auto s = S0[0];
  const auto hS0 = f(S0);
  const auto hS1 = f(S0 ^ delta);

  auto dd0 = delta ^ e0.delta;
  auto dd1 = delta ^ e1.delta;
  // some extra random bits to mask the definitely zero lsb of delta ^ delta0/delta1
  const auto mask = seed();
  dd0[0] = mask[0];
  dd1[0] = mask[1];

  mat[0] ^= hS0 ^ dd0;
  mat[1] ^= hS1 ^ dd1;

  auto bad_diff0 = hS1 ^ mat[0];
  auto bad_diff1 = hS0 ^ mat[1];
  // we know that the lsb of delta ^ delta0/delta1 should be 0
  bad_diff0[0] = 0;
  bad_diff1[0] = 0;

  mat = mat.subspan(2);

  const auto gbDem1 = [&](const Label& A0) -> std::pair<Label, Label> {
    const auto a = A0[0];

    const auto hA0 = f(A0);
    const auto hA1 = f(A0 ^ delta);

    const auto row0 = hA0 ^ ((a != s) ? A0 : 0);
    const auto row1 = hA1 ^ ((a != s) ? delta : A0);

    const auto good = a ? row1 : row0;
    const auto bad = good ^ A0;
    mat[0] ^= row0 ^ row1;
    return { good, bad };
  };

  Labelling bad0(n);
  Labelling bad1(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto [g, b] = gbDem1(zeros[i]);
    mat = mat.subspan(1);
    const auto good0 = g ^ (g[0] ? (delta ^ e0.delta) : 0);
    const auto good1 = g ^ (g[0] ? (delta ^ e1.delta) : 0);
    bad0[i] = b ^ (b[0] ? bad_diff0 : 0);
    bad1[i] = b ^ (b[0] ? bad_diff1 : 0);
    mat[0] ^= hS0 ^ good0 ^ e0.zeros[i];
    mat[1] ^= hS1 ^ good1 ^ e1.zeros[i];
    bad0[i] ^= hS1 ^ mat[0];
    bad1[i] ^= hS0 ^ mat[1];
    mat = mat.subspan(2);
  }

  return { bad0, bad1 };
}


std::pair<Labelling, Labelling> evDem(
    const PRF& f,
    const Label& S,
    std::span<const Label> inp,
    std::span<Label> mat) {

  const auto n = inp.size();

  const auto s = S[0];
  const auto hS = f(S);
  auto diff0 = hS ^ mat[0];
  auto diff1 = hS ^ mat[1];
  diff0[0] = 0;
  diff1[0] = 0;

  const auto& evDem1 = [&](const Label& A) -> std::pair<Label, Label> {
    const auto a = A[0];
    const auto hA = f(A);
    const auto r0 = hA ^ (a ? mat[0] : 0) ^ ((a != s) ? A : 0);
    const auto r1 = hA ^ (a ? mat[0] : 0) ^ ((a != s) ? 0 : A);
    return { r0, r1 };
  };

  mat = mat.subspan(2);

  Labelling inp0(n);
  Labelling inp1(n);
  for (std::size_t i = 0; i < n; ++i) {
    auto [i0, i1] = evDem1(inp[i]);
    inp0[i] = i0 ^ (i0[0] ? diff0 : 0);
    inp1[i] = i1 ^ (i1[0] ? diff1 : 0);
    mat = mat.subspan(1);

    inp0[i] ^= hS ^ mat[0];
    inp1[i] ^= hS ^ mat[1];
    mat = mat.subspan(2);
  }

  return { inp0, inp1 };
}


Encoding gbMux(
    PRG& prg,
    const PRF& f,
    const Label& delta,
    const Label& S0,
    const Encoding& good0,
    const Encoding& good1,
    const Labelling& bad0,
    const Labelling& bad1,
    std::span<Label> mat) {

  const auto n = good0.zeros.size();

  const auto s = S0[0];
  const auto hS0 = f(S0);
  const auto hS1 = f(S0 ^ delta);

  auto dd0 = delta ^ good0.delta;
  auto dd1 = delta ^ good1.delta;
  // some extra random bits to mask the definitely zero lsb of delta ^ delta0/delta1
  const auto mask = prg();
  dd0[0] = mask[0];
  dd1[0] = mask[1];

  mat[0] ^= hS0 ^ dd0;
  mat[1] ^= hS1 ^ dd1;

  auto bad_diff0 = hS1 ^ mat[0];
  auto bad_diff1 = hS0 ^ mat[1];
  bad_diff0[0] = 0;
  bad_diff1[0] = 0;

  mat = mat.subspan(2);

  Encoding e;
  e.delta = delta;
  e.zeros.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    auto b0 = bad0[i];
    auto b1 = bad1[i];
    auto g0 = good0.zeros[i];
    auto g1 = good1.zeros[i];
    b0 ^= b0[0] ? bad_diff0 : 0;
    b1 ^= b1[0] ? bad_diff1 : 0;
    g0 ^= g0[0] ? (delta ^ good0.delta) : 0;
    g1 ^= g1[0] ? (delta ^ good1.delta) : 0;

    // implement the following garbled truth table for arbitrary X:
    // S0 -> Good0 ^ Bad1 ^ X
    // S1 -> Good1 ^ Bad0 ^ X
    mat[0] = hS0 ^ hS1 ^ g0 ^ g1 ^ b0 ^ b1;

    e.zeros[i] = s ? (hS1 ^ g1 ^ b0) : (hS0 ^ g0 ^ b1);
    mat = mat.subspan(1);
  }

  return e;
}


Labelling evMux(
    const PRF& f,
    const Label& S,
    const Labelling& Xs,
    const Labelling& Ys,
    std::span<Label> mat) {

  const auto n = Xs.size();

  const auto s = S[0];
  const auto hS = f(S);
  auto diff0 = hS ^ mat[0];
  auto diff1 = hS ^ mat[1];
  diff0[0] = 0;
  diff1[0] = 0;

  mat = mat.subspan(2);

  Labelling out(n);
  for (std::size_t i = 0; i < n; ++i) {
    auto X = Xs[i];
    auto Y = Ys[i];
    X ^= X[0] ? diff0 : 0;
    Y ^= Y[0] ? diff1 : 0;
    out[i] = X ^ Y ^ hS ^ (s ? mat[0] : 0);
    mat = mat.subspan(1);
  }
  return out;
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

      ctxt.material[0] ^= hA0 ^ hA1 ^ B0D;
      ctxt.material[1] ^= hB0 ^ hB1 ^ A0D;

      ctxt.nonce += 2;
      ctxt.material = ctxt.material.subspan(2);

      C0 = A0&B0 ^ X ^ Y;
      break;
    }

    case GateType::XOR:
      C0 = A0 ^ B0;
      break;

    case GateType::NOT:
      C0 = A0 ^ delta;
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

    case GateType::NOT:
      C = A;
      break;
  }
}


Encoding gbCond_(const PRF& f, std::span<const Circuit> cs, const Label& seed, std::span<Label> mat) {
  const auto n = cs.size();
  PRG prg(seed);

  if (n == 1) {
    return garble(prg, f, cs[0], mat).inputEncoding;
  } else {
    // Generate a fresh encoding.
    // All branches should have same number of inputs.
    const auto n = cs[0].nInp + ilog2(cs.size());
    auto e = genEncoding(prg, n);
    const auto S0 = e.zeros[0];
    const auto S1 = e.zeros[0] ^ e.delta;

    // split the vector of circuits
    const std::span<const Circuit> cs0 = cs.subspan(0, n/2);
    const std::span<const Circuit> cs1 = cs.subspan(n/2);

    // skip past the demux material
    const auto demMat = mat;
    mat = mat.subspan(3*(n-1) + 2);

    const auto e0 = gbCond_(f, cs0, S1, mat);
    const auto e1 = gbCond_(f, cs1, S0, mat);

    std::span<Label> zeros { e.zeros };
    zeros = zeros.subspan(1);
    gbDem(prg, f, e.delta, S0, zeros, e0, e1, demMat);

    return e;
  }
}


Labelling evCond(
    const PRF& f,
    std::span<const Circuit> cs,
    const Labelling& inp,
    std::span<Label> mat,
    std::span<Label> muxMat) {

  const auto n = cs.size();

  if (n == 1) {
    return ev(f, cs[0], inp, mat);
  } else {
    // split input into the branch condition and the remaining wires
    const auto S = inp[0];
    std::span<const Label> inp_(inp);
    inp_ = inp_.subspan(1);

    // evaluate the demux to compute inputs for both branches
    const auto [inp0, inp1] = evDem(f, S, inp_, mat);
    // move past the demux material
    mat = mat.subspan(3*inp_.size() + 2);

    // split the vector of circuits
    const std::span<const Circuit> cs0 = cs.subspan(0, n/2);
    const std::span<const Circuit> cs1 = cs.subspan(n/2);
    const auto muxMatSize0 = (cs0.size()-1) * (cs0[0].nOut + 2);
    const auto muxMatSize1 = (cs1.size()-1) * (cs1[0].nOut + 2);
    const auto muxMat0 = muxMat.subspan(0, muxMatSize0);
    const auto muxMat1 = muxMat.subspan(muxMatSize0, muxMatSize1);
    const auto muxMatNow = muxMat.subspan(muxMatSize0 + muxMatSize1);

    // we must make one copy of the material so we can symmetrically try second branch
    const auto mat0 = mat;
    Material mat1 { mat.begin(), mat.end() };

    // garble straight into material so as to unstack
    gbCond_(f, cs1, S, mat0);
    const auto l0 = evCond(f, cs0, inp0, mat0, muxMat0);

    gbCond_(f, cs0, S, mat1);
    const auto l1 = evCond(f, cs1, inp1, mat1, muxMat1);

    return evMux(f, S, l0, l1, muxMatNow);
  }
}


// Garble a vector of conditionally composed circuits, starting from a seed,
// into the specified material.
//
// While the material for most of the circuit is stacked with XOR,
// the multiplexer that collects garbage is not.
// Thus, we reference the main material and the mux material by different spans.
Interface gbCond(
    const PRF& f,
    std::span<const Circuit> cs,
    const Label& seed,
    std::span<Label> mat,
    std::span<Label> muxMat) {
  const auto n = cs.size();
  PRG prg(seed);

  if (n == 1) {
    return garble(prg, f, cs[0], mat);
  } else {
    const auto n = cs[0].nInp + ilog2(cs.size());
    auto e = genEncoding(prg, n);
    const auto S0 = e.zeros[0];
    const auto S1 = e.zeros[0] ^ e.delta;

    // split the vector of circuits
    const std::span<const Circuit> cs0 = cs.subspan(0, n/2);
    const std::span<const Circuit> cs1 = cs.subspan(n/2);

    // skip the demux material for now
    const auto demSize = 3*(n-1) + 2;
    const auto demMat = mat.subspan(0, demSize);
    const auto branchmat = mat.subspan(demSize);

    const auto muxMatSize0 = (cs0.size()-1) * (cs0[0].nOut + 2);
    const auto muxMatSize1 = (cs1.size()-1) * (cs1[0].nOut + 2);
    const auto muxMat0 = muxMat.subspan(0, muxMatSize0);
    const auto muxMat1 = muxMat.subspan(muxMatSize0, muxMatSize1);
    const auto muxMatNow = muxMat.subspan(muxMatSize0 + muxMatSize1);

    Interface i0, i1;
    // because garbling the conditional recursively involves unstacking and evaluating,
    // we cannot garble the material in place, and instead must garble into a fresh buffer.
    {
      Material mat0(branchmat.size());
      i0 = gbCond(f, cs0, S1, mat0, muxMat0);
      branchmat ^= mat0;
    }
    {
      Material mat1(branchmat.size());
      i1 = gbCond(f, cs1, S0, mat1, muxMat1);
      branchmat ^= mat1;
    }

    // now, with the input encodings available, we can garble the demux
    // into the front of the material
    std::span<Label> zeros { e.zeros };
    zeros = zeros.subspan(1);
    auto [bad0, bad1] = gbDem(prg, f, e.delta, S0, zeros, i0.inputEncoding, i1.inputEncoding, demMat);

    Encoding e0_, e1_;
    // copy the stacked material so as not to trash it
    // then garble a branch incorrectly to emulate bad evaluation of the other
    // branch (note garbling in place is safe because gbCond_ does not
    // recursively involve evaluation).
    {
      Material mat0(branchmat.begin(), branchmat.end());
      e0_ = gbCond_(f, cs1, S1, mat0);
      bad0 = evCond(f, cs0, bad0, mat0, muxMat0);
    }
    {
      Material mat1(branchmat.begin(), branchmat.end());
      e1_ = gbCond_(f, cs0, S0, mat1);
      bad1 = evCond(f, cs1, bad1, mat1, muxMat1);
    }


    const auto eout = gbMux(
        prg, f, e.delta, S0,
        i0.outputEncoding,
        i1.outputEncoding,
        bad0,
        bad1,
        muxMatNow);

    return { e, eout };
  }
}


void gbTrans(PRG& seed, const Encoding& inp, const Encoding& out, std::span<Label> mat) {
  const auto n = inp.zeros.size();
  assert(out.zeros.size() == n);

  const auto dd = inp.delta ^ out.delta;
  const auto mask = seed();
  mat[0] = dd;
  mat[0][0] = mask[0];
  mat = mat.subspan(1);

  for (std::size_t i = 0; i < n; ++i) {
    const auto in = inp.zeros[i] ^ (inp.zeros[i][0] ? dd : 0);
    mat[0] = in ^ out.zeros[i];
    mat = mat.subspan(1);
  }
}


Labelling evTrans(const Labelling& inp, std::span<Label> mat) {
  auto dd = mat[0];
  dd[0] = 0;
  mat = mat.subspan(1);

  const auto n = inp.size();
  Labelling out(n);
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = inp[i] ^ mat[0] ^ (inp[i][0] ? dd : 0);
    mat = mat.subspan(1);
  }
  return out;
}


Encoding gb(const PRF& f, const Circuit& c, const Encoding& inputEncoding, std::span<Label> mat) {
  return std::visit(overloaded {
    [&](const Netlist& n) {
      const auto delta = inputEncoding.delta;
      Encoding outputEncoding;
      outputEncoding.delta = delta;
      outputEncoding.zeros.resize(c.nOut);

      NetlistCtxt ctxt;
      ctxt.w = Wiring(n.size() - c.nOut + 2);
      ctxt.material = mat;
      ctxt.inp = std::span<const Label>(inputEncoding.zeros);
      ctxt.out = std::span<Label>(outputEncoding.zeros);
      ctxt.nonce = 0;

      for (const auto& g: n) {
        gbGate(f, g, delta, ctxt);
      }

      return outputEncoding;
    },
    [&](const Conditional& cond) {
      PRG prg;
      const auto seed = prg();

      const auto b = cond.cs.size();
      const auto n = c.nInp;
      const auto m = c.nOut;

      const auto transSize = n+1;
      const auto transMat = mat.subspan(0, transSize);
      const auto muxSize = (b-1)*(m+2);
      const auto bodyMat = mat.subspan(
          transSize,
          mat.size() - transSize - muxSize);
      const auto muxMat = mat.subspan(mat.size() - muxSize);
      const auto interface = gbCond(f, cond.cs, seed, bodyMat, muxMat);

      gbTrans(prg, inputEncoding, interface.inputEncoding, transMat);

      return interface.outputEncoding;
    },
    [&](const Sequence& seq) {
      auto encoding = inputEncoding;
      for (const auto& c: seq) {
        encoding = gb(f, c, encoding, mat);
      }
      return encoding;
    },
  }, c.content);
}


Interface garble(PRG& seed, const PRF& f, const Circuit& c, std::span<Label> mat) {
  Interface g;
  g.inputEncoding = genEncoding(seed, c.nInp);
  g.outputEncoding = gb(f, c, g.inputEncoding, mat);
  return g;
}


Labelling ev(const PRF& f, const Circuit& c, const Labelling& input, std::span<Label> mat) {
  return std::visit(overloaded {
    [&](const Netlist& n) {
      Labelling output(c.nOut);

      NetlistCtxt ctxt;
      ctxt.w = Wiring(n.size() - c.nOut + 2);
      ctxt.material = mat;
      ctxt.inp = std::span<const Label>(input);
      ctxt.out = std::span<Label>(output);
      ctxt.nonce = 0;

      for (const auto& g: n) {
        evGate(f, g, ctxt);
      }

      return output;
    },
    [&](const Conditional& cond) {
      const auto b = cond.cs.size();
      const auto n = c.nInp;
      const auto m = c.nOut;

      const auto transSize = n+1;
      const auto transMat = mat.subspan(0, transSize);
      const auto muxSize = (b-1)*(m+2);
      const auto bodyMat = mat.subspan(
          transSize,
          mat.size() - transSize - muxSize);
      const auto muxMat = mat.subspan(mat.size() - muxSize);

      const auto translated = evTrans(input, transMat);
      return evCond(f, std::span { cond.cs }, translated, bodyMat, muxMat);
    },
    [&](const Sequence& seq) {
      auto labelling = input;
      for (const auto& c: seq) {
        labelling = ev(f, c, labelling, mat);
        mat = mat.subspan(c.nRow);
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
