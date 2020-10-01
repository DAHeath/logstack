#include "scheme.h"
#include "gate.h"
#include <thread>

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


Labelling operator^(const Labelling& x, const Labelling& y) {
  const std::size_t n = x.size();
  assert(y.size() == n);

  Labelling out(n);
  for (std::size_t i = 0; i < n; ++i) { out[i] = x[i] ^ y[i]; }
  return out;
}


Encoding operator^(const Encoding& x, const Encoding& y) {
  return { x.zeros ^ y.zeros, x.delta ^ y.delta };
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


void gbGadget_rec(
    std::size_t b,
    const PRF& f,
    const Label& delta,
    std::span<const Label> inseeds,
    std::span<Label> outseeds,
    std::span<const Label> index,
    const Label& isAncestor0,
    std::span<Label>& mat) {
  if (b > 1) {
    std::size_t nonce = 0;
    const auto isAncestorR0 = gbAnd(f, delta, index[0], isAncestor0, mat, nonce);
    const auto isAncestorL0 = isAncestorR0 ^ isAncestor0;

    const auto bL = b/2;
    const auto bR = b - bL;

    mat[0] ^= isAncestorR0 ^ inseeds[0] ^ delta;
    outseeds[0] = inseeds[0] ^ delta;
    mat[1] ^= isAncestorL0 ^ inseeds[2*bL-1] ^ delta;
    outseeds[2*bL-1] = inseeds[2*bL-1] ^ delta;
    mat = mat.subspan(2);

    gbGadget_rec(bL, f, delta, inseeds.subspan(1, 2*bL-1), outseeds.subspan(1, 2*bL-1), index.subspan(1), isAncestorL0, mat);
    gbGadget_rec(bR, f, delta, inseeds.subspan(2*bL), outseeds.subspan(2*bL), index.subspan(1), isAncestorR0, mat);
  }
}


std::vector<Label> gbGadget(
    std::size_t b,
    const PRF& f,
    const Label& delta,
    std::span<const Label> inseeds,
    std::span<const Label> index,
    std::span<Label>& mat) {
  std::vector<Label> out(inseeds.size());

  std::span<Label> outseeds(out);

  const auto isAncestorR0 = index[0];
  const auto isAncestorL0 = index[0] ^ delta;

  const auto bL = b/2;
  const auto bR = b - bL;
  mat[0] ^= isAncestorR0 ^ inseeds[0] ^ delta;
  outseeds[0] = inseeds[0] ^ delta;
  mat[1] ^= isAncestorL0 ^ inseeds[2*bL-1] ^ delta;
  outseeds[2*bL-1] = inseeds[2*bL-1] ^ delta;
  mat = mat.subspan(2);

  gbGadget_rec(bL, f, delta, inseeds.subspan(1, 2*bL-1), outseeds.subspan(1, 2*bL-1), index.subspan(1), isAncestorL0, mat);
  gbGadget_rec(bR, f, delta, inseeds.subspan(2*bL), outseeds.subspan(2*bL), index.subspan(1), isAncestorR0, mat);

  return out;
}


void evGadget_rec(
    std::size_t b,
    const PRF& f,
    std::span<Label> outseeds,
    std::span<const Label> index,
    const Label& isAncestor,
    std::span<Label>& mat) {
  if (b > 1) {
    std::size_t nonce = 0;
    const auto isAncestorR = evAnd(f, index[0], isAncestor, mat, nonce);
    const auto isAncestorL = isAncestorR ^ isAncestor;

    const auto bL = b/2;
    const auto bR = b - bL;

    outseeds[0] = mat[0] ^ isAncestorR;
    outseeds[2*bL-1] = mat[1] ^ isAncestorL;
    mat = mat.subspan(2);

    evGadget_rec(bL, f, outseeds.subspan(1, 2*bL-1), index.subspan(1), isAncestorL, mat);
    evGadget_rec(bR, f, outseeds.subspan(2*bL), index.subspan(1), isAncestorR, mat);
  }
}


std::vector<Label> evGadget(
    std::size_t b,
    const PRF& f,
    std::span<const Label> index,
    std::span<Label>& mat) {
  std::vector<Label> out(2*b - 2);

  std::span<Label> outseeds(out);

  const auto isAncestor1 = index[0];
  const auto isAncestor0 = index[0];

  const auto b0 = b/2;
  const auto b1 = b - b0;

  outseeds[0] = mat[0] ^ isAncestor1;
  outseeds[2*b0-1] = mat[1] ^ isAncestor0;
  mat = mat.subspan(2);

  evGadget_rec(b0, f, outseeds.subspan(1, 2*b0-1), index.subspan(1), isAncestor0, mat);
  evGadget_rec(b1, f, outseeds.subspan(2*b1), index.subspan(1), isAncestor1, mat);

  return out;
}



Material gbCond_(const PRF& f, std::span<const Circuit> cs, const Encoding& e, PRG& prg) {
  const auto b = cs.size();
  if (b == 1) {
    Material mat (cs[0].nRow);
    gb(f, cs[0], e, mat);
    return mat;
  } else {
    const auto n = cs[0].nInp + ilog2(cs.size());

    PRG prg0 = prg.child();
    PRG prg1 = prg.child();

    const auto e0 = genEncoding(prg0, n-1);
    const auto e1 = genEncoding(prg1, n-1);

    const auto b0 = b/2;
    const auto cs0 = cs.subspan(0, b0);
    const auto cs1 = cs.subspan(b0);

    const auto mat0 = gbCond_(f, cs0, e0, prg0);
    const auto mat1 = gbCond_(f, cs1, e1, prg1);

    Material material { mat0.size() + 3*n + 2 };

    std::span<Label> mat { material };

    std::span<const Label> zeros { e.zeros };
    gbDem(prg, f, e.delta, zeros[0], zeros.subspan(1), e0, e1, mat);
    mat.subspan(3*n+2) ^= mat0;
    mat.subspan(3*n+2) ^= mat1;
    return material;
  }
}


Labelling evCond(
    const PRF& f,
    std::span<const Circuit> cs,
    std::span<const Label> seeds,
    const Labelling& inp,
    std::span<Label> mat,
    std::span<Label> muxMat) {

  const auto b = cs.size();
  if (b == 1) {
    return ev(f, cs[0], inp, mat);
  } else {
    const auto n = cs[0].nInp + ilog2(cs.size());
    const auto b0 = b/2;
    const auto cs0 = cs.subspan(0, b0);
    const auto cs1 = cs.subspan(b0);

    const auto s0 = seeds[0];
    const auto seeds0 = seeds.subspan(1, 2*b0-1);
    const auto s1 = seeds[2*b0-1];
    const auto seeds1 = seeds.subspan(2*b0);

    const auto S = inp[0];
    std::span<const Label> inp_(inp);
    inp_ = inp_.subspan(1);
    const auto demOut = evDem(f, S, inp_, mat);
    const auto inp0 = demOut.first;
    const auto inp1 = demOut.second;

    // move past the demux material
    mat = mat.subspan(3*inp_.size() + 2);

    const auto muxMatSize0 = (cs0.size()-1) * (cs0[0].nOut + 2);
    const auto muxMatSize1 = (cs1.size()-1) * (cs1[0].nOut + 2);
    const auto muxMat0 = muxMat.subspan(0, muxMatSize0);
    const auto muxMat1 = muxMat.subspan(muxMatSize0, muxMatSize1);
    const auto muxMatNow = muxMat.subspan(muxMatSize0 + muxMatSize1);

    Material material1 { mat.begin(), mat.end() };
    std::span<Label> mat0 = mat;
    std::span<Label> mat1 { material1 };


    Labelling out0, out1;
    {
      PRG prg(s1);
      mat0 ^= gbCond_(f, cs1, genEncoding(prg, n-1), prg);
      out0 = evCond(f, cs0, seeds0, inp0, mat0, muxMat0);
    }
    {
      PRG prg(s0);
      mat1 ^= gbCond_(f, cs0, genEncoding(prg, n-1), prg);
      out1 = evCond(f, cs1, seeds1, inp1, mat1, muxMat1);
    }

    return evMux(f, S, out0, out1, muxMatNow);
  }
}


CondGarbling gbCond(
    const PRF& f,
    std::span<const Circuit> cs,
    PRG& prg,
    const Encoding& e,
    std::span<const Label> badSeeds,
    std::vector<Labelling>& badInps,
    std::vector<std::span<Label>> badSiblings,
    std::span<Label> muxMat) {
  const auto b = cs.size();

  if (b == 1) {
    Material material (cs[0].nRow);
    const auto d = gb(f, cs[0], e, material);


    std::vector<Labelling> badouts(badInps.size());
    Material scratch = material;
    for (int i = badInps.size()-1; i >= 0; --i) {
      scratch ^= badSiblings[i];
      badouts[i] = ev(f, cs[0], badInps[i], scratch);
    }

    return { material, d, badouts };
  } else {
    const auto n = cs[0].nInp + ilog2(cs.size());

    PRG prg0 = prg.child();
    PRG prg1 = prg.child();

    const auto e0 = genEncoding(prg0, n-1);
    const auto e1 = genEncoding(prg1, n-1);

    Material material(3*(n-1) + 2);
    std::span<const Label> zeros { e.zeros };
    const auto [bad0, bad1] = gbDem(prg, f, e.delta, zeros[0], zeros.subspan(1), e0, e1, material);

    Material scratch = material;
    std::vector<Labelling>& badInps0 = badInps;
    std::vector<Labelling> badInps1(badInps.size());
    for (int i = badInps.size()-1; i >= 0; --i) {
      scratch ^= badSiblings[i].subspan(0, 3*(n-1) + 2);
      badSiblings[i] = badSiblings[i].subspan(3*(n-1) + 2);
      auto demuxed = evDem(f, badInps[i][0], std::span { badInps[i] }.subspan(1), scratch);
      badInps0[i] = demuxed.first;
      badInps1[i] = demuxed.second;
    }
    badInps0.emplace_back(bad0);
    badInps1.emplace_back(bad1);


    const auto b0 = b/2;
    const auto cs0 = cs.subspan(0, b0);
    const auto cs1 = cs.subspan(b0);

    const auto bads0 = badSeeds[0];
    const auto badSeeds0 = badSeeds.subspan(1, 2*b0-1);
    const auto bads1 = badSeeds[2*b0-1];
    const auto badSeeds1 = badSeeds.subspan(2*b0);

    const auto muxMatSize0 = (cs0.size()-1) * (cs0[0].nOut + 2);
    const auto muxMatSize1 = (cs1.size()-1) * (cs1[0].nOut + 2);
    const auto muxMat0 = muxMat.subspan(0, muxMatSize0);
    const auto muxMat1 = muxMat.subspan(muxMatSize0, muxMatSize1);
    const auto muxMatNow = muxMat.subspan(muxMatSize0 + muxMatSize1);


    // immediately garble all of the right branches to compute garbage material
    const auto m1 = gbCond_(f, cs1, e1, prg1);

    // set up bad seed and expand
    PRG prg1_ { bads1 };
    auto m1_ = gbCond_(f, cs1, genEncoding(prg1_, n-1), prg1_);
    m1_ ^= m1;

    // now, recursively garble left branches, using garbage from right
    auto badSiblings0 = badSiblings;
    badSiblings0.emplace_back(m1_);
    const auto [m0, d0, badOuts0] = gbCond(f, cs0, prg0, e0, badSeeds0, badInps0, badSiblings0, muxMat0);

    // now that garbage from left is available,
    // recursively garble right branches
    PRG prg0_ { bads0 };
    auto m0_ = gbCond_(f, cs0, genEncoding(prg0_, n-1), prg0_);
    m0_ ^= m0;
    auto badSiblings1 = badSiblings;
    badSiblings1.emplace_back(m0_);
    const auto [m1_re, d1, badOuts1] = gbCond(f, cs1, prg1, e1, badSeeds1, badInps1, badSiblings1, muxMat1);

    material.resize(material.size() + m0.size());

    std::span<Label> mat = material;
    mat.subspan(3*(n-1) + 2);
    mat ^= m0;
    mat ^= m1;
    const auto d = gbMux(prg, f, e.delta, e.zeros[0], d0, d1, badOuts0.back(), badOuts1.back(), muxMatNow);

    std::vector<Labelling> badOuts(badOuts0.size()-1);

    for (std::size_t i = 0; i < badOuts.size(); ++i) {
      badOuts[i] = evMux(f, badInps[i][0], badOuts0[i], badOuts1[i], muxMatNow);
    }

    return { material, d, badOuts };
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


std::vector<Label> seedTree(const Label& root, std::size_t b) {
  if (b <= 1) {
    return { root };
  } else {
    PRG prg(root);
    const auto l = prg.child();
    const auto r = prg.child();

    const auto bL = b/2;
    const auto bR = b - bL;

    const auto lSeeds = seedTree(l, bL);
    const auto rSeeds = seedTree(r, bR);

    std::vector<Label> out;
    out.push_back(root);
    std::copy(lSeeds.begin(), lSeeds.end(), std::back_inserter(out));
    std::copy(rSeeds.begin(), rSeeds.end(), std::back_inserter(out));
    return out;
  }
}


Encoding gb(const PRF& f, const Circuit& c, const Encoding& inpEnc, std::span<Label> mat) {
  return std::visit(overloaded {
    [&](const Netlist& n) {
      return netlistgb(f, n, inpEnc, mat);
    },
    [&](const Conditional& cond) {
      const auto b = cond.cs.size();
      const auto m = c.nOut;
      PRG prg;
      const auto seed = prg();


      const auto gadgetSize = 3*b;
      const auto muxSize = (b-1)*(m+2);

      auto gadgetMat = mat.subspan(0, gadgetSize);
      const auto bodyMat = mat.subspan(
          gadgetSize,
          mat.size() - muxSize - gadgetSize);
      const auto muxMat = mat.subspan(mat.size() - muxSize);


      std::span<const Label> inp = inpEnc.zeros;
      const auto goodSeeds = seedTree(seed, b);
      const auto badSeeds = gbGadget(b, f, inpEnc.delta, goodSeeds, inp.subspan(log2(b)), gadgetMat);


      PRG seedPRG(seed);
      std::vector<Labelling> badInps;
      const auto [material, d, bo_] = gbCond(f, cond.cs, seedPRG, inpEnc, badSeeds, badInps, { }, muxMat);
      bodyMat ^= material;
      return d;
    },
    [&](const Sequence& seq) {
      auto encoding = inpEnc;
      for (const auto& c: seq) {
        encoding = gb(f, c, encoding, mat);
      }
      return encoding;
    },
  }, c.content);
}


Labelling ev(const PRF& f, const Circuit& c, const Labelling& input, std::span<Label> mat) {
  return std::visit(overloaded {
    [&](const Netlist& n) {
      return netlistev(f, n, input, mat);
    },
    [&](const Conditional& cond) {
      const auto b = cond.cs.size();
      const auto m = c.nOut;

      const auto gadgetSize = 3*b;
      const auto muxSize = (b-1)*(m+2);
      auto gadgetMat = mat.subspan(0, gadgetSize);
      const auto bodyMat = mat.subspan(
          gadgetSize,
          mat.size() - muxSize - gadgetSize);
      const auto muxMat = mat.subspan(mat.size() - muxSize);

      std::span<const Label> inps (input);
      const auto seeds = evGadget(b, f, inps.subspan(0, log2(b)), gadgetMat);


      return evCond(f, cond.cs, seeds, input, bodyMat, muxMat);
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
