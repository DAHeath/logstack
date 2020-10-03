#include "scheme.h"
#include "gate.h"
#include <thread>

#include <iostream>

std::size_t n_netlistgb = 0;


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

    const auto b0 = b/2;
    const auto b1 = b - b0;

    mat[0] ^= isAncestorR0 ^ inseeds[0] ^ delta;
    outseeds[0] = inseeds[0] ^ delta;
    mat[1] ^= isAncestorL0 ^ inseeds[2*b0-1] ^ delta;
    outseeds[2*b0-1] = inseeds[2*b0-1] ^ delta;
    mat = mat.subspan(2);

    gbGadget_rec(b0, f, delta, inseeds.subspan(1, 2*b0-2), outseeds.subspan(1, 2*b0-2), index.subspan(1), isAncestorL0, mat);
    gbGadget_rec(b1, f, delta, inseeds.subspan(2*b0), outseeds.subspan(2*b0), index.subspan(1), isAncestorR0, mat);
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

  const auto b0 = b/2;
  const auto b1 = b - b0;
  mat[0] ^= isAncestorR0 ^ inseeds[0] ^ delta;
  outseeds[0] = inseeds[0] ^ delta;
  mat[1] ^= isAncestorL0 ^ inseeds[2*b0-1] ^ delta;
  outseeds[2*b0-1] = inseeds[2*b0-1] ^ delta;
  mat = mat.subspan(2);

  gbGadget_rec(b0, f, delta, inseeds.subspan(1, 2*b0-2), outseeds.subspan(1, 2*b0-2), index.subspan(1), isAncestorL0, mat);
  gbGadget_rec(b1, f, delta, inseeds.subspan(2*b0), outseeds.subspan(2*b0), index.subspan(1), isAncestorR0, mat);

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

    const auto b0 = b/2;
    const auto b1 = b - b0;

    outseeds[0] = mat[0] ^ isAncestorR;
    outseeds[2*b0-1] = mat[1] ^ isAncestorL;
    mat = mat.subspan(2);

    evGadget_rec(b0, f, outseeds.subspan(1, 2*b0-2), index.subspan(1), isAncestorL, mat);
    evGadget_rec(b1, f, outseeds.subspan(2*b0), index.subspan(1), isAncestorR, mat);
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

  evGadget_rec(b0, f, outseeds.subspan(1, 2*b0-2), index.subspan(1), isAncestor0, mat);
  evGadget_rec(b1, f, outseeds.subspan(2*b0), index.subspan(1), isAncestor1, mat);

  return out;
}



Material gbCond_(const PRF& f, std::span<const Circuit> cs, EncodingView e, PRG& prg, std::size_t toFill) {
  const auto b = cs.size();

  if (e.size() > ilog2(b) + cs[0].nInp) {
    // special case, remove a branch condition non-full tree
    for (int i = ilog2(b)-1; i >= 0; --i) {
      e[i+1] = e[i];
    }
    e = e.subview(1);
  }

  if (b == 1) {
    Material mat (toFill);
    gb(f, cs[0], e, mat);
    for (std::size_t i = cs[0].nRow; i < toFill; ++i) {
      mat[i] = prg();
    }
    return mat;
  } else {
    const auto n = e.size();

    PRG prg0 = prg.child();
    PRG prg1 = prg.child();

    const auto b0 = b/2;

    auto e0 = genEncoding(prg0, n-1);
    auto e1 = genEncoding(prg1, n-1);

    const auto cs0 = cs.subspan(0, b0);
    const auto cs1 = cs.subspan(b0);
    const auto R = std::max(condSize(cs0), condSize(cs1));

    const auto mat0 = gbCond_(f, cs0, e0, prg0, R);
    const auto mat1 = gbCond_(f, cs1, e1, prg1, R);

    Material material(toFill);

    std::span<Label> mat { material };

    std::span<const Label> zeros { e.zeros };
    zeros = zeros.subspan(1);
    gbDem(prg, f, e.delta, e.zeros[0], zeros, e0, e1, mat);

    mat.subspan(3*(n-1)+2) ^= mat0;
    mat.subspan(3*(n-1)+2) ^= mat1;

    for (std::size_t i = R + 3*(n-1) + 2; i < toFill; ++i) {
      material[i] = prg();
    }

    return material;
  }
}


Labelling evCond(
    const PRF& f,
    std::span<const Circuit> cs,
    std::span<const Label> seeds,
    std::span<Label> inp,
    std::span<Label> mat,
    std::span<Label> muxMat) {

  const auto b = cs.size();

  if (inp.size() > ilog2(b) + cs[0].nInp) {
    // special case, remove a branch condition non-full tree
    for (int i = ilog2(b)-1; i >= 0; --i) {
      inp[i+1] = inp[i];
    }
    inp = inp.subspan(1);
  }

  if (b == 1) {
    // special case
    const auto out = ev(f, cs[0], inp, mat);
    return out;
  } else {
    const auto n = inp.size();
    const auto m = cs[0].nOut;
    const auto b0 = b/2;
    const auto b1 = b - b0;
    const auto cs0 = cs.subspan(0, b0);
    const auto cs1 = cs.subspan(b0);

    const auto s0 = seeds[0];
    const auto seeds0 = seeds.subspan(1, 2*b0-2);
    const auto s1 = seeds[2*b0-1];
    const auto seeds1 = seeds.subspan(2*b0);

    const auto S = inp[0];
    const auto demOut = evDem(f, S, inp.subspan(1), mat);
    auto inp0 = demOut.first;
    auto inp1 = demOut.second;

    // move past the demux material
    mat = mat.subspan(3*(n-1) + 2);

    const auto muxMatSize0 = (b0 - 1) * (m + 2);
    const auto muxMatSize1 = (b1 - 1) * (m + 2);

    const auto muxMat0 = muxMat.subspan(0, muxMatSize0);
    const auto muxMat1 = muxMat.subspan(muxMatSize0, muxMatSize1);
    const auto muxMatNow = muxMat.subspan(muxMatSize0 + muxMatSize1);
    assert(muxMatNow.size() == m + 2);

    Material material1 { mat.begin(), mat.end() };
    std::span<Label> mat0 = mat;
    std::span<Label> mat1 { material1 };

    const auto R = std::max(condSize(cs0), condSize(cs1));


    Labelling out0, out1;
    {
    /* std::thread th { [&] { */
      PRG prg(s1);
      auto e1 = genEncoding(prg, n-1);
      mat0 ^= gbCond_(f, cs1, e1, prg, R);
      out0 = evCond(f, cs0, seeds0, inp0, mat0, muxMat0);
    /* }}; */
    }
    {
      PRG prg(s0);
      auto e1 = genEncoding(prg, n-1);
      mat1 ^= gbCond_(f, cs0, e1, prg, R);
      out1 = evCond(f, cs1, seeds1, inp1, mat1, muxMat1);
    }
    /* th.join(); */

    return evMux(f, S, out0, out1, muxMatNow);
  }
}


CondGarbling gbCond(
    const PRF& f,
    std::span<const Circuit> cs,
    PRG& prg,
    EncodingView e,
    std::span<const Label> badSeeds,
    std::vector<Labelling> badInps,
    std::vector<std::span<Label>> badSiblings,
    std::span<Label> muxMat,
    std::size_t toFill) {
  const auto b = cs.size();

  if (e.size() > ilog2(b) + cs[0].nInp) {
    // special case, remove a branch condition non-full tree
    for (int i = ilog2(b)-1; i >= 0; --i) {
      e[i+1] = e[i];
    }
    e = e.subview(1);
    for (auto& bi: badInps) {
      for (int i = ilog2(b)-1; i >= 0; --i) {
        bi[i+1] = bi[i];
      }
      bi.erase(bi.begin());
    }
  }

  if (b == 1) {
    Material material (toFill);
    const auto d = gb(f, cs[0], e, material);

    std::vector<Labelling> badouts(badInps.size());
    Material scratch = material;
    for (int i = badInps.size()-1; i >= 0; --i) {
      scratch ^= badSiblings[i];
      badouts[i] = ev(f, cs[0], badInps[i], scratch);
    }

    for (std::size_t i = cs[0].nRow; i < toFill; ++i) {
      material[i] = prg();
    }

    return { material, d, badouts };
  } else {
    const auto n = e.zeros.size();
    const auto m = cs[0].nOut;

    // derive children of the node seed
    const auto seed0 = prg.child();
    const auto seed1 = prg.child();
    PRG prg0 = seed0;
    PRG prg1 = seed1;

    const auto b0 = b/2;
    const auto b1 = b - b0;

    auto e0 = genEncoding(prg0, n-1);
    auto e1 = genEncoding(prg1, n-1);
    const auto cs0 = cs.subspan(0, b0);
    const auto cs1 = cs.subspan(b0);

    Material material(3*(n-1) + 2);
    std::span<const Label> zeros { e.zeros };
    zeros = zeros.subspan(1);
    const auto [bad0, bad1] = gbDem(prg, f, e.delta, e.zeros[0], zeros, e0, e1, material);

    Material scratch = material;
    std::vector<Labelling> badInps0(badInps.size());
    std::vector<Labelling> badInps1(badInps.size());
    for (int i = badInps.size()-1; i >= 0; --i) {
      // peel off the demux portion of the material of each bad sibling
      // and add it to the scratch material
      scratch ^= badSiblings[i].subspan(0, 3*(n-1) + 2);
      badSiblings[i] = badSiblings[i].subspan(3*(n-1) + 2);
      auto demuxed = evDem(f, badInps[i][0], std::span { badInps[i] }.subspan(1), scratch);
      badInps0[i] = demuxed.first;
      badInps1[i] = demuxed.second;
    }
    badInps0.emplace_back(bad0);
    badInps1.emplace_back(bad1);


    const auto bads0 = badSeeds[0];
    const auto badSeeds0 = badSeeds.subspan(1, 2*b0-2);
    const auto bads1 = badSeeds[2*b0-1];
    const auto badSeeds1 = badSeeds.subspan(2*b0);

    const auto muxMatSize0 = (b0-1) * (m + 2);
    const auto muxMatSize1 = (b1-1) * (m+2);
    const auto muxMat0 = muxMat.subspan(0, muxMatSize0);
    const auto muxMat1 = muxMat.subspan(muxMatSize0, muxMatSize1);
    const auto muxMatNow = muxMat.subspan(muxMatSize0 + muxMatSize1);


    const auto R = std::max(condSize(cs0), condSize(cs1));

    // set up bad seeds and expand
    Material m0_, m1_, m1;
    /* std::thread th { [&] { */
    {
      PRG prg1_ { bads1 };
      auto e1_ = genEncoding(prg1_, n-1);
      m1_ = gbCond_(f, cs1, e1_ , prg1_, R);
      m1 = gbCond_(f, cs1, e1, prg1, R);
      m1_ ^= m1;
    /* }}; */
    }
    {
      PRG prg0_ { bads0 };
      auto e0_ = genEncoding(prg0_, n-1);
      m0_ = gbCond_(f, cs0, e0_, prg0_, R);
    }
    /* th.join(); */

    // immediately garble all of the right branches to compute good material

    // now, recursively garble left branches, using garbage from right
    auto badSiblings0 = badSiblings;
    badSiblings0.push_back(m1_);
    const auto [m0, d0, badOuts0] = gbCond(f, cs0, prg0, e0, badSeeds0, badInps0, badSiblings0, muxMat0, R);

    // now that garbage from left is available,
    // recursively garble right branches
    m0_ ^= m0;
    auto badSiblings1 = badSiblings;
    badSiblings1.push_back(m0_);

    PRG prg1_re = seed1;
    const auto [m1_re, d1, badOuts1] = gbCond(f, cs1, prg1_re, e1, badSeeds1, badInps1, badSiblings1, muxMat1, R);

    material.resize(toFill);

    std::span<Label> mat = material;
    mat = mat.subspan(3*(n-1) + 2);
    mat ^= m0;
    mat ^= m1;
    assert(muxMatNow.size() == m + 2);


    const auto d = gbMux(prg, f, e.delta, e.zeros[0], d0, d1, badOuts0.back(), badOuts1.back(), muxMatNow);

    std::vector<Labelling> badOuts(badOuts0.size()-1);

    for (std::size_t i = 0; i < badOuts.size(); ++i) {
      badOuts[i] = evMux(f, badInps[i][0], badOuts0[i], badOuts1[i], muxMatNow);
    }

    for (std::size_t i = R + 3*(n-1) + 2; i < toFill; ++i) {
      material[i] = prg();
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

    const auto b0 = b/2;
    const auto b1 = b - b0;

    const auto lSeeds = seedTree(l, b0);
    const auto rSeeds = seedTree(r, b1);

    std::vector<Label> out;
    out.push_back(root);
    std::copy(lSeeds.begin(), lSeeds.end(), std::back_inserter(out));
    std::copy(rSeeds.begin(), rSeeds.end(), std::back_inserter(out));
    return out;
  }
}


Encoding gb(const PRF& f, const Circuit& c, EncodingView inpEnc, std::span<Label> mat) {
  return std::visit(overloaded {
    [&](const Netlist& n) {
      ++n_netlistgb;
      return netlistgb(f, n, inpEnc, mat);
    },
    [&](const Conditional& cond) {
      n_netlistgb = 0;
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

      std::span<const Label> gSeeds = goodSeeds;
      // skip past the root seed
      gSeeds = gSeeds.subspan(1);

      const auto badSeeds = gbGadget(
          b, f, inpEnc.delta, gSeeds, inp.subspan(0, ilog2(b)), gadgetMat);




      PRG seedPRG(seed);
      std::vector<Labelling> badInps;

      const auto [material, d, bo_] =
        gbCond(f, cond.cs, seedPRG, inpEnc, badSeeds, badInps, { }, muxMat, condSize(cond.cs));
      bodyMat ^= material;
      return d;
    },
    [&](const Sequence& seq) {
      Encoding e;
      auto encoding = inpEnc;
      for (const auto& c: seq) {
        e = gb(f, c, encoding, mat);
        encoding = e;
      }
      return e;
    },
  }, c.content);
}


Labelling ev(const PRF& f, const Circuit& c, std::span<Label> input, std::span<Label> mat) {
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
      const auto seeds = evGadget(b, f, inps.subspan(0, ilog2(b)), gadgetMat);

      return evCond(f, cond.cs, seeds, input, bodyMat, muxMat);
    },
    [&](const Sequence& seq) {
      auto labelling = input;
      Labelling lab;
      for (const auto& c: seq) {
        lab = ev(f, c, labelling, mat);
        labelling = lab;
        mat = mat.subspan(c.nRow);
      }
      return lab;
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
