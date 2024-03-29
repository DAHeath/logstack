#include "scheme.h"
#include "gate.h"
#include <thread>

#include <iostream>

std::size_t n_netlistgb = 0;


std::atomic<std::size_t> available_threads = 1;


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

HalfGarbling operator^(const HalfGarbling& x, const HalfGarbling& y) {
  return { x.material ^ y.material, x.outEnc ^ y.outEnc };
}

HalfGarbling operator^(const HalfGarbling& x, const HalfGarblingView& y) {
  HalfGarbling out;
  out.material = x.material;
  out.material ^= y.material;
  out.outEnc = x.outEnc ^ y.outEnc;
  return out;
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



HalfGarbling gbCond_(const PRF& f, std::span<const Circuit> cs, EncodingView e, PRG& prg, std::size_t toFill) {
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
    const auto d = gb(f, cs[0], e, mat);
    for (std::size_t i = cs[0].nRow; i < toFill; ++i) {
      mat[i] = prg();
    }
    return { mat, d };
  } else {
    const auto n = e.size();
    const auto demSize = 3*(n-1)+2;

    PRG prg0 = prg.child();
    PRG prg1 = prg.child();

    const auto b0 = b/2;

    const auto cs0 = cs.subspan(0, b0);
    const auto cs1 = cs.subspan(b0);
    const auto R = toFill - demSize;

    auto e0 = genEncoding(prg0, n-1);
    auto e1 = genEncoding(prg1, n-1);


    Material material(toFill);
    std::span<Label> mat { material };

    std::span<const Label> zeros { e.zeros };
    zeros = zeros.subspan(1);
    gbDem(prg, f, e.delta, e.zeros[0], zeros, e0, e1, mat);


    Material mat0, mat1;
    Encoding d0, d1;

    if (available_threads > 0) {
      --available_threads;
      std::thread th { [&] {
        const auto gb0 = gbCond_(f, cs0, e0, prg0, R);
        mat0 = gb0.material;
        d0 = gb0.outEnc;
      }};
      const auto gb1 = gbCond_(f, cs1, e1, prg1, R);
      mat1 = gb1.material;
      d1 = gb1.outEnc;
      th.join();
      ++available_threads;
    } else {
      const auto gb0 = gbCond_(f, cs0, e0, prg0, R);
      const auto gb1 = gbCond_(f, cs1, e1, prg1, R);
      mat0 = gb0.material;
      mat1 = gb1.material;
      d0 = gb0.outEnc;
      d1 = gb1.outEnc;
    }
    mat.subspan(demSize) ^= mat0;
    mat.subspan(demSize) ^= mat1;
    return { material, d0 ^ d1 };

  }
}


std::tuple<std::span<Label>, std::span<Label>, std::span<Label>>
parseMux(std::size_t b, std::size_t m, std::span<Label> muxMat) {
  const auto b0 = b / 2;
  const auto b1 = b - b0;
  const auto muxMatSize0 = (b0 - 1) * (m + 2);
  const auto muxMatSize1 = (b1 - 1) * (m + 2);

  const auto muxMat0 = muxMat.subspan(0, muxMatSize0);
  const auto muxMat1 = muxMat.subspan(muxMatSize0, muxMatSize1);
  const auto muxMatNow = muxMat.subspan(muxMatSize0 + muxMatSize1);
  assert(muxMatNow.size() == m + 2);

  return { muxMat0, muxMat1, muxMatNow };
}


std::tuple<Label&, std::span<Label>, Label&, std::span<Label>>
parseSeeds(std::size_t b, std::span<Label> seeds) {
  const auto b0 = b/2;
  return { seeds[0], seeds.subspan(1, 2*b0-2), seeds[2*b0-1], seeds.subspan(2*b0) };
}


template <typename T>
std::pair<std::span<T>, std::span<T>> halves(std::span<T> xs) {
  return { xs.subspan(0, xs.size()/2), xs.subspan(xs.size()/2) };
}



Labelling evCond(
    const PRF& f,
    std::span<const Circuit> cs,
    std::span<Label> seeds,
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
    const auto [cs0_, cs1_] = halves(cs);
    const auto cs0 = cs0_;
    const auto cs1 = cs1_;

    const auto [s0_, seeds0_, s1_, seeds1_] = parseSeeds(b, seeds);
    const auto seeds0 = seeds0_;
    const auto seeds1 = seeds1_;
    const auto s0 = s0_;
    const auto s1 = s1_;

    const auto S = inp[0];
    const auto demOut = evDem(f, S, inp.subspan(1), mat);
    auto inp0 = demOut.first;
    auto inp1 = demOut.second;

    // move past the demux material
    mat = mat.subspan(3*(n-1) + 2);

    const auto [muxMat0_, muxMat1_, muxMatNow] = parseMux(b, m, muxMat);
    auto muxMat0 = muxMat0_;
    auto muxMat1 = muxMat1_;

    Material material1 { mat.begin(), mat.end() };
    std::span<Label> mat0 = mat;
    std::span<Label> mat1 { material1 };

    const auto R = std::max(condSize(cs0), condSize(cs1));


    Labelling out0, out1;
    /* { */

    const auto F = [&] {
      PRG prg(s1);
      auto e1 = genEncoding(prg, n-1);
      const auto [m1, d1] = gbCond_(f, cs1, e1, prg, R);
      mat0 ^= m1;
      out0 = evCond(f, cs0, seeds0, inp0, mat0, muxMat0);
    };
    const auto G = [&] {
      PRG prg(s0);
      auto e0 = genEncoding(prg, n-1);
      const auto [m0, d0] = gbCond_(f, cs0, e0, prg, R);
      mat1 ^= m0;
      out1 = evCond(f, cs1, seeds1, inp1, mat1, muxMat1);
    };

    if (available_threads > 0) {
      --available_threads;
      std::thread th { F };
      G();
      th.join();
      ++available_threads;
    } else {
      F();
      G();
    }

    return evMux(f, S, out0, out1, muxMatNow);
  }
}


CondGarbling gbCond(
    const PRF& f,
    std::span<const Circuit> cs,
    PRG& prg,
    EncodingView e,
    HalfGarblingView garbling,
    std::span<Label> badSeeds,
    std::vector<std::span<Label>> badInps,
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
      bi = bi.subspan(1);
    }
  }

  if (b == 1) {
    std::vector<Labelling> badouts(badInps.size());
    Material scratch { garbling.material.begin(), garbling.material.end() };
    for (int i = badInps.size()-1; i >= 0; --i) {
      scratch ^= badSiblings[i];
      badouts[i] = ev(f, cs[0], badInps[i], scratch);
    }

    return { garbling.outEnc, badouts };
  } else {
    const auto n = e.zeros.size();
    const auto m = cs[0].nOut;
    const auto demSize = 3*(n-1)+2;

    // derive children of the node seed
    const auto seed0 = prg.child();
    const auto seed1 = prg.child();
    PRG prg0 = seed0;
    PRG prg1 = seed1;

    const auto b0 = b/2;

    auto e0 = genEncoding(prg0, n-1);
    auto e1 = genEncoding(prg1, n-1);
    const auto cs0 = cs.subspan(0, b0);
    const auto cs1 = cs.subspan(b0);

    Material material(demSize);
    std::span<const Label> zeros { e.zeros };
    zeros = zeros.subspan(1);
    auto [bad0, bad1] = gbDem(prg, f, e.delta, e.zeros[0], zeros, e0, e1, material);

    Material scratch = material;
    std::vector<Labelling> badInpsStore0(badInps.size());
    std::vector<Labelling> badInpsStore1(badInps.size());
    std::vector<std::span<Label>> badInps0(badInps.size());
    std::vector<std::span<Label>> badInps1(badInps.size());
    for (int i = badInps.size()-1; i >= 0; --i) {
      // peel off the demux portion of the material of each bad sibling
      // and add it to the scratch material
      scratch ^= badSiblings[i].subspan(0, demSize);
      badSiblings[i] = badSiblings[i].subspan(demSize);
      auto demuxed = evDem(f, badInps[i][0], std::span { badInps[i] }.subspan(1), scratch);
      badInpsStore0[i] = demuxed.first;
      badInpsStore1[i] = demuxed.second;
      badInps0[i] = badInpsStore0[i];
      badInps1[i] = badInpsStore1[i];
    }
    badInps0.emplace_back(bad0);
    badInps1.emplace_back(bad1);

    const auto [bads0_, badSeeds0_, bads1_, badSeeds1_] = parseSeeds(b, badSeeds);
    auto badSeeds0 = badSeeds0_;
    auto badSeeds1 = badSeeds1_;
    auto bads0 = bads0_;
    auto bads1 = bads1_;

    const auto [muxMat0_, muxMat1_, muxMatNow] = parseMux(b, m, muxMat);
    const auto muxMat0 = muxMat0_;
    const auto muxMat1 = muxMat1_;

    const auto R = toFill - demSize;

    garbling.material = garbling.material.subspan(demSize, R);
    auto garbling0 = gbCond_(f, cs0, e0, prg0, R);
    auto garbling1 = garbling0 ^ garbling;

    // set up bad seeds and expand
    Material m0, m1;
    Encoding d0, d1;
    std::vector<Labelling> badOuts0, badOuts1;


    const auto F = [&] {
      PRG prg1_ { bads1 };
      auto e1_ = genEncoding(prg1_, n-1);
      auto [m1_, d1_] = gbCond_(f, cs1, e1_ , prg1_, R);
      m1_ ^= garbling1.material;
      auto badSiblings0 = badSiblings;
      badSiblings0.push_back(m1_);
      PRG prg0_re = seed0;
      auto gb0 = gbCond(f, cs0, prg0_re, e0, garbling0, badSeeds0, badInps0, badSiblings0, muxMat0, R);
      badOuts0 = gb0.badOuts;
      d0 = gb0.outEnc;
    };
    const auto G = [&] {
      PRG prg0_ { bads0 };
      auto e0_ = genEncoding(prg0_, n-1);
      auto [m0_, d0_] = gbCond_(f, cs0, e0_, prg0_, R);

      m0_ ^= garbling0.material;
      auto badSiblings1 = badSiblings;
      badSiblings1.push_back(m0_);

      PRG prg1_re = seed1;
      auto gb1 = gbCond(f, cs1, prg1_re, e1, garbling1, badSeeds1, badInps1, badSiblings1, muxMat1, R);
      badOuts1 = gb1.badOuts;
      d1 = gb1.outEnc;
    };

    if (available_threads > 0) {
      --available_threads;
      std::thread th { F };
      G();
      th.join();
      ++available_threads;
    } else {
      F();
      G();
    }

    // immediately garble all of the right branches to compute good material

    // now, recursively garble left branches, using garbage from right

    // now that garbage from left is available,
    // recursively garble right branches

    const auto d = gbMux(
        prg, f, e.delta, e.zeros[0], d0, d1, badOuts0.back(), badOuts1.back(), muxMatNow);

    std::vector<Labelling> badOuts(badOuts0.size()-1);
    for (std::size_t i = 0; i < badOuts.size(); ++i) {
      badOuts[i] = evMux(f, badInps[i][0], badOuts0[i], badOuts1[i], muxMatNow);
    }

    return { d, badOuts };
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


      const auto gadgetSize = 4*b;
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

      auto badSeeds = gbGadget(
          b, f, inpEnc.delta, gSeeds, inp.subspan(0, ilog2(b)), gadgetMat);




      PRG seedPRG(seed);

      const auto [mat, dmid] = gbCond_(f, cond.cs, inpEnc, seedPRG, condSize(cond.cs));
      bodyMat ^= mat;


      std::vector<std::span<Label>> badInps(0);

      PRG seedPRG2(seed);
      HalfGarblingView view { bodyMat, dmid };
      const auto [d, bo_] =
        gbCond(f, cond.cs, seedPRG2,
            inpEnc, view, badSeeds, badInps, { }, muxMat, condSize(cond.cs));
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

      const auto gadgetSize = 4*b;
      const auto muxSize = (b-1)*(m+2);
      auto gadgetMat = mat.subspan(0, gadgetSize);
      const auto bodyMat = mat.subspan(
          gadgetSize,
          mat.size() - muxSize - gadgetSize);
      const auto muxMat = mat.subspan(mat.size() - muxSize);

      std::span<const Label> inps (input);
      auto seeds = evGadget(b, f, inps.subspan(0, ilog2(b)), gadgetMat);

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
