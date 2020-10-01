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


/* Encoding gbMux( */
/*     PRG& prg, */
/*     const PRF& f, */
/*     const Label& delta, */
/*     const Label& S0, */
/*     const Encoding& good0, */
/*     const Encoding& good1, */
/*     const Labelling& bad0, */
/*     const Labelling& bad1, */
/*     std::span<Label> mat) { */

/*   const auto n = good0.zeros.size(); */

/*   const auto s = S0[0]; */
/*   const auto hS0 = f(S0); */
/*   const auto hS1 = f(S0 ^ delta); */

/*   auto dd0 = delta ^ good0.delta; */
/*   auto dd1 = delta ^ good1.delta; */
/*   // some extra random bits to mask the definitely zero lsb of delta ^ delta0/delta1 */
/*   const auto mask = prg(); */
/*   dd0[0] = mask[0]; */
/*   dd1[0] = mask[1]; */

/*   mat[0] ^= hS0 ^ dd0; */
/*   mat[1] ^= hS1 ^ dd1; */

/*   auto bad_diff0 = hS1 ^ mat[0]; */
/*   auto bad_diff1 = hS0 ^ mat[1]; */
/*   bad_diff0[0] = 0; */
/*   bad_diff1[0] = 0; */

/*   mat = mat.subspan(2); */

/*   Encoding e; */
/*   e.delta = delta; */
/*   e.zeros.resize(n); */
/*   for (std::size_t i = 0; i < n; ++i) { */
/*     auto b0 = bad0[i]; */
/*     auto b1 = bad1[i]; */
/*     auto g0 = good0.zeros[i]; */
/*     auto g1 = good1.zeros[i]; */
/*     b0 ^= b0[0] ? bad_diff0 : 0; */
/*     b1 ^= b1[0] ? bad_diff1 : 0; */
/*     g0 ^= g0[0] ? (delta ^ good0.delta) : 0; */
/*     g1 ^= g1[0] ? (delta ^ good1.delta) : 0; */

/*     // implement the following garbled truth table for arbitrary X: */
/*     // S0 -> Good0 ^ Bad1 ^ X */
/*     // S1 -> Good1 ^ Bad0 ^ X */
/*     mat[0] = hS0 ^ hS1 ^ g0 ^ g1 ^ b0 ^ b1; */

/*     e.zeros[i] = s ? (hS1 ^ g1 ^ b0) : (hS0 ^ g0 ^ b1); */
/*     mat = mat.subspan(1); */
/*   } */

/*   return e; */
/* } */


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

  const auto isAncestorR = index[0];
  const auto isAncestorL = index[0];

  const auto bL = b/2;
  const auto bR = b - bL;

  outseeds[0] = mat[0] ^ isAncestorR;
  outseeds[2*bL-1] = mat[1] ^ isAncestorL;
  mat = mat.subspan(2);

  evGadget_rec(bL, f, outseeds.subspan(1, 2*bL-1), index.subspan(1), isAncestorL, mat);
  evGadget_rec(bR, f, outseeds.subspan(2*bL), index.subspan(1), isAncestorR, mat);

  return out;
}


Interface gbBranches(
    const PRF& f,
    std::span<const Circuit> cs,
    const Label& seed,
    std::span<Label> mat) {
  const auto b = cs.size();
  PRG prg(seed);

  if (b == 1) {
    return garble(prg, f, cs[0], mat);
  } else {
    const auto seedl = prg();
    const auto seedr = prg();

    // skip past the demux material
    const auto demMat = mat;
    // All branches should have same number of inputs.
    const auto n = cs[0].nInp + ilog2(cs.size());
    mat = mat.subspan(3*(n-1) + 2);

    const auto i0 = gbBranches(f, cs.subspan(0, b/2), seedl, mat);
    const auto i1 = gbBranches(f, cs.subspan(b/2), seedr, mat);

    // Generate a fresh encoding.
    auto e = genEncoding(prg, n);
    std::span<Label> zeros { e.zeros };
    gbDem(prg, f, e.delta, zeros[0], zeros.subspan(1), i0.inpEnc, i1.inpEnc, demMat);

    return { e, i0.outEnc ^ i1.outEnc };
  }
}


Labelling evBranches(
    const PRF& f,
    std::span<const Circuit> cs,
    std::span<const Label> seeds,
    const Labelling& inp,
    std::span<Label> mat,
    std::span<Label> muxMat) {
  const auto b = cs.size();

  if (b == 1) {
    return { ev(f, cs[0], inp, mat) };
  } else {
    const auto bL = b/2;

    // split input into the branch condition and the remaining wires
    const auto S = inp[0];
    std::span<const Label> inp_(inp);
    inp_ = inp_.subspan(1);

    // evaluate the demux to compute inputs for both branches
    const auto demOut = evDem(f, S, inp_, mat);
    const auto inp0 = demOut.first;
    const auto inp1 = demOut.second;
    // move past the demux material
    mat = mat.subspan(3*inp_.size() + 2);

    // we must make one copy of the material so we can symmetrically try second branch
    const auto mat0 = mat;
    Material mat1 { mat.begin(), mat.end() };

    Labelling l0, l1;
    const auto cs0 = cs.subspan(0, b/2);
    const auto cs1 = cs.subspan(b/2);
    const auto muxMatSize0 = (bL-1) * (cs0[0].nOut + 2);
    const auto muxMatSize1 = ((b - bL)-1) * (cs1[0].nOut + 2);
    const auto muxMat0 = muxMat.subspan(0, muxMatSize0);
    const auto muxMat1 = muxMat.subspan(muxMatSize0, muxMatSize1);
    const auto muxMatNow = muxMat.subspan(muxMatSize0 + muxMatSize1);



    // garble straight into material so as to unstack
    {
      gbBranches(f, cs1, seeds[2*bL-1], mat0);
      l0 = evBranches(f, cs0, seeds.subspan(1, 2*bL-1), inp0, mat0, muxMat0);
    }
    {
      gbBranches(f, cs0, seeds[0], mat1);
      l1 = evBranches(f, cs1, seeds.subspan(2*bL), inp1, mat1, muxMat1);
    }

    return evMux(f, S, l0, l1, muxMatNow);
  }
}


Interface gbCond(
    const PRF& f,
    std::span<const Circuit> cs,
    std::span<const Label> goodSeeds,
    std::span<const Label> badSeeds,
    std::span<Label> mat,
    std::span<Label> muxMat) {
  const auto b = cs.size();
  PRG prg(goodSeeds[0]);
  if (b == 1) {
    return garble(prg, f, cs[0], mat);
  } else {
    // generate child seeds and discard to advance PRG
    prg();
    prg();


    const auto bL = b/2;
    const auto cs0 = cs.subspan(0, bL);
    const auto cs1 = cs.subspan(bL);

    const auto n = cs[0].nInp + ilog2(cs.size());

    // TODO mat management
    const auto demMat = mat;

    gbBranches(f, cs0, badSeeds[0], mat);
    gbBranches(f, cs1, badSeeds[2*bL-1], mat);

    const auto i0 = gbCond(f, cs0, goodSeeds.subspan(1, 2*bL-1), badSeeds.subspan(1, 2*bL-1), mat, muxMat);
    const auto i1 = gbCond(f, cs1, goodSeeds.subspan(2*bL), badSeeds.subspan(2*b), mat, muxMat);

    auto e = genEncoding(prg, n);
    std::span<Label> zeros { e.zeros };
    gbDem(prg, f, e.delta, zeros[0], zeros.subspan(1), i0.inpEnc, i1.inpEnc, demMat);
  }
}








/* Labelling evCond( */
/*     const PRF& f, */
/*     std::span<const Circuit> cs, */
/*     std::span<const Label> seeds, */
/*     const Labelling& inp, */
/*     std::span<Label> mat, */
/*     std::span<Label> muxMat) { */

/*   const auto b = cs.size(); */

/*   if (b == 1) { */
/*     return ev(f, cs[0], inp, mat); */
/*   } else { */
/*     const auto bL = b/2; */

/*     // split input into the branch condition and the remaining wires */
/*     const auto S = inp[0]; */
/*     std::span<const Label> inp_(inp); */
/*     inp_ = inp_.subspan(1); */

/*     // evaluate the demux to compute inputs for both branches */
/*     const auto demOut = evDem(f, S, inp_, mat); */
/*     const auto inp0 = demOut.first; */
/*     const auto inp1 = demOut.second; */
/*     // move past the demux material */
/*     mat = mat.subspan(3*inp_.size() + 2); */

/*     // split the vector of circuits */
/*     const std::span<const Circuit> cs0 = cs.subspan(0, b/2); */
/*     const std::span<const Circuit> cs1 = cs.subspan(b/2); */
/*     const auto muxMatSize0 = (cs0.size()-1) * (cs0[0].nOut + 2); */
/*     const auto muxMatSize1 = (cs1.size()-1) * (cs1[0].nOut + 2); */
/*     const auto muxMat0 = muxMat.subspan(0, muxMatSize0); */
/*     const auto muxMat1 = muxMat.subspan(muxMatSize0, muxMatSize1); */
/*     const auto muxMatNow = muxMat.subspan(muxMatSize0 + muxMatSize1); */

/*     // we must make one copy of the material so we can symmetrically try second branch */
/*     const auto mat0 = mat; */
/*     Material mat1 { mat.begin(), mat.end() }; */

/*     Labelling l0, l1; */

/*     // garble straight into material so as to unstack */
/*     std::thread th0 ([&] { */
/*     /1* { *1/ */
/*       gbCond_(f, cs1, seeds[0], mat0); */
/*       l0 = evCond(f, cs0, seeds.subspan(1, 2*bL-1), inp0, mat0, muxMat0); */
/*     }); */
/*     /1* } *1/ */

/*     { */
/*       gbCond_(f, cs0, seeds[2*bL-1], mat1); */
/*       l1 = evCond(f, cs1, seeds.subspan(2*bL), inp1, mat1, muxMat1); */
/*     } */

/*     th0.join(); */

/*     return evMux(f, S, l0, l1, muxMatNow); */
/*   } */
/* } */


// Garble a vector of conditionally composed circuits, starting from a seed,
// into the specified material.
//
// While the material for most of the circuit is stacked with XOR,
// the multiplexer that collects garbage is not.
// Thus, we reference the main material and the mux material by different spans.
/* Interface gbCond( */
/*     const PRF& f, */
/*     std::span<const Circuit> cs, */
/*     std::span<const Label> goodSeeds, */
/*     std::span<const Label> badSeeds, */
/*     std::span<Label> mat, */
/*     std::span<Label> muxMat) { */
/*   const auto b = cs.size(); */
/*   PRG prg(goodSeeds[0]); */

/*   if (b == 1) { */
/*     return garble(prg, f, cs[0], mat); */
/*   } else { */
/*     const auto bL = b/2; */

/*     // construct and dispose child seeds (which are already in goodSeeds) */
/*     prg(); */
/*     prg(); */

/*     const auto n = cs[0].nInp + ilog2(cs.size()); */
/*     auto e = genEncoding(prg, n); */
/*     const auto S0 = e.zeros[0]; */

/*     // split the vector of circuits */
/*     const std::span<const Circuit> cs0 = cs.subspan(0, b/2); */
/*     const std::span<const Circuit> cs1 = cs.subspan(b/2); */

/*     // skip the demux material for now */
/*     const auto demSize = 3*(n-1) + 2; */
/*     const auto demMat = mat.subspan(0, demSize); */
/*     const auto branchmat = mat.subspan(demSize); */

/*     const auto muxMatSize0 = (cs0.size()-1) * (cs0[0].nOut + 2); */
/*     const auto muxMatSize1 = (cs1.size()-1) * (cs1[0].nOut + 2); */
/*     const auto muxMat0 = muxMat.subspan(0, muxMatSize0); */
/*     const auto muxMat1 = muxMat.subspan(muxMatSize0, muxMatSize1); */
/*     const auto muxMatNow = muxMat.subspan(muxMatSize0 + muxMatSize1); */

/*     Interface i0, i1; */
/*     Material mat0(branchmat.size()); */
/*     Material mat1(branchmat.size()); */
/*     // because garbling the conditional recursively involves unstacking and evaluating, */
/*     // we cannot garble the material in place, and instead must garble into a fresh buffer. */
/*     std::thread th ([&] { */
/*     /1* { *1/ */
/*       i0 = gbCond(f, cs0, goodSeeds.subspan(1, 2*bL-1), badSeeds.subspan(1, 2*bL-1), mat0, muxMat0); */
/*     }); */
/*     /1* } *1/ */
/*     { */
/*       i1 = gbCond(f, cs1, goodSeeds.subspan(2*bL), badSeeds.subspan(2*bL), mat1, muxMat1); */
/*     } */
/*     th.join(); */

/*     branchmat ^= mat0; */
/*     branchmat ^= mat1; */

/*     // now, with the input encodings available, we can garble the demux */
/*     // into the front of the material */
/*     std::span<Label> zeros { e.zeros }; */
/*     zeros = zeros.subspan(1); */
/*     auto bad = gbDem(prg, f, e.delta, S0, zeros, i0.inpEnc, i1.inpEnc, demMat); */
/*     auto bad0 = bad.first; */
/*     auto bad1 = bad.second; */

/*     Encoding e0_, e1_; */
/*     // copy the stacked material so as not to trash it */
/*     // then garble a branch incorrectly to emulate bad evaluation of the other */
/*     // branch (note garbling in place is safe because gbCond_ does not */
/*     // recursively involve evaluation). */
/*     std::thread th2 ([&] { */
/*     /1* { *1/ */
/*       Material mat0(branchmat.begin(), branchmat.end()); */
/*       e0_ = gbCond_(f, cs1, badSeeds[0], mat0); */
/*       bad0 = evCond(f, cs0, badSeeds.subspan(1, 2*bL-1), bad0, mat0, muxMat0); */
/*     }); */
/*   /1* } *1/ */
/*     { */
/*       Material mat1(branchmat.begin(), branchmat.end()); */
/*       e1_ = gbCond_(f, cs0, badSeeds[2*bL-1], mat1); */
/*       bad1 = evCond(f, cs1, badSeeds.subspan(1, 2*bL-1), bad1, mat1, muxMat1); */
/*     } */
/*     th2.join(); */


/*     const auto eout = gbMux( */
/*         prg, f, e.delta, S0, */
/*         i0.outEnc, */
/*         i1.outEnc, */
/*         bad0, */
/*         bad1, */
/*         muxMatNow); */

/*     return { e, eout }; */
/*   } */
/* } */


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
    const auto l = prg();
    const auto r = prg();

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
      PRG prg;
      const auto seed = prg();

      const auto b = cond.cs.size();
      const auto n = c.nInp;
      const auto m = c.nOut;

      const auto gadgetSize = 3*b;
      const auto transSize = n+1;
      const auto muxSize = (b-1)*(m+2);

      auto gadgetMat = mat.subspan(0, gadgetSize);
      const auto transMat = mat.subspan(gadgetSize, transSize);
      const auto bodyMat = mat.subspan(
          transSize + gadgetSize,
          mat.size() - transSize - muxSize - gadgetSize);
      const auto muxMat = mat.subspan(mat.size() - muxSize);


      std::span<const Label> inp = inpEnc.zeros;
      const auto goodSeeds = seedTree(seed, b);
      const auto badSeeds = gbGadget(b, f, inpEnc.delta, goodSeeds, inp.subspan(log2(b)), gadgetMat);

      const auto interface = gbBranches(f, cond.cs, seed, mat);

      gbTrans(prg, inpEnc, interface.inpEnc, transMat);

      return Encoding { };
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


Interface garble(PRG& seed, const PRF& f, const Circuit& c, std::span<Label> mat) {
  Interface g;
  g.inpEnc = genEncoding(seed, c.nInp);
  g.outEnc = gb(f, c, g.inpEnc, mat);
  return g;
}


Labelling ev(const PRF& f, const Circuit& c, const Labelling& input, std::span<Label> mat) {
  return std::visit(overloaded {
    [&](const Netlist& n) {
      return netlistev(f, n, input, mat);
    },
    [&](const Conditional& cond) {
      const auto b = cond.cs.size();
      const auto n = c.nInp;
      const auto m = c.nOut;

      const auto gadgetSize = 3*b;
      const auto transSize = n+1;
      const auto muxSize = (b-1)*(m+2);
      const auto gadgetMat = mat.subspan(0, gadgetSize);
      const auto transMat = mat.subspan(gadgetSize, transSize);
      const auto bodyMat = mat.subspan(
          transSize + gadgetSize,
          mat.size() - transSize - muxSize - gadgetSize);
      const auto muxMat = mat.subspan(mat.size() - muxSize);

      std::span<const Label> inps (input);
      const auto seeds = evGadget(b, f, inps.subspan(0, log2(b)), mat);


      const auto translated = evTrans(input, transMat);
      evBranches(f, cond.cs, seeds, translated, bodyMat, muxMat);

      //TODO
      return Labelling { };
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
