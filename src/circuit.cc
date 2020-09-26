#include "circuit.h"


thread_local Gen::Ctxt Gen::ctxt;
thread_local Eval::Ctxt Eval::ctxt;
thread_local CircuitDesc Count::ctxt;


Gen::Rep Gen::Rep::operator&(const Gen::Rep& o) const {
  const auto& delta = ctxt.delta;
  const auto& f = ctxt.f;
  const auto& A0 = label;
  const auto& B0 = o.label;

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

  return { A0&B0 ^ X ^ Y };
}


Eval::Rep Eval::Rep::operator&(const Eval::Rep& o) const {
  const auto& f = ctxt.f;
  const auto& A = label;
  const auto& B = o.label;

  const auto nonce0 = Label { ctxt.nonce };
  const auto nonce1 = Label { ctxt.nonce + 1 };
  const auto hA = f(A ^ nonce0);
  const auto hB = f(B ^ nonce1);
  const auto X = A[0] ? hA ^ ctxt.material[0] : hA;
  const auto Y = B[0] ? hB ^ ctxt.material[1] : hB;

  ctxt.material = ctxt.material.subspan(2);
  ctxt.nonce += 2;

  return { (A&B) ^ X ^ Y };
}


Interface gb(PRG& prg, const PRF& prf, const Circuit& c, std::span<Label> mat) {
  Interface interface;

  interface.inpEnc.delta = prg();
  interface.inpEnc.delta[0] = 1;
  interface.inpEnc.zeros.resize(c.desc.nInp);

  for (auto& i: interface.inpEnc.zeros) { i = prg(); }

  interface.outEnc.zeros.resize(c.desc.nOut);

  const auto stash = Gen::ctxt;

  Gen::ctxt.material = mat;
  Gen::ctxt.delta = interface.inpEnc.delta;
  Gen::ctxt.inps = interface.inpEnc.zeros;
  Gen::ctxt.outs = interface.outEnc.zeros;
  Gen::ctxt.nonce = 0;
  Gen::ctxt.f = prf;

  c.gb();

  interface.outEnc.delta = Gen::ctxt.delta;

  Gen::ctxt = stash;

  return interface;
}


Labelling ev(const PRF& prf, const Circuit& c, const Labelling& inp, std::span<Label> mat) {
  Labelling out(c.desc.nOut);

  const auto stash = Eval::ctxt;

  Eval::ctxt.material = mat;
  Eval::ctxt.inps = inp;
  Eval::ctxt.outs = out;
  Eval::ctxt.nonce = 0;
  Eval::ctxt.f = prf;

  c.ev();

  Eval::ctxt = stash;
  return out;
}
