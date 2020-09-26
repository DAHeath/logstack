#include "gen.h"


thread_local Gen::Ctxt Gen::ctxt;
thread_local Eval::Ctxt Eval::ctxt;
thread_local Count::Ctxt Count::ctxt;


Gen::Bool Gen::Bool::operator^(const Gen::Bool& o) const {
  if (isConstant && o.isConstant) {
    return { label ^ o.label, true };
  } else if (isConstant && label == 0) {
    return o;
  } else if (isConstant) {
    return { o.label ^ ctxt.delta, false };
  } else if (o.isConstant && o.label == 0) {
    return *this;
  } else if (o.isConstant) {
    return { label ^ ctxt.delta, false };
  } else {
    return { label ^ o.label, false };
  }
}


Gen::Bool Gen::Bool::operator&(const Gen::Bool& o) const {
  if (isConstant && o.isConstant) {
    return { label & o.label, true };
  } else if (isConstant && label == 0) {
    return { false };
  } else if (isConstant) {
    return o;
  } else if (o.isConstant && o.label == 0) {
    return { false };
  } else if (o.isConstant) {
    return *this;
  } else {
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

    return { A0&B0 ^ X ^ Y, false };
  }
};


Eval::Bool Eval::Bool::operator^(const Eval::Bool& o) const {
  if (isConstant && o.isConstant) {
    return { label ^ o.label, true };
  } else if (isConstant) {
    return o;
  } else if (o.isConstant) {
    return *this;
  } else {
    return { label ^ o.label, false };
  }
}


Eval::Bool Eval::Bool::operator&(const Eval::Bool& o) const {
  if (isConstant && o.isConstant) {
    return { label & o.label, true };
  } else if (isConstant && label == 0) {
    return { false };
  } else if (isConstant) {
    return o;
  } else if (o.isConstant && o.label == 0) {
    return { false };
  } else if (o.isConstant) {
    return *this;
  } else {
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

    return { (A&B) ^ X ^ Y, false };
  }
}


Count::Bool Count::Bool::operator^(const Count::Bool& o) const {
}
