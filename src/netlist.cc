#include "netlist.h"
#include "gate.h"


thread_local Gen::Ctxt Gen::ctxt;
thread_local Eval::Ctxt Eval::ctxt;
thread_local NetlistDesc Count::ctxt;


Gen::Rep Gen::Rep::operator&(const Gen::Rep& o) const {
  const auto& delta = ctxt.delta;
  const auto& f = ctxt.f;
  const auto& A0 = label;
  const auto& B0 = o.label;

  return gbAnd(f, delta, A0, B0, ctxt.material, ctxt.nonce);
}


Eval::Rep Eval::Rep::operator&(const Eval::Rep& o) const {
  const auto& f = ctxt.f;
  const auto& A = label;
  const auto& B = o.label;

  return evAnd(f, A, B, ctxt.material, ctxt.nonce);
}


Encoding netlistgb(const PRF& prf, const Netlist& c, const Encoding& inpEnc, std::span<Label> mat) {
  Encoding outEnc;
  outEnc.zeros.resize(c.desc.nOut);

  const auto stash = Gen::ctxt;

  Gen::ctxt.material = mat;
  Gen::ctxt.delta = inpEnc.delta;
  Gen::ctxt.inps = inpEnc.zeros;
  Gen::ctxt.outs = outEnc.zeros;
  Gen::ctxt.nonce = 0;
  Gen::ctxt.f = prf;

  c.gb();

  outEnc.delta = Gen::ctxt.delta;

  Gen::ctxt = stash;

  return outEnc;
}


Labelling netlistev(const PRF& prf, const Netlist& c, const Labelling& inp, std::span<Label> mat) {
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
