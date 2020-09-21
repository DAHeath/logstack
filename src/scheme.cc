#include "scheme.h"


template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;


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
  g.inputEncoding.delta = seed();
  g.inputEncoding.delta[0] = 1;
  g.inputEncoding.zeros.resize(c.nInp);
  for (auto& z: g.inputEncoding.zeros) { z = seed(); }

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
      Labelling out;
      // TODO
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
