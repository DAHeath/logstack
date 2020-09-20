#include <vector>
#include <variant>
#include <memory>
#include <bitset>
#include <span>

template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

enum class GateType { INPUT, OUTPUT, AND, XOR };


struct Gate {
  GateType type;
  std::size_t inp0, inp1, out;
};


using Netlist = std::vector<Gate>;


struct Circuit;


struct Conditional {
  std::vector<Circuit> cs;
};


using Sequence = std::vector<Circuit>;


struct Circuit {
  std::variant<Netlist, Conditional, Sequence> content;
  std::size_t nWires;
  std::size_t nInp;
  std::size_t nOut;
  std::size_t nRow;
};



using Label = std::bitset<128>;
using Labelling = std::vector<Label>;


struct Encoding {
  Labelling zeros;
  Label delta;
};


using Material = std::vector<Label>;


struct Garbling {
  Encoding inputEncoding;
  Encoding outputEncoding;
  Material material;
};


using Wiring = std::vector<Label>;


void gbGate(
    const Gate& g,
    Wiring& w,
    const Label& delta,
    std::span<Label>& material,
    std::span<const Label>& inp,
    std::span<Label>& out) {
  switch (g.type) {
    case GateType::INPUT:
      // read in the next label from the input labelling
      w[g.out] = inp[0];
      inp = inp.subspan(1);
      break;
    case GateType::OUTPUT:
      // write out to the end of the output labelling
      out[0] = w[g.inp0];
      out = out.subspan(1);
      break;
    case GateType::AND:
      // TODO
      break;
    case GateType::XOR:
      w[g.out] = w[g.inp0] ^ w[g.inp1];
      break;
  }
}


void evGate(
    const Gate& g,
    Wiring& w,
    std::span<Label>& material,
    std::span<const Label>& inp,
    std::span<Label>& out) {
  switch (g.type) {
    case GateType::INPUT:
      // read in the next label from the input labelling
      w[g.out] = inp[0];
      inp = inp.subspan(1);
      break;
    case GateType::OUTPUT:
      // write out to the end of the output labelling
      out[0] = w[g.inp0];
      out = out.subspan(1);
      break;
    case GateType::AND:
      // TODO
      break;
    case GateType::XOR:
      w[g.out] = w[g.inp0] ^ w[g.inp1];
      break;
  }
}

Garbling garble(const Circuit&, const Label& seed);

Encoding gb(const Circuit& c, const Encoding& inputEncoding, std::span<Label>& material) {
  return std::visit(overloaded {
    [&](const Netlist& n) {
      const auto delta = inputEncoding.delta;
      Encoding outputEncoding;
      outputEncoding.delta = delta;
      outputEncoding.zeros.resize(c.nOut);

      Wiring wiring(c.nWires);
      std::span<const Label> inps(inputEncoding.zeros);
      std::span<Label> outs(outputEncoding.zeros);

      for (const auto& g: n) {
        gbGate(g, wiring, delta, material, inps, outs);
      }
      return outputEncoding;
    },
    [&](const Conditional& cond) {
      Encoding e;

      // TODO gen demux
      // TODO gen mux
      std::vector<Garbling> gs;
      gs.reserve(cond.cs.size());
      for (const auto& c: cond.cs) {
        gs.push_back(garble(c, 0)); // TODO properly define seed based on condition
      }
      return e;
    },
    [&](const Sequence& seq) {
      auto encoding = inputEncoding;
      for (const auto& c: seq) {
        encoding = gb(c, encoding, material);
      }
      return encoding;
    },
  }, c.content);
}


Garbling garble(const Circuit& c, const Label& seed) {
  Garbling g;
  // TODO generate by pulling from PRG
  // g.inputEncoding.delta = PRG
  // g.inputEncoding[0] = 1;
  // g.inputEncoding.zeros = std::generate 

  g.material.resize(c.nRow);
  std::span<Label> mat(g.material);
  g.outputEncoding = gb(c, g.inputEncoding, mat);
  return g;
}


Labelling ev(const Circuit& c, const Labelling& input, std::span<Label>& material) {
  return std::visit(overloaded {
    [&](const Netlist& n) {
      Labelling output(c.nOut);

      Wiring wiring(c.nWires);
      std::span<const Label> inps(input);
      std::span<Label> outs(output);

      for (const auto& g: n) {
        evGate(g, wiring, material, inps, outs);
      }
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
        labelling = ev(c, labelling, material);
      }
      return labelling;
    },
  }, c.content);
}


int main() {
}
