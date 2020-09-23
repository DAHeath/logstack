#ifndef SCHEME_H__
#define SCHEME_H__


#include <vector>
#include <variant>
#include <memory>
#include <bitset>
#include <span>

#include "prf.h"
#include "prg.h"

inline constexpr std::size_t ilog2(std::size_t n) {
  std::size_t out = 0;
  n *= 2;
  n -= 1;
  while (n >>= 1) { ++out; }
  return out;
}


enum class GateType { INPUT, OUTPUT, AND, XOR, NOT };


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


struct Interface {
  Encoding inputEncoding;
  Encoding outputEncoding;
};


using Wiring = std::vector<Label>;


struct NetlistCtxt {
  Wiring w;
  std::span<const Label> inp;
  std::span<Label> out;
  std::span<Label> material;
  std::size_t nonce;
};


Interface garble(PRG& seed, const PRF&, const Circuit&, std::span<Label>);
Encoding gb(const PRF&, const Circuit&, const Encoding&, std::span<Label>);
Labelling ev(const PRF&, const Circuit&, const Labelling&, std::span<Label>);

void gbGate(const PRF&, const Gate&, const Label& delta, NetlistCtxt&);
void evGate(const PRF&, const Gate&, NetlistCtxt&);
Encoding genEncoding(PRG&, std::size_t);

Labelling evCond(
    const PRF&,
    std::span<Circuit>,
    const Labelling&,
    std::span<Label> mat,
    std::span<Label> muxMat);

Encoding gbCond_(
    const PRF&,
    std::span<const Circuit>,
    const Label& seed,
    std::span<Label>);

Interface gbCond(
    const PRF&,
    std::span<const Circuit>,
    const Label& seed,
    std::span<Label> mat,
    std::span<Label> muxMat);

Encoding gbMux(
    PRG& prg,
    const PRF& f,
    const Label& delta,
    const Label& S0,
    const Encoding& good0,
    const Encoding& good1,
    const Labelling& bad0,
    const Labelling& bad1,
    std::span<Label> mat);

Labelling evMux(
    const PRF& f,
    const Label& S,
    const Labelling& X,
    const Labelling& Y,
    std::span<Label> mat);

std::pair<Labelling, Labelling> gbDem(
    PRG&,
    const PRF&,
    const Label& delta,
    const Label& S0,
    std::span<const Label>,
    const Encoding&,
    const Encoding&,
    std::span<Label>);

std::pair<Labelling, Labelling> evDem(
    const PRF& f,
    const Label& S,
    std::span<const Label>,
    std::span<Label>);

#endif
