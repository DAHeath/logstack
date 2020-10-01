#ifndef SCHEME_H__
#define SCHEME_H__


#include <vector>
#include <variant>
#include <memory>
#include <bitset>
#include <span>

#include "netlist.h"


template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;


inline constexpr std::size_t ilog2(std::size_t n) {
  std::size_t out = 0;
  n *= 2;
  n -= 1;
  while (n >>= 1) { ++out; }
  return out;
}


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



using Material = std::vector<Label>;


struct Garbling {
  Material material;
  Interface interface;
};

struct CondGarbling {
  Material material;
  Encoding outEnc;
  std::vector<Labelling> badOuts;
};

Encoding gb(const PRF&, const Circuit&, const Encoding&, std::span<Label>);
Labelling ev(const PRF&, const Circuit&, const Labelling&, std::span<Label>);

Encoding genEncoding(PRG&, std::size_t);

Labelling evCond(
    const PRF&,
    std::span<Circuit>,
    std::span<const Label> seeds,
    const Labelling&,
    std::span<Label> mat,
    std::span<Label> muxMat);

Interface gbCond(
    const PRF&,
    std::span<const Circuit>,
    std::span<const Label> goodSeeds,
    std::span<const Label> badSeeds,
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


std::vector<Label> gbGadget(
    std::size_t,
    const PRF&,
    const Label& delta,
    std::span<const Label> inseeds,
    std::span<const Label> index,
    std::span<Label>& mat);


std::vector<Label> evGadget(
    std::size_t b,
    const PRF& f,
    std::span<const Label> index,
    std::span<Label>& mat);


#endif
