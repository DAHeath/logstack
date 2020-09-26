#ifndef GEN_H__
#define GEN_H__


#include "prf.h"
#include "prg.h"
#include <vector>
#include <span>


using Label = std::bitset<128>;


template <typename Rep>
struct Boolean {
public:
  Boolean() : isConstant(false), constant(false) { }
  Boolean(bool b) : isConstant(true), constant(b) { }


  static Boolean input() { return Boolean(Rep::input()); }
  void output() const {
    assert(!isConstant);
    val.output();
  }

  Boolean operator^(const Boolean& o) const {
    if (isConstant && o.isConstant) {
      return { constant != o.constant };
    } else if (isConstant && !constant) {
      return o;
    } else if (isConstant) {
      return { ~o.val };
    } else if (o.isConstant && !o.constant) {
      return *this;
    } else if (o.isConstant) {
      return { ~val };
    } else {
      return { val ^ o.val };
    }
  }

  Boolean operator&(const Boolean& o) const {
    if (isConstant && o.isConstant) {
      return { constant && o.constant };
    } else if (isConstant && !constant) {
      return { false };
    } else if (isConstant) {
      return o;
    } else if (o.isConstant && !o.constant) {
      return { false };
    } else if (o.isConstant) {
      return *this;
    } else {
      return { val & o.val };
    }
  }

  Boolean operator~() const { return *this ^ Boolean { true }; }

  Boolean operator|(const Boolean& o) const { return ~(~(*this) & ~o); }

  Boolean& operator^=(const Boolean& o) {
    *this = *this ^ o;
    return *this;
  }

  Boolean& operator&=(const Boolean& o) {
    *this = *this & o;
    return *this;
  }

  Boolean& operator|=(const Boolean& o) {
    *this = *this | o;
    return *this;
  }


private:
  Boolean(const Rep& val) : isConstant(false), val(val) { }

  bool isConstant;
  bool constant;
  Rep val;
};


struct Gen {
  struct Ctxt {
    std::span<Label> material;
    std::span<Label> inps;
    std::span<Label> outs;
    Label delta;
    std::size_t nonce;
    PRF f;
  };

  static thread_local Ctxt ctxt;

  struct Rep {
  public:
    static Rep input() {
      Label l = ctxt.inps[0];
      ctxt.inps = ctxt.inps.subspan(1);
      return { l };
    }

    void output() const {
      ctxt.outs[0] = label;
      ctxt.outs = ctxt.outs.subspan(1);
    }

    Rep() : label(0) { }

    Rep operator&(const Rep&) const;
    Rep operator^(const Rep& o) const { return { label ^ o.label }; }
    Rep operator~() const { return { label ^ ctxt.delta }; }

  private:
    Rep(const Label& label) : label(label) { }

    Label label;
  };

  using Bool = Boolean<Rep>;
};


struct Eval {
  struct Ctxt {
    std::span<Label> material;
    std::span<const Label> inps;
    std::span<Label> outs;
    std::size_t nonce;
    PRF f;
  };

  static thread_local Ctxt ctxt;

  struct Rep {
  public:
    static Rep input() {
      Label l = ctxt.inps[0];
      ctxt.inps = ctxt.inps.subspan(1);
      return { l };
    }

    void output() const {
      ctxt.outs[0] = label;
      ctxt.outs = ctxt.outs.subspan(1);
    }

    Rep() : label(0) { }

    Rep operator&(const Rep&) const;
    Rep operator^(const Rep& o) const { return { label ^ o.label }; }
    Rep operator~() const { return { label }; }

  private:
    Rep(const Label& label) : label(label) { }

    Label label;
  };

  using Bool = Boolean<Rep>;
};


struct CircuitDesc {
  std::size_t nInp;
  std::size_t nOut;
  std::size_t nRow;
};


struct Count {
  static thread_local CircuitDesc ctxt;

  struct Rep {
    static Rep input() { ++ctxt.nInp; return Rep { }; }
    void output() const { ++ctxt.nOut; }

    Rep operator&(const Rep&) const { ctxt.nRow += 2; return Rep { }; }
    Rep operator^(const Rep&) const { return Rep { }; }
    Rep operator~() const { return Rep { }; }
  };

  using Bool = Boolean<Rep>;
};


struct Circuit {
  void (*gb)(void);
  void (*ev)(void);
  CircuitDesc desc;
};


template <typename Gb, typename Ev, typename Desc>
Circuit compile(const Gb& gb, const Ev& ev, const Desc& desc) {
  CircuitDesc tmp = Count::ctxt;
  Count::ctxt = { };
  desc();
  const auto d = Count::ctxt;
  Count::ctxt = tmp;

  return { gb, ev, d };
}


#define CIRCUIT(name) \
  template <typename Bool> \
  void name##IMPL() {

#define END_CIRCUIT(name) \
  } \
  const Circuit name = compile(name##IMPL<Gen::Bool>, name##IMPL<Eval::Bool>, name##IMPL<Count::Bool>);


using Labelling = std::vector<Label>;

struct Encoding {
  Labelling zeros;
  Label delta;
};


struct Interface {
  Encoding inpEnc;
  Encoding outEnc;
};


Interface circuitgb(PRG&, const PRF&, const Circuit&, std::span<Label>);
Labelling circuitev(const PRF&, const Circuit&, const Labelling&, std::span<Label>);


#endif
