#ifndef GEN_H__
#define GEN_H__


#include "prf.h"
#include <vector>


using Label = std::bitset<128>;


struct Gen {
  struct Ctxt {
    std::vector<Label> material;
    Label delta;
    std::size_t nonce;
    PRF f;
  };

  static thread_local Ctxt ctxt;

  struct Bool {
  public:
    Bool() : label(0) { }
    Bool(bool b) : label(b ? ctxt.delta : 0) { }
    Bool(const Label& label) : label(label) { }

    Bool operator&(const Bool&) const;

    Bool operator^(const Bool& o) const { return { label ^ o.label }; }
    Bool operator~() const { return { label ^ ctxt.delta }; }
    Bool operator|(const Bool& o) const { return ~(~(*this) & ~o); }

    Bool& operator^=(const Bool& o) {
      *this = *this ^ o;
      return *this;
    }

    Bool& operator&=(const Bool& o) {
      *this = *this & o;
      return *this;
    }

    Bool& operator|=(const Bool& o) {
      *this = *this | o;
      return *this;
    }

  private:
    Label label;
  };
};


#endif
