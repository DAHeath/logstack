#ifndef GEN_H__
#define GEN_H__


#include "prf.h"
#include <vector>
#include <span>


using Label = std::bitset<128>;


template <typename Val>
struct Boolean {
public:
  Boolean() : isConstant(false), constant(false) { }
  Boolean(bool b) : isConstant(true), constant(b) { }


  static Boolean input() { return Boolean(Val::input()); }
  void output() {
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
      return { false };
    } else if (o.isConstant) {
      return *this;
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

  Bool operator~() const { return *this ^ Bool { true }; }

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
  Boolean(const Val& val) : isConstant(false), val(val) { }

  bool isConstant;
  bool constant;
  Val val;
};


struct Gen {
  struct Ctxt {
    std::span<Label> material;
    Label delta;
    std::size_t nonce;
    PRF f;
  };

  static thread_local Ctxt ctxt;

  struct Bool {
  public:
    Bool() : label(0), isConstant(false) { }
    Bool(bool b) : label(b ? 1 : 0), isConstant(true) { }

    Bool operator&(const Bool&) const;
    Bool operator^(const Bool&) const;

    Bool operator~() const { return *this ^ Bool { true }; }
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
    Bool(const Label& label, bool isConstant) : label(label), isConstant(isConstant) { }

    Label label;
    bool isConstant;
  };
};


struct Eval {
  struct Ctxt {
    std::span<Label> material;
    std::size_t nonce;
    PRF f;
  };

  static thread_local Ctxt ctxt;

  struct Bool {
  public:
    Bool() : label(0), isConstant(false) { }
    Bool(bool b) : label(b ? 1 : 0), isConstant(true) { }

    Bool operator&(const Bool&) const;
    Bool operator^(const Bool&) const;

    Bool operator~() const { return *this ^ Bool { true }; }

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
    Bool(const Label& label, bool isConstant) : label(label), isConstant(isConstant) { }

    Label label;
    bool isConstant;
  };
};


struct Count {
  struct Ctxt {
    std::size_t nInp;
    std::size_t nOut;
    std::size_t nRow;
  };

  static thread_local Ctxt ctxt;

  struct Bool {
  public:
    Bool() : isConstant(false), constant(false) { }
    Bool(bool b) : isConstant(true), constant(b) { }

    Bool operator&(const Bool&) const;
    Bool operator^(const Bool&) const;

    Bool operator~() const { return *this ^ Bool { true }; }

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
    Bool(bool isConstant, bool constant) : isConstant(isConstant), constant(constant) { }

    bool isConstant;
    bool constant;
  };
};


#endif
