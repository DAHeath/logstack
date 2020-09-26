#ifndef COMPILER_H__
#define COMPILER_H__

#include "scheme.h"
#include <array>
#include <bitset>


struct Bool {
public:
  Bool() : wire(0) { }
  Bool(bool b) : wire(b ? 1 : 0) { }

  static Bool input() {
    Bool out;
    out.wire = counter;
    ++counter;
    gates->push_back(Gate { GateType::INPUT, 0, 0, out.wire });
    return out;
  }

  Bool operator&(Bool o) const {
    if (wire < 2 && o.wire < 2) {
      return (bool)wire && (bool)o.wire;
    } else if (wire < 2) {
      return o & *this;
    } else if (o.wire == 1) {
      return *this;
    } else if (o.wire == 0) {
      return false;
    } else {
      Bool out;
      out.wire = counter;
      ++counter;
      gates->push_back(Gate { GateType::AND, wire, o.wire, out.wire });
      return out;
    }
  }

  Bool operator^(Bool o) const {
    if (wire < 2 && o.wire < 2) {
      return (bool)wire != (bool)o.wire;
    } else if (wire < 2) {
      return o ^ *this;
    } else if (o.wire == 1) {
      return ~(*this);
    } else if (o.wire == 0) {
      return *this;
    } else {
      Bool out;
      out.wire = counter;
      ++counter;
      gates->push_back(Gate { GateType::XOR, wire, o.wire, out.wire });
      return out;
    }
  }

  Bool operator~() const {
    if (wire == 0) {
      return true;
    } else if (wire == 1) {
      return false;
    } else {
      Bool out;
      out.wire = counter;
      ++counter;
      gates->push_back(Gate { GateType::NOT, wire, 0, out.wire });
      return out;
    }
  }

  Bool& operator&=(Bool o) {
    *this = *this & o;
    return *this;
  }

  Bool& operator^=(Bool o) {
    *this = *this ^ o;
    return *this;
  }

  Bool operator|(Bool o) const {
    return ~(~(*this) & ~o);
  }

  Bool& operator|=(Bool o) {
    *this = *this | o;
    return *this;
  }

  void output() const {
    gates->push_back(Gate { GateType::OUTPUT, wire, 0, 0 });
  }

  static Circuit compile() {
    Circuit c;
    c.content = *gates;
    c.nInp = 0;
    c.nOut = 0;
    c.nRow = 0;
    for (const auto& g: *gates) {
      switch (g.type) {
        case GateType::INPUT: ++c.nInp; break;
        case GateType::OUTPUT: ++c.nOut; break;
        case GateType::AND: c.nRow += 2; break;
        case GateType::XOR: break;
        case GateType::NOT: break;
      }
    }

    gates = new std::vector<Gate>();
    counter = 2;
    return c;
  }

private:
  static std::vector<Gate>* gates;
  static std::uint32_t counter;

  // we reserve the special wire values 0 and 1 for false/true respectively
  std::uint32_t wire;
};


#endif
