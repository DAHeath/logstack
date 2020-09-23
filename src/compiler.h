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


struct U32 {
public:
  U32() { }
  U32(std::uint32_t x) {
    std::bitset<32> xx = x;
    for (std::size_t i = 0; i < 32; ++i) { bits[i] = Bool { xx[i] }; }
  }

  static U32 input() {
    U32 out;
    for (auto& b: out.bits) { b = Bool::input(); }
    return out;
  }

  void output() const {
    for (const auto& b: bits) { b.output(); }
  }

  U32 operator&(const U32& o) const {
    U32 out;
    for (std::size_t i = 0; i < 32; ++i) { out.bits[i] = bits[i] & o.bits[i]; }
    return out;
  }

  U32 operator^(const U32& o) const {
    U32 out;
    for (std::size_t i = 0; i < 32; ++i) { out.bits[i] = bits[i] ^ o.bits[i]; }
    return out;
  }

  U32 operator|(const U32& o) const {
    U32 out;
    for (std::size_t i = 0; i < 32; ++i) { out.bits[i] = bits[i] | o.bits[i]; }
    return out;
  }

  U32 operator~() const {
    U32 out;
    for (std::size_t i = 0; i < 32; ++i) { out.bits[i] = ~bits[i]; }
    return out;
  }

  U32& operator&=(const U32& o) {
    for (std::size_t i = 0; i < 32; ++i) { bits[i] &= o.bits[i]; }
    return *this;
  }

  U32& operator^=(const U32& o) {
    for (std::size_t i = 0; i < 32; ++i) { bits[i] ^= o.bits[i]; }
    return *this;
  }

  U32& operator|=(const U32& o) {
    for (std::size_t i = 0; i < 32; ++i) { bits[i] |= o.bits[i]; }
    return *this;
  }

  Bool& operator[](std::size_t ix) { return bits[ix]; }
  const Bool& operator[](std::size_t ix) const { return bits[ix]; }

  U32 operator+(const U32& o) const {
    Bool carry { false };
    U32 out;
    for (size_t i = 0; i < 31; ++i) {
      out[i] = o[i] ^ bits[i] ^ carry;
      carry = carry ^ ((bits[i] ^ carry) & (o[i] ^ carry));
    }
    out[31] = o[31] ^ bits[31] ^ carry;
    return out;
  }

  U32& operator+=(const U32& o) {
    *this = (*this + o);
    return *this;
  }

  U32 operator<<(size_t sh) const {
    U32 out;
    if (sh > 32) {
      for (size_t i = 0; i < 32; ++i) { out[i] = false; }
    } else {
      for (int i = 31; i >= sh; --i) { out[i] = bits[i-sh]; }
      for (int i = sh-1; i >= 0; --i) { out[i] = false; }
    }
    return out;
  }

  U32 operator>>(size_t sh) const {
    U32 out;
    if (sh > 32) {
      for (size_t i = 0; i < 32; ++i) { out[i] = false; }
    } else {
      for (size_t i = sh; i < 32; ++i) { out[i-sh] = bits[i]; }
      for (size_t i = 32-sh; i < 32; ++i) { out[i] = false; }
    }
    return out;
  }

private:
  std::array<Bool, 32> bits;
};


#endif
