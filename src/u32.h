#ifndef U32_H__
#define U32_H__

template <typename Bool>
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
