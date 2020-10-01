#ifndef PRG_H__
#define PRG_H__


#include "prf.h"


struct PRG {
public:
  PRG() : nonce(0) { }
  PRG(PRF prf) : prf(std::move(prf)), nonce(0), child_nonce(0-1) { }
  PRG(std::bitset<128> seed) : prf(std::move(seed)), nonce(0) { }

  std::bitset<128> operator()() { return prf(nonce++); }

  std::bitset<128> child() { return prf(child_nonce--); }

  void setNonce(std::size_t n) { nonce = n; }

private:
  PRF prf;
  std::size_t nonce;
  std::size_t child_nonce;
};


#endif
