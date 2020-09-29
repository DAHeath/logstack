#include "gate.h"


Label gbAnd(
    const PRF& f,
    const Label& delta,
    const Label& A0,
    const Label& B0,
    std::span<Label>& mat,
    std::size_t& nonce) {
  const auto nonce0 = Label { nonce };
  const auto nonce1 = Label { nonce + 1 };
  const auto hA0 = f(A0 ^ nonce0);
  const auto hA1 = f(A0 ^ delta ^ nonce0);
  const auto hB0 = f(B0 ^ nonce1);
  const auto hB1 = f(B0 ^ delta ^ nonce1);
  const auto A0D = A0 & delta;
  const auto B0D = B0 & delta;

  const auto X = A0[0] ? (hA1 ^ B0D) : hA0;
  const auto Y = B0[0] ? (hB1 ^ A0D) : hB0;

  mat[0] ^= hA0 ^ hA1 ^ B0D;
  mat[1] ^= hB0 ^ hB1 ^ A0D;

  nonce += 2;
  mat = mat.subspan(2);

  return A0&B0 ^ X ^ Y;
}


Label evAnd(
    const PRF& f,
    const Label& A,
    const Label& B,
    std::span<Label>& mat,
    std::size_t& nonce) {
  const auto nonce0 = Label { nonce };
  const auto nonce1 = Label { nonce + 1 };
  const auto hA = f(A ^ nonce0);
  const auto hB = f(B ^ nonce1);
  const auto X = A[0] ? hA ^ mat[0] : hA;
  const auto Y = B[0] ? hB ^ mat[1] : hB;
  mat = mat.subspan(2);
  nonce += 2;
  return (A&B) ^ X ^ Y;
}


