#ifndef SHA256_H__
#define SHA256_H__

#include <cstdint>
#include <cstring>
#include <cstdio>
#include "compiler.h"

std::array<U32, 8> sha256(std::array<U32, 16> input) {
  const std::array<U32, 64> k = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
  };

  const std::array<U32, 16> final_chunk = {
    0x80000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000200
  };

  const auto right_rot = [](U32 value, size_t count) {
    return value >> count | value << (32 - count);
  }; 
  const auto compress = [&](std::array<U32, 16> w, std::array<U32, 8> h) {
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 16; ++j) {
        if (i > 0) {
          /* Extend the first 16 words into the remaining 48 words w[16..63] of the message schedule std::array: */
          const auto index1 = w[(j + 1) & 0xf];
          const auto index14 = w[(j + 14) & 0xf];
          const U32 s0 = right_rot(index1, 7) ^ right_rot(index1, 18) ^ (index1 >> 3);
          const U32 s1 = right_rot(index14, 17) ^ right_rot(index14, 19) ^ (index14 >> 10);
          w[j] += w[(j + 9) & 0xf] + (s0 + s1);
        }
        const U32 s1 = right_rot(h[4], 6) ^ right_rot(h[4], 11) ^ right_rot(h[4], 25);
        const U32 ch = (h[4] & h[5]) ^ (~h[4] & h[6]);
        const U32 temp1 = h[7] + s1 + ch + k[i << 4 | j] + w[j];

        const U32 s0 = right_rot(h[0], 2) ^ right_rot(h[0], 13) ^ right_rot(h[0], 22);
        const U32 maj = (h[0] & (h[1] ^ h[2])) ^ (h[1] & h[2]);
        const U32 temp2 = s0 + maj;

        h[7] = h[6];
        h[6] = h[5];
        h[5] = h[4];
        h[4] = h[3] + temp1;
        h[3] = h[2];
        h[2] = h[1];
        h[1] = h[0];
        h[0] = temp1 + temp2;
      }
    }
    return h;
  };

  const auto add = [](const std::array<U32, 8>& xs, const std::array<U32, 8>& ys) {
    std::array<U32, 8> out;
    for (size_t i = 0; i < 8; ++i) { out[i] = xs[i] + ys[i]; }
    return out;
  };

  const std::array<U32, 8> h_init = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 };
  const auto h1 = add(h_init, compress(input, h_init));
  return add(h1, compress(final_chunk, h1));
}

#endif