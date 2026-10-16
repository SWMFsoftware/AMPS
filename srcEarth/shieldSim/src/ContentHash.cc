#include "ContentHash.hh"

#include <array>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <vector>

namespace shield_sim {
namespace {

constexpr std::array<std::uint32_t, 64> k = {{
  0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,
  0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,
  0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
  0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,
  0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,
  0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
  0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,
  0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u
}};

std::uint32_t RotateRight(std::uint32_t value, unsigned count) {
  return (value >> count) | (value << (32u - count));
}

} // namespace

std::string Sha256(const std::string& bytes) {
  std::vector<std::uint8_t> data(bytes.begin(), bytes.end());
  const std::uint64_t bitLength = static_cast<std::uint64_t>(data.size()) * 8u;
  data.push_back(0x80u);
  while ((data.size() % 64u) != 56u) data.push_back(0u);
  for (int shift = 56; shift >= 0; shift -= 8)
    data.push_back(static_cast<std::uint8_t>((bitLength >> shift) & 0xffu));

  std::array<std::uint32_t, 8> hash = {{
    0x6a09e667u,0xbb67ae85u,0x3c6ef372u,0xa54ff53au,
    0x510e527fu,0x9b05688cu,0x1f83d9abu,0x5be0cd19u
  }};

  for (std::size_t block = 0; block < data.size(); block += 64u) {
    std::array<std::uint32_t, 64> w{};
    for (std::size_t i = 0; i < 16; ++i) {
      const std::size_t p = block + i * 4u;
      w[i] = (static_cast<std::uint32_t>(data[p]) << 24) |
             (static_cast<std::uint32_t>(data[p + 1]) << 16) |
             (static_cast<std::uint32_t>(data[p + 2]) << 8) |
              static_cast<std::uint32_t>(data[p + 3]);
    }
    for (std::size_t i = 16; i < 64; ++i) {
      const std::uint32_t s0 = RotateRight(w[i - 15], 7) ^ RotateRight(w[i - 15], 18) ^ (w[i - 15] >> 3);
      const std::uint32_t s1 = RotateRight(w[i - 2], 17) ^ RotateRight(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }

    std::uint32_t a = hash[0], b = hash[1], c = hash[2], d = hash[3];
    std::uint32_t e = hash[4], f = hash[5], g = hash[6], h = hash[7];
    for (std::size_t i = 0; i < 64; ++i) {
      const std::uint32_t big1 = RotateRight(e, 6) ^ RotateRight(e, 11) ^ RotateRight(e, 25);
      const std::uint32_t choose = (e & f) ^ ((~e) & g);
      const std::uint32_t temp1 = h + big1 + choose + k[i] + w[i];
      const std::uint32_t big0 = RotateRight(a, 2) ^ RotateRight(a, 13) ^ RotateRight(a, 22);
      const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t temp2 = big0 + majority;
      h = g; g = f; f = e; e = d + temp1;
      d = c; c = b; b = a; a = temp1 + temp2;
    }
    hash[0] += a; hash[1] += b; hash[2] += c; hash[3] += d;
    hash[4] += e; hash[5] += f; hash[6] += g; hash[7] += h;
  }

  std::ostringstream out;
  out << "sha256:" << std::hex << std::setfill('0');
  for (const std::uint32_t value : hash) out << std::setw(8) << value;
  return out.str();
}

} // namespace shield_sim
