#include "murmurhash32.hpp"

int main() {
  uint64_t data = 123;
  uint32_t seed = 0;
  uint32_t hash32;

  // For 32-bit hash

  hash32 = murmurhash(&data, seed);
  std::cout << "Valor de hash de '"<< data << "': " << hash32 << "\n";
  return 0;
}
