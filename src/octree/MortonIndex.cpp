#include "oktal/octree/MortonIndex.hpp"

namespace oktal {

MortonIndex MortonIndex::safeChild(morton_bits_t idx) const {
  if (level() >= MAX_DEPTH) {
    throw std::logic_error("Non-Safe Child: Maximum depth (" +
                           std::to_string(MAX_DEPTH) + ") would be exceeded.");
  }
  return child(idx);
}

bool MortonIndex::operator>(const MortonIndex &y) const {
  const MortonIndex &x = *this;

  if (x == y) {
    return false;
  }

  const size_t x_level = x.level();
  const size_t y_level = y.level();

  if (x_level >= y_level) {
    return false;
  }

  const size_t shift_count = (y_level - x_level) * 3;

  return x._m_idx == (y._m_idx >> shift_count);
}

bool MortonIndex::operator<(const MortonIndex &y) const { return y > *this; }

Vec<std::size_t, 3> MortonIndex::gridCoordinates() const {
  Vec<std::size_t, 3> coords{};
  coords[0] = 0; // x
  coords[1] = 0; // y
  coords[2] = 0; // z

  morton_bits_t tmp_idx = this->_m_idx;
  std::size_t bit_pos = 0;

  // NOTE: Skip the leading 1, we're reading from right to left here.
  while (tmp_idx > 1) {
    coords[0] |= ((tmp_idx & 0b1) << bit_pos); // x bit
    tmp_idx >>= 1;
    coords[1] |= ((tmp_idx & 0b1) << bit_pos); // y bit
    tmp_idx >>= 1;
    coords[2] |= ((tmp_idx & 0b1) << bit_pos); // z bit
    tmp_idx >>= 1;
    bit_pos++;
  }

  return coords;
}

MortonIndex MortonIndex::fromGridCoordinates(std::size_t level,
                                      const Vec<std::size_t, 3> &coords) {
  if (level > MAX_DEPTH) {
    throw std::invalid_argument(
        "Level exceeds maximum depth of" + std::to_string(MAX_DEPTH) +
        ", requested level is " + std::to_string(level));
  }

  morton_bits_t result = 1;

	//Loop: Most significant to least significant
  for (std::size_t i = level; i > 0; --i) {
    const std::size_t bit_idx = i - 1;
    const morton_bits_t x_bit = (coords[0] >> bit_idx) & 1;
    const morton_bits_t y_bit = (coords[1] >> bit_idx) & 1;
    const morton_bits_t z_bit = (coords[2] >> bit_idx) & 1;

    // Shift result lefts and OR-ing interleaved bits
    result = (result << 3) | (z_bit << 0b10) | (y_bit << 0b1) | x_bit;
  }
  return {result};
}

} // namespace oktal
