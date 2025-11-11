#pragma once

#include "oktal/geometry/Vec.hpp"
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace oktal {
using morton_bits_t = std::uint64_t;

class MortonIndex {
private:
  morton_bits_t _m_idx;
  // Private constructor used by the factory function
  // MortonIndex(morton_bits_t index) : _m_idx(index) {}

public:
  MortonIndex() : _m_idx(1) {}

  MortonIndex(morton_bits_t index) : _m_idx(index) {}

  [[nodiscard]] morton_bits_t getBits() const { return _m_idx; }

  static constexpr unsigned int MAX_DEPTH =
      ((std::numeric_limits<morton_bits_t>::digits - 1) / 3);

  static auto fromPath(const std::vector<morton_bits_t> &path) {
    if (path.size() > MAX_DEPTH) {
      throw std::invalid_argument(
          "Path is too long. Maximum depth for a morton_bits_t is " +
          std::to_string(MAX_DEPTH) + ", but path length is " +
          std::to_string(path.size()));
    }
    morton_bits_t res_idx = 1;
    for (const morton_bits_t &p : path) {
      if (p > 7) {
        throw std::invalid_argument("Invalid child choice " +
                                    std::to_string(p) +
                                    ". Must be between 0 and 7.");
      }
      res_idx = (res_idx << 3) | p;
    }
    return MortonIndex(res_idx);
  }

  [[nodiscard]] std::vector<morton_bits_t> getPath() const {
    std::vector<morton_bits_t> path;
    morton_bits_t tmp_idx = this->_m_idx;
    while (tmp_idx > 1) {
      const morton_bits_t p = tmp_idx & 0b111;
      path.insert(path.begin(), p);
      tmp_idx >>= 3;
    }
    return path;
  }

  [[nodiscard]] std::size_t level() const {
    morton_bits_t tmp_idx = this->_m_idx;
    size_t count = 0;
    while (tmp_idx > 1) {
      count++;
      tmp_idx >>= 3;
    }
    return count;
  }
  [[nodiscard]] bool isRoot() const { return this->_m_idx == 1; }

  [[nodiscard]] std::size_t siblingIndex() const {
    if (isRoot()) {
      return 0;
    }
    return static_cast<size_t>(this->_m_idx & 0b111);
  }

  [[nodiscard]] bool isFirstSibling() const {
    return isRoot() || (siblingIndex() == 0);
  }
  [[nodiscard]] bool isLastSibling() const {
    return isRoot() || (siblingIndex() == 7);
  }

  [[nodiscard]] auto parent() const { return MortonIndex(this->_m_idx >> 3); }

  [[nodiscard]] auto safeParent() const {
    if (isRoot()) {
      throw std::logic_error("Root cell has no parent.");
    }
    return parent();
  }

  [[nodiscard]] auto child(morton_bits_t idx) const {
    // Note: No check for depth limit, returns undefined value if exceeded.
    return MortonIndex((this->_m_idx << 3) | (idx & 0b111));
  }

  // safeChild is complex due to validation -> Declared here, defined in .cpp
  [[nodiscard]] auto safeChild(morton_bits_t idx) const -> MortonIndex;

  // --- Operators (Defined Inline) ---

  bool operator==(const MortonIndex &other) const {
    return this->_m_idx == other._m_idx;
  }
  bool operator!=(const MortonIndex &other) const { return !(*this == other); }

  // Partial Order Operators (Declared here, defined in .cpp due to complexity)
  bool operator>(const MortonIndex &y) const;
  bool operator<(const MortonIndex &y) const;
  bool operator>=(const MortonIndex &y) const {
    return (*this == y) || (*this > y);
  }
  bool operator<=(const MortonIndex &y) const {
    return (*this == y) || (*this < y);
  }

  [[nodiscard]] Vec<std::size_t, 3> gridCoordinates() const;

  static MortonIndex fromGridCoordinates(std::size_t level,
                                         const Vec<std::size_t, 3> &coords);
};

} // namespace oktal
