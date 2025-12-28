#pragma once
#include "oktal/geometry/Vec.hpp"
#include "oktal/octree/CellOctree.hpp"
#include <array>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace oktal {

template <typename T>
concept PeriodicityMapper =
    std::semiregular<T> &&
    requires(T const mapper, Vec<std::ptrdiff_t, 3> coords, std::size_t level) {
      {
        mapper(coords, level)
      } -> std::same_as<std::optional<Vec<std::size_t, 3>>>;
    };

class NoPeriodicity {
public:
  NoPeriodicity() = default;

  [[nodiscard]] std::optional<Vec<std::size_t, 3>>
  operator()(Vec<std::ptrdiff_t, 3> coords, std::size_t level) const {
    const std::size_t max_coord = (1ULL << level);
    if (coords[0] >= 0 && coords[0] < static_cast<std::ptrdiff_t>(max_coord) &&
        coords[1] >= 0 && coords[1] < static_cast<std::ptrdiff_t>(max_coord) &&
        coords[2] >= 0 && coords[2] < static_cast<std::ptrdiff_t>(max_coord)) {
      Vec<std::size_t, 3> result;
      result[0] = static_cast<std::size_t>(coords[0]);
      result[1] = static_cast<std::size_t>(coords[1]);
      result[2] = static_cast<std::size_t>(coords[2]);
      return result;
    }
    return std::nullopt;
  }
};

class Torus {
public:
  Torus() : periodicity_({false, false, false}) {}

  explicit Torus(std::array<bool, 3> periodicity) : periodicity_(periodicity) {}

  [[nodiscard]] std::optional<Vec<std::size_t, 3>>
  operator()(Vec<std::ptrdiff_t, 3> coords, std::size_t level) const {
    const auto max_coord = static_cast<std::ptrdiff_t>(1ULL << level);
    Vec<std::size_t, 3> result;
    for (std::size_t i = 0; i < 3; i++) {
      if (periodicity_[i]) {
        std::ptrdiff_t wrapped = coords[i] % max_coord;
        if (wrapped < 0) {
          wrapped += max_coord;
        }
        result[i] = static_cast<std::size_t>(wrapped);
      } else {
        if (coords[i] < 0 || coords[i] >= max_coord) {
          return std::nullopt;
        }
        result[i] = static_cast<std::size_t>(coords[i]);
      }
    }
    return result;
  }

private:
  std::array<bool, 3> periodicity_;
};

class CellGrid {
public:
  static constexpr std::size_t NOT_ENUMERATED =
      std::numeric_limits<std::size_t>::max();
  static constexpr std::size_t NO_NEIGHBOR = NOT_ENUMERATED;

  class CellGridBuilder;
  class CellView;
  class CellIterator;

  [[nodiscard]] const CellOctree &octree() const;
  [[nodiscard]] std::span<const MortonIndex> mortonIndices() const;
  static CellGridBuilder create(std::shared_ptr<const CellOctree> octree);

  [[nodiscard]] std::size_t getEnumerationIndex(std::size_t streamIndex) const;
  [[nodiscard]] std::size_t
  getEnumerationIndex(const CellOctree::CellView &cell) const;
  [[nodiscard]] std::span<const std::size_t>
  neighborIndices(Vec<std::ptrdiff_t, 3> offset) const;

  [[nodiscard]] CellView operator[](std::size_t index) const;

  [[nodiscard]] CellIterator begin() const;
  [[nodiscard]] CellIterator end() const;
  [[nodiscard]] std::size_t size() const;

private:
  CellGrid() = default;
  std::shared_ptr<const CellOctree> octree_;
  std::vector<MortonIndex> morton_indices_;
  std::map<std::size_t, std::size_t> stream_to_enum_;
  std::map<Vec<std::ptrdiff_t, 3>, std::vector<std::size_t>> adjacency_lists_;

  friend class CellGridBuilder;
};

using CellGridBuilder = CellGrid::CellGridBuilder;

class CellGrid::CellGridBuilder {
public:
  explicit CellGridBuilder(std::shared_ptr<const CellOctree> octree)
      : octree_(std::move(octree)), periodicity_mapper_(NoPeriodicity{}) {}

  CellGridBuilder &levels(std::span<const std::size_t> lvls);
  CellGridBuilder &levels(std::initializer_list<std::size_t> lvls);
  CellGridBuilder &neighborhood(std::span<const Vec<std::ptrdiff_t, 3>> offset);
  CellGridBuilder &
  neighborhood(std::initializer_list<Vec<std::ptrdiff_t, 3>> offset);

  [[nodiscard]] CellGrid build() const;

  template <PeriodicityMapper TMapper>
  CellGridBuilder &periodicityMapper(TMapper mapper) {
    periodicity_mapper_ = std::move(mapper);
    return *this;
  }

private:
  std::shared_ptr<const CellOctree> octree_;
  std::vector<std::size_t> levels_;
  std::vector<Vec<std::ptrdiff_t, 3>> neighborhood_;
  std::function<std::optional<Vec<std::size_t, 3>>(Vec<std::ptrdiff_t, 3>,
                                                   std::size_t)>
      periodicity_mapper_;
};

inline CellGrid::CellGridBuilder
CellGrid::create(std::shared_ptr<const CellOctree> octree) {
  return CellGridBuilder(std::move(octree));
}

class CellGrid::CellView {
public:
  CellView() = default;
  CellView(std::size_t enum_index, const CellGrid *grid)
      : enum_index_(enum_index), grid_(grid) {}

  [[nodiscard]] std::size_t enumerationIndex() const { return enum_index_; }

  [[nodiscard]] bool isValid() const {
    return enum_index_ != NOT_ENUMERATED && grid_ != nullptr;
  }

  [[nodiscard]] operator std::size_t() const { return enum_index_; }
  [[nodiscard]] explicit operator bool() const { return isValid(); }

  [[nodiscard]] Vec<double, 3> center() const;
  [[nodiscard]] Box<double> boundingBox() const;

  [[nodiscard]] MortonIndex mortonIndex() const;
  [[nodiscard]] std::size_t level() const;
  [[nodiscard]] CellView neighbor(Vec<std::ptrdiff_t, 3> offset) const;

private:
  std::size_t enum_index_ = NOT_ENUMERATED;
  const CellGrid *grid_ = nullptr;
};

class CellGrid::CellIterator {
public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = CellGrid::CellView;
  using difference_type = std::ptrdiff_t;
  using pointer = const CellGrid::CellView *;
  using reference = const CellGrid::CellView &;

  CellIterator() = default;

  CellIterator(const CellGrid *grid, std::size_t index)
      : grid_(grid), index_(index) {}

  [[nodiscard]] CellGrid::CellView operator*() const { return {index_, grid_}; }

  CellIterator &operator++() {
    ++index_;
    return *this;
  }

  CellIterator operator++(int) {
    CellIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  [[nodiscard]] bool operator==(const CellIterator &other) const {
    return grid_ == other.grid_ && index_ == other.index_;
  }

  [[nodiscard]] bool operator!=(const CellIterator &other) const {
    return !(*this == other);
  }

private:
  const CellGrid *grid_ = nullptr;
  std::size_t index_ = 0;
};

} // namespace oktal
