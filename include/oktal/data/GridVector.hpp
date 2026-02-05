#pragma once
#include "oktal/octree/CellGrid.hpp"
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <experimental/mdspan>
#include <memory>
#include <vector>

namespace oktal {

using namespace std::experimental;
// Remove the const so for both T and const T the concept is furfilled
template <typename T>
concept GridVectorElement =
    std::semiregular<std::remove_const_t<T>> && !std::is_reference_v<T>;

// For Q>1: Vector Fields
template <GridVectorElement T, std::size_t Q>
  requires(Q > 0)
struct GridVectorViewMeta {
  using type = mdspan<T, extents<std::size_t, dynamic_extent, Q>, layout_left>;
};

// For Q=1: Scalar Fields
template <GridVectorElement T> struct GridVectorViewMeta<T, 1> {
  using type = mdspan<T, extents<std::size_t, dynamic_extent>>;
};

// Assign proper Type depending on Q
template <GridVectorElement T, std::size_t Q>
  requires(Q > 0)
using GridVectorView = typename GridVectorViewMeta<T, Q>::type;

template <std::semiregular T, std::size_t Q>
  requires(Q > 0)
class GridVector {
private:
  std::vector<T> data_;
  std::size_t num_cells_ = 0;

public:
  using value_type = T;
  static constexpr std::size_t NUM_COMPONENTS = Q;

  explicit GridVector(const CellGrid &cells)
      : data_(cells.size() * NUM_COMPONENTS), num_cells_(cells.size()) {}
  // std::vector automatically manages memory
  ~GridVector() = default;

  GridVector(const GridVector &other)
      : data_(other.data_), num_cells_(other.num_cells_) {}

  GridVector &operator=(const GridVector &other) {
    if (this != &other) {
      data_ = other.data_;
      num_cells_ = other.num_cells_;
    }
    return *this;
  }

  GridVector(GridVector &&other) noexcept = default;

  GridVector &operator=(GridVector &&other) noexcept = default;

  [[nodiscard]] std::size_t allocSize() const noexcept {
    return num_cells_ * NUM_COMPONENTS;
  }

  // vector::data() returns raw pointer
  // The contract is:
  // 1. The caller must not delete the returned pointer
  // 2. The pointer is valid only while the GridVector exists and hasn't been
  // moved from
  [[nodiscard]] T *data() noexcept { return data_.data(); }

  [[nodiscard]] const T *data() const noexcept { return data_.data(); }

  [[nodiscard]] GridVectorView<T, NUM_COMPONENTS> view() {
    return GridVectorView<T, NUM_COMPONENTS>(data_.data(), num_cells_);
  }

  [[nodiscard]] GridVectorView<const T, NUM_COMPONENTS> view() const {
    return GridVectorView<const T, NUM_COMPONENTS>(data_.data(), num_cells_);
  }

  [[nodiscard]] GridVectorView<const T, NUM_COMPONENTS> const_view() const {
    return GridVectorView<const T, NUM_COMPONENTS>(data_.data(), num_cells_);
  }
  operator GridVectorView<T, NUM_COMPONENTS>() { return view(); }
  operator GridVectorView<const T, NUM_COMPONENTS>() const {
    return const_view();
  }

  // element access
    T& operator[](std::size_t i, std::size_t q) {
        if constexpr (Q == 1) {
            return view()[i]; // 1D
        } else {
            return view()[i, q]; // 2D
        }
    }

    const T& operator[](std::size_t i, std::size_t q) const {
        if constexpr (Q == 1) {
            return const_view()[i]; // 1D
        } else {
            return const_view()[i, q]; // 2D
        }
    }

};

} // namespace oktal
