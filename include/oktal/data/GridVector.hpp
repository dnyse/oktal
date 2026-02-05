#pragma once
#include "oktal/octree/CellGrid.hpp"
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <experimental/mdspan>
#include <memory>

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
  std::unique_ptr<T[]>
      data_; // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
  std::size_t num_cells_ = 0;

public:
  using value_type = T;
  static constexpr std::size_t NUM_COMPONENTS = Q;

  explicit GridVector(const CellGrid &cells)
      : data_(new T[cells.size() * NUM_COMPONENTS]), num_cells_(cells.size()) {}
  // std::unique_ptr<T[]> automatically calls delete[] on destruction
  ~GridVector() = default;

  GridVector(const GridVector &other)
      : data_(new T[other.num_cells_ * NUM_COMPONENTS]),
        num_cells_(other.num_cells_) {
    std::copy_n(other.data_.get(), num_cells_ * NUM_COMPONENTS, data_.get());
  }

  GridVector &operator=(const GridVector &other) {
    if (this != &other) {
      if (num_cells_ != other.num_cells_) {
        data_.reset(new T[other.num_cells_ * NUM_COMPONENTS]);
        num_cells_ = other.num_cells_;
      }
      std::copy_n(other.data_.get(), num_cells_ * NUM_COMPONENTS, data_.get());
    }
    return *this;
  }

  GridVector(GridVector &&other) noexcept
      : data_(std::move(other.data_)), num_cells_(other.num_cells_) {
    other.num_cells_ = 0;
  }

  GridVector &operator=(GridVector &&other) noexcept {
    if (this != &other) {
      data_ = std::move(other.data_);
      num_cells_ = other.num_cells_;
      other.num_cells_ = 0;
    }
    return *this;
  }

  [[nodiscard]] std::size_t allocSize() const noexcept {
    return num_cells_ * NUM_COMPONENTS;
  }

  // unique_ptr returns raw poiter
  // The contract is:
  // 1. The caller must not delete the returned pointer
  // 2. The pointer is valid only while the GridVector exists and hasn't been
  // moved from
  [[nodiscard]] T *data() noexcept { return data_.get(); }

  [[nodiscard]] const T *data() const noexcept { return data_.get(); }

  [[nodiscard]] GridVectorView<T, NUM_COMPONENTS> view() {
    return GridVectorView<T, NUM_COMPONENTS>(data_.get(), num_cells_);
  }

  [[nodiscard]] GridVectorView<const T, NUM_COMPONENTS> view() const {
    return GridVectorView<const T, NUM_COMPONENTS>(data_.get(), num_cells_);
  }

  [[nodiscard]] GridVectorView<const T, NUM_COMPONENTS> const_view() const {
    return GridVectorView<const T, NUM_COMPONENTS>(data_.get(), num_cells_);
  }
  operator GridVectorView<T, NUM_COMPONENTS>() { return view(); }
  operator GridVectorView<const T, NUM_COMPONENTS>() const {
    return const_view();
  }

  T &operator[](std::size_t cell_idx)
    requires(NUM_COMPONENTS == 1)
  {
    return view()[cell_idx];
  }

  const T &operator[](std::size_t cell_idx) const
    requires(NUM_COMPONENTS == 1)
  {
    return view()[cell_idx];
  }

  T &operator[](std::size_t cell_idx, std::size_t comp_idx)
    requires(NUM_COMPONENTS > 1)
  {
    return view()[cell_idx, comp_idx];
  }

  const T &operator[](std::size_t cell_idx, std::size_t comp_idx) const
    requires(NUM_COMPONENTS > 1)
  {
    return view()[cell_idx, comp_idx];
  }
};

} // namespace oktal
