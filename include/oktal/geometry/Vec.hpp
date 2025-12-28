#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace oktal {

template <typename T, std::size_t DIM> class Vec {
private:
  std::array<T, DIM> v_;

public:
  // Constructors
  Vec();
  Vec(const T &value);
  Vec(std::initializer_list<T> list);

  // Converting Constructor
  template <typename S> Vec(const Vec<S, DIM> &other);

  // Element access
  T &operator[](std::size_t i);             // Write
  const T &operator[](std::size_t i) const; // Read

  // Size
  [[nodiscard]] std::size_t size() const;

  // Data access
  T *data();
  [[nodiscard]] const T *data() const;

  // Equality Operator overload
  friend bool operator==(const Vec &lhs, const Vec &rhs) {
    return lhs.v_ == rhs.v_;
  };
  friend bool operator!=(const Vec &lhs, const Vec &rhs) {
    return !(lhs == rhs);
  }

  friend bool operator<(const Vec &lhs, const Vec &rhs) {
    for (std::size_t i = 0; i < DIM; ++i) {
      if (lhs[i] < rhs[i]) {
        return true;
      }
      if (lhs[i] > rhs[i]) {
        return false;
      }
    }
    return false;
  }

  // Range Interface
  T *begin();
  [[nodiscard]] const T *begin() const;
  T *end();
  [[nodiscard]] const T *end() const;

  // Vector Space Operations
  Vec operator+(const Vec &rhs) const;
  Vec operator-() const;
  Vec operator-(const Vec &rhs) const;
  Vec &operator+=(const Vec &rhs);
  Vec &operator-=(const Vec &rhs);
  Vec &operator*=(const T &scalar);
  Vec &operator/=(const T &scalar);

  // Euclidean magnitude
  [[nodiscard]] T sqrMagnitude() const;
  [[nodiscard]] T magnitude() const;

  // Dot product
  [[nodiscard]] T dot(const Vec &other) const;

  // Normalized copy
  [[nodiscard]] Vec normalized() const;

  // Normalizes  in-place
  Vec &normalize();
};

template <typename T, std::size_t DIM, typename S>
Vec<T, DIM> operator*(const S &scalar, const Vec<T, DIM> &vec);

template <typename T, std::size_t DIM, typename S>
Vec<T, DIM> operator*(const Vec<T, DIM> &vec, const S &scalar);

template <typename T, std::size_t DIM, typename S>
Vec<T, DIM> operator/(const Vec<T, DIM> &vec, const S &scalar);

// ====== Implementations ======

// Constructors
template <typename T, std::size_t DIM> Vec<T, DIM>::Vec() {
  std::fill(v_.begin(), v_.end(), T{});
}

template <typename T, std::size_t DIM> Vec<T, DIM>::Vec(const T &value) {
  std::fill(v_.begin(), v_.end(), value);
}

template <typename T, std::size_t DIM>
Vec<T, DIM>::Vec(std::initializer_list<T> list) {
  const std::size_t copy_count = std::min(DIM, list.size());
  std::copy(list.begin(), list.begin() + copy_count, v_.begin());
  if (copy_count < DIM) {
    std::fill(v_.begin() + copy_count, v_.end(), T{});
  }
}

// Converting Constructor
template <typename T, std::size_t DIM>
template <typename S>
Vec<T, DIM>::Vec(const Vec<S, DIM> &other) {
  std::transform(other.begin(), other.end(), v_.begin(),
                 [](const S &val) { return static_cast<T>(val); });
}

// Element access
template <typename T, std::size_t DIM>
T &Vec<T, DIM>::operator[](std::size_t i) {
  // .at() has bound checking with exception, but worse in performance.
  return v_.at(i);
}

template <typename T, std::size_t DIM>
const T &Vec<T, DIM>::operator[](std::size_t i) const {
  return v_.at(i);
}

// Size
template <typename T, std::size_t DIM> std::size_t Vec<T, DIM>::size() const {
  return DIM;
}

// Data access
template <typename T, std::size_t DIM> T *Vec<T, DIM>::data() {
  return v_.data();
}

template <typename T, std::size_t DIM> const T *Vec<T, DIM>::data() const {
  return v_.data();
}

// Range Interface
template <typename T, std::size_t DIM> T *Vec<T, DIM>::begin() {
  return v_.begin();
}
template <typename T, std::size_t DIM> const T *Vec<T, DIM>::begin() const {
  return v_.begin();
}

template <typename T, std::size_t DIM> T *Vec<T, DIM>::end() {
  return v_.end();
}
template <typename T, std::size_t DIM> const T *Vec<T, DIM>::end() const {
  return v_.end();
}

// Vector Space Operations
template <typename T, std::size_t DIM>
Vec<T, DIM> Vec<T, DIM>::operator+(const Vec &rhs) const {
  Vec<T, DIM> result = *this;
  result += rhs;
  return result;
}

template <typename T, std::size_t DIM>
Vec<T, DIM> Vec<T, DIM>::operator-() const {
  Vec<T, DIM> result;
  for (std::size_t i = 0; i < DIM; ++i) {
    result[i] = -(*this)[i];
  }
  return result;
}

template <typename T, std::size_t DIM>
Vec<T, DIM> Vec<T, DIM>::operator-(const Vec &rhs) const {
  Vec<T, DIM> result = *this;
  result -= rhs;
  return result;
}

template <typename T, std::size_t DIM>
Vec<T, DIM> &Vec<T, DIM>::operator+=(const Vec &rhs) {
  for (std::size_t i = 0; i < DIM; ++i) {
    (*this)[i] += rhs[i];
  }
  return *this;
}

template <typename T, std::size_t DIM>
Vec<T, DIM> &Vec<T, DIM>::operator-=(const Vec &rhs) {
  for (std::size_t i = 0; i < DIM; ++i) {
    (*this)[i] -= rhs[i];
  }
  return *this;
}

template <typename T, std::size_t DIM>
Vec<T, DIM> &Vec<T, DIM>::operator*=(const T &scalar) {
  for (std::size_t i = 0; i < DIM; ++i) {
    (*this)[i] *= scalar;
  }
  return *this;
}

template <typename T, std::size_t DIM>
Vec<T, DIM> &Vec<T, DIM>::operator/=(const T &scalar) {
  if (scalar == T{}) { // Using T{} for generic zero check
    throw std::runtime_error("Division by zero error");
  }

  for (std::size_t i = 0; i < DIM; ++i) {
    (*this)[i] /= scalar;
  }
  return *this;
}

// Euclidean magnitude
template <typename T, std::size_t DIM> T Vec<T, DIM>::sqrMagnitude() const {
  T sum_of_squares = T{};
  for (std::size_t i = 0; i < DIM; ++i) {
    sum_of_squares += (*this)[i] * (*this)[i];
  }
  return sum_of_squares;
}

template <typename T, std::size_t DIM> T Vec<T, DIM>::magnitude() const {
  return static_cast<T>(std::sqrt(sqrMagnitude()));
}

// Dot product
template <typename T, std::size_t DIM>
T Vec<T, DIM>::dot(const Vec &other) const {
  T result = T{};
  for (std::size_t i = 0; i < DIM; ++i) {
    result += (*this)[i] * other[i];
  }
  return result;
}

// Normalization
template <typename T, std::size_t DIM>
Vec<T, DIM> Vec<T, DIM>::normalized() const {
  Vec<T, DIM> result = *this;
  result.normalize();
  return result;
}

template <typename T, std::size_t DIM> Vec<T, DIM> &Vec<T, DIM>::normalize() {
  T mag = magnitude();

  if (mag != T{}) {
    (*this) /= mag;
  } else {
    throw std::runtime_error("Magnitude -> Division by zero error");
  }
  return *this;
}

// Non-member Vector Space Operators

// Scalar * Vec (using generic scalar type S)
template <typename T, std::size_t DIM, typename S>
Vec<T, DIM> operator*(const S &scalar, const Vec<T, DIM> &vec) {
  Vec<T, DIM> result = vec;
  // Convert scalar type S to T before calling the internal operator*=
  result *= static_cast<T>(scalar);
  return result;
}

// Vec * Scalar (using generic scalar type S)
template <typename T, std::size_t DIM, typename S>
Vec<T, DIM> operator*(const Vec<T, DIM> &vec, const S &scalar) {
  Vec<T, DIM> result = vec;
  // Convert scalar type S to T before calling the internal operator*=
  result *= static_cast<T>(scalar);
  return result;
}

// Vec / Scalar (using generic scalar type S)
template <typename T, std::size_t DIM, typename S>
Vec<T, DIM> operator/(const Vec<T, DIM> &vec, const S &scalar) {
  Vec<T, DIM> result = vec;
  // Convert scalar type S to T before calling the internal operator/=
  result /= static_cast<T>(scalar);
  return result;
}

using Vec3F = Vec<float, 3>;

using Vec3D = Vec<double, 3>;

} // namespace oktal
