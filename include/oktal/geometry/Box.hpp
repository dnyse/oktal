#pragma once

#include "oktal/geometry/Vec.hpp"

namespace oktal {

template<typename T = double> 
class Box {
public:
    using vector_type = Vec<T, 3>;

    Box() : minCorner_{0, 0, 0}, maxCorner_{0, 0, 0} {}

    Box(const vector_type& minCorner, const vector_type& maxCorner)
        : minCorner_(minCorner), maxCorner_(maxCorner) {}

    vector_type& minCorner() {
        return minCorner_;
    }

    [[nodiscard]] const vector_type& minCorner() const {
        return minCorner_;
    }

    vector_type& maxCorner() {
        return maxCorner_;
    }

    [[nodiscard]] const vector_type& maxCorner() const {
        return maxCorner_;
    }

    [[nodiscard]] vector_type center() const {
        return (minCorner_ + maxCorner_) / static_cast<T>(2);
    }

    [[nodiscard]] vector_type extents() const {
        return maxCorner_ - minCorner_;
    }

    [[nodiscard]] T volume() const {
        vector_type sidelengths = extents();
        return sidelengths[0] * sidelengths[1] * sidelengths[2];
    }

private:
    vector_type minCorner_;
    vector_type maxCorner_;
};

} // namespace oktal
