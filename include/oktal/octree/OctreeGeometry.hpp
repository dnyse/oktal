#pragma once

#include "oktal/geometry/Vec.hpp"
#include "oktal/geometry/Box.hpp"
#include "oktal/octree/MortonIndex.hpp"
#include <cmath>

namespace oktal {

class OctreeGeometry {
public:
    OctreeGeometry(Vec3D origin, double sidelength)
        : origin_(origin), sidelength_(sidelength) {}


    OctreeGeometry() : origin_(Vec3D(0.)), sidelength_(1.) {}

    [[nodiscard]] double sidelength() const {
        return sidelength_;
    }

    [[nodiscard]] Vec3D origin() const {
        return origin_;
    }

    [[nodiscard]] double dx(std::size_t level) const {
        return sidelength_ / std::pow(2.0, static_cast<double>(level));
    }

    [[nodiscard]] Vec3D cellExtents(std::size_t level) const {
        const double sidelength = dx(level);
        Vec3D result{sidelength, sidelength, sidelength};
        return result;
    }

    [[nodiscard]] Vec3D cellMinCorner(MortonIndex m) const;

    [[nodiscard]] Vec3D cellMaxCorner(MortonIndex m) const;

    [[nodiscard]] Box<> cellBoundingBox(MortonIndex m) const;
    
    [[nodiscard]] Vec3D cellCenter(MortonIndex m) const;

private:
    Vec3D origin_;
    double sidelength_;
};

} // namespace oktal
