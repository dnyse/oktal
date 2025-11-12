#include "oktal/octree/OctreeGeometry.hpp"
#include <stdexcept>
#include <algorithm>

namespace oktal {

    Vec3D OctreeGeometry::cellMinCorner(MortonIndex m) const {
        const Vec3D gridCoords = m.gridCoordinates();
        const double sidelength = dx(m.level());
        return origin_ + ( gridCoords * sidelength );
    }

    Vec3D OctreeGeometry::cellMaxCorner(MortonIndex m) const {
        const Vec3D gridCoords = m.gridCoordinates();
        const double sidelength = dx(m.level());
        return origin_ + ( (gridCoords + Vec3D(1.)) * sidelength );
    }

    Box<> OctreeGeometry::cellBoundingBox(MortonIndex m) const {
        return Box<>{cellMinCorner(m), cellMaxCorner(m)};
    }

    Vec3D OctreeGeometry::cellCenter(MortonIndex m) const {
        return (cellMinCorner(m) + cellMaxCorner(m)) / 2.0;

    }

} // namespace oktal
