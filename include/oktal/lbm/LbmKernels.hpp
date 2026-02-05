#pragma once

#include "oktal/octree/CellGrid.hpp"
#include "oktal/data/GridVector.hpp"
#include "oktal/lbm/D3Q19.hpp"


namespace oktal::lbm {

// initialization kernels
class InitializePdfs {
public:
    void operator() (CellGrid::CellView cell,
      D3Q19LatticeView pdfs,
      GridVectorView<const double, 1> rho,
      GridVectorView<const double, 3> u) {
        
        const std::size_t cellIdx{cell};
        
        // density
        const double density = rho[cellIdx];

        // velocity
        const double ux = u[cellIdx, 0];
        const double uy = u[cellIdx, 1];
        const double uz = u[cellIdx, 2];

        const double u_sq = ux * ux + uy * uy + uz * uz;

        const auto& CS = D3Q19::CS;
        const auto& W = D3Q19::W;

        for (std::size_t q = 0; q < D3Q19::Q; ++q) {
            const auto cx = static_cast<double>(CS.at(q)[0]);
            const auto cy = static_cast<double>(CS.at(q)[1]);
            const auto cz = static_cast<double>(CS.at(q)[2]);  
            
            const double ci_dot_u = cx * ux + cy * uy + cz * uz;
            const double ci_dot_u_sq = ci_dot_u * ci_dot_u;

            // discrete equilibrium distribution function (eq. 19)
            pdfs[cellIdx, q] = W.at(q) * density * ( 
                1.0 +
                (2.0 * ci_dot_u - u_sq) / (2.0 * D3Q19::C_S_SQ) +
                (ci_dot_u_sq) / (2.0 * D3Q19::C_S_SQ * D3Q19::C_S_SQ)
                );
        }

      }

};

// Macroscopic quantities kernels
class ComputeMacroscopicQuantities {
public:
    void operator() (CellGrid::CellView cell,
      D3Q19LatticeConstView pdfs,
      GridVectorView<double, 1> rho,
      GridVectorView<double, 3> u) {

        const std::size_t cellIdx{cell};
        
        // rho(x,t) (eq. 17)
        double density = 0.0;
        for (std::size_t q = 0; q < D3Q19::Q; ++q) {
            density += pdfs[cellIdx, q];
        }
        rho[cellIdx] = density;

        // u(x,t) (eq. 18)
        double ux = 0.0;
        double uy = 0.0;
        double uz = 0.0;
        const auto& CS = D3Q19::CS;

        for (std::size_t q = 0; q < D3Q19::Q; ++q) {
            const double f_i = pdfs[cellIdx, q];
            ux += f_i * static_cast<double>(CS.at(q)[0]);
            uy += f_i * static_cast<double>(CS.at(q)[1]);
            uz += f_i * static_cast<double>(CS.at(q)[2]);    
        }

        u[cellIdx, 0] = ux / density;
        u[cellIdx, 1] = uy / density;
        u[cellIdx, 2] = uz / density;
    }
};

// Collision kernels
class Collide {
public:
    // constructor
    Collide(double omega) : omega_(omega) {}

    void operator() (CellGrid::CellView cell,
      D3Q19LatticeView pdfs,
      GridVectorView<const double, 1> rho,
      GridVectorView<const double, 3> u) const {
        const std::size_t cellIdx{cell};
        
        const double density = rho[cellIdx];
        const double ux = u[cellIdx, 0];
        const double uy = u[cellIdx, 1];
        const double uz = u[cellIdx, 2];

        const double u_sq = ux * ux + uy * uy + uz * uz;

        const auto& CS = D3Q19::CS;
        const auto& W = D3Q19::W;

        for (std::size_t q = 0; q < D3Q19::Q; ++q) {
            const auto cx = static_cast<double>(CS.at(q)[0]);
            const auto cy = static_cast<double>(CS.at(q)[1]);
            const auto cz = static_cast<double>(CS.at(q)[2]);  
            
            const double ci_dot_u = cx * ux + cy * uy + cz * uz;
            const double ci_dot_u_sq = ci_dot_u * ci_dot_u;

            // discrete equilibrium distribution function (eq. 19)
            const double f_eq = W.at(q) * density * ( 
                1.0 +
                (2.0 * ci_dot_u - u_sq) / (2.0 * D3Q19::C_S_SQ) +
                (ci_dot_u_sq) / (2.0 * D3Q19::C_S_SQ * D3Q19::C_S_SQ)
                );

            const double f_i = pdfs[cellIdx, q];

            // collision (eq. 21)
            pdfs[cellIdx, q] = f_i + omega_ * (f_eq - f_i);
        }
      }

private:
      double omega_;
};

// streaming kernels
class Stream {
public:
    void operator() (CellGrid::CellView cell,
      D3Q19LatticeView pdfsDst,
      D3Q19LatticeConstView pdfsSrc) {
        const std::size_t cellIdx{cell};

        const auto& CS_NO_CENTER = D3Q19::CS_NO_CENTER;

        // stream center direction
        pdfsDst[cellIdx, 0] = pdfsSrc[cellIdx, 0];

        // stream all directions except the center direction
        for (std::size_t q = 1; q < D3Q19::Q;  ++q) {
            const double f_i_star = pdfsSrc[cellIdx, q];

            // find neighbor cell in direction of CS[q]
            const auto neighbor_cell = cell.neighbor(CS_NO_CENTER.at(q - 1));
            const std::size_t neighbor_idx{neighbor_cell};

            // stream to neighbor cell
            pdfsDst[neighbor_idx, q] = f_i_star;
        }
      }
};

} // namespace oktal::lbm