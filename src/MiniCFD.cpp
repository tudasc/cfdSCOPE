
#include "IO.h"
#include "Grid.h"
#include "Matrix.h"
#include "Solver.h"
#include <memory>

using ScalarT = double;

template<typename T>
inline SparseMatrix<T> evalTransportJacobian(const VelocityField<T>& v) {

}

template<typename T>
inline Vector<T> evalTransportEquation(const VelocityField<T>& v) {
    Vector<T> transportEq(v.getNumValues());
    size_t idx = 0;
    for (size_t i = 0; i < v.getWidth(); i++) {
        for (size_t j = 0; j < v.getHeight(); i++) {
            for (size_t k = 0; k < v.getDepth(); i++) {
                auto f_x = - v.getLeftU(i, j, k) * v.dudx(i, j, k) - v.getTopV(i, j, k) * v.dudy(i, j, k) - v.getFrontW(i, j, k) * v.dudz(i,j,k);
                auto f_y = - v.getLeftU(i, j, k) * v.dvdx(i, j, k) - v.getTopV(i, j, k) * v.dvdy(i, j, k) - v.getFrontW(i, j, k) * v.dvdz(i,j,k);
                auto f_z = - v.getLeftU(i, j, k) * v.dwdx(i, j, k) - v.getTopV(i, j, k) * v.dwdy(i, j, k) - v.getFrontW(i, j, k) * v.dwdz(i,j,k);
                transportEq[idx++] = f_x;
                transportEq[idx++] = f_y;
                transportEq[idx++] = f_z;
            }
        }
    }
    return transportEq;
}

template<typename T>
inline VelocityField<T> solveAdvection(const VelocityField<T>& v0, double dt) {
    auto J = evalTransportJacobian(v0);

    // Assemble full matrix for implicit Euler: I - dt*J
    auto A = J;
    
    auto b = evalTransportEquation(v0) * dt;

    auto v_adv = pcg(A, b);

    return {v0.getGrid(), v_adv};
}

template<typename T>
inline void applyPressureCorrection(const VelocityField<T>& v_adv, const PressureField<T>& p) {
    auto size = v_adv.getGrid()->getCellCount();
    Vector<T> v_div(size);
    Vector<T> p_laplace(size);
    size_t idx = 0;
    for (size_t i = 0; i < v_adv.getWidth(); i++) {
        for (size_t j = 0; j < v_adv.getHeight(); i++) {
            for (size_t k = 0; k < v_adv.getDepth(); i++) {
                auto v_div_ijk = v_adv.div(i, j, k);
                auto p_laplace_ijk = p.laplace(i, j, k);
                v_div[idx++] = v_div_ijk;
                p_laplace[idx++] = p_laplace_ijk;
            }
        }
    }
}


int main(int argc, char** argv) {
    // Command line parsing (later)

    // Set up grid
    auto grid = std::make_shared<Grid<ScalarT>>(100, 100, 100, 1);

    // Create velovity and pressure fields
    auto v = std::make_unique<VelocityField<ScalarT>>(grid);
    auto p = std::make_unique<PressureField<ScalarT>>(grid);

    // Init values

    
    double endTime = 10.0;
    double dt = 0.1;

    // Time loop
    for (double t = 0; t < endTime; t+=dt) {
        // - Solve advection: implicit euler
                

        // - Pressure corection
        // - Write field (later)
    }

    // - Validation

    // - Visualization

    return 0;
}
