#ifndef GRID_H
#define GRID_H

#include "Vector.h"
#include <memory>

/**
 * TODO: Implement.
 * For Cavity problem: Dirichlet boundary condiditions for velocity, Neumann for
 * pressure.
 */
class BoundaryCondition {
  public:
};

enum BCType { DIRICHLET, NEUMANN };

enum Boundaries {
    LEFT,
    RIGHT,
    TOP,
    DOWN,
    FRONT,
    BACK,
    NUM_BOUNDARIES = BACK + 1
};

// "indexing helper"
template <typename T>
class Grid {
    size_t width;
    size_t height;
    size_t depth;

    T cellSize;

    // Vector<T> p;
    // Vector<T> v;

  public:
    Grid(size_t width, size_t height, size_t depth, T cellSize)
        : width(width), height(height), depth(depth), cellSize(cellSize) {}

    size_t getWidth() const { return width; }

    size_t getHeight() const { return height; }

    size_t getDepth() const { return depth; }

    size_t cellIndex(size_t x, size_t y, size_t z) const {
        return z * width * height + y * height + x;
    }

    size_t getCellCount() const { return width * height * depth; }

    T getCellSize() const { return cellSize; }

    bool inBounds(size_t i, size_t j, size_t k) const {
        return i >= 0 && i < getWidth() && j >= 0 && j < getHeight() &&
               k >= 0 && k < getDepth();
    }

    Vec3<T> getNearestInsidePos(Vec3<T> pos) {
        auto maxX = width * cellSize;
        auto maxY = height * cellSize;
        auto maxZ = depth * cellSize;
        if (pos.x < 0) {
            pos.x = 0;
        }
        if (pos.x > maxX) {
            pos.x = maxX;
        }

        if (pos.y < 0) {
            pos.y = 0;
        }
        if (pos.y > maxY) {
            pos.y = maxY;
        }

        if (pos.z < 0) {
            pos.z = 0;
        }
        if (pos.z > maxZ) {
            pos.z = maxZ;
        }
        return pos;
    }
};

// template <typename T>
// class Field {
//     std::shared_ptr<Grid<T>> grid;
//     Vector<T> field;
// public:
//     Field(std::shared_ptr<Grid<T>> grid)
//         : grid(grid), field(grid->getCellCount()) {}

//     const T& getPressure(size_t x, size_t y, size_t z) const {
//         return field[grid->cellIndex(x, y, z)];
//     }

//     void setPressure(size_t x, size_t y, size_t z, const T& val) {
//         field[grid->cellIndex(x, y, z)] = val;
//     }

//     size_t getWidth() const { return grid->getWidth(); }

//     size_t getHeight() const { return grid->getHeight(); }

//     size_t getDepth() const { return grid->getDepth(); }

//     T getCellSize() const { return grid->getCellSize(); }

//     Vec3<T> div(size_t x, size_t y) const {}
// }

template <typename T>
class PressureField {
    std::shared_ptr<Grid<T>> grid;
    Vector<T> field;

  public:
    PressureField(std::shared_ptr<Grid<T>> grid)
        : grid(grid), field(grid->getCellCount()) {}

    PressureField(std::shared_ptr<Grid<T>> grid, Vector<T> values)
        : grid(grid), field(std::move(values)) {
        assert(values.getSize() == grid->getCellCount() &&
               "Size does not match");
    }

    T getPressure(size_t x, size_t y, size_t z) const {
        if (!grid->inBounds(x, y, z)) {
            // TODO: BC handling
            return 0;
        }
        return field[grid->cellIndex(x, y, z)];
    }

    void setPressure(size_t x, size_t y, size_t z, const T& val) {
        assert(grid->inBounds(x, y, z) && "Unhandled out of bounds");
        field[grid->cellIndex(x, y, z)] = val;
    }

    size_t getWidth() const { return grid->getWidth(); }

    size_t getHeight() const { return grid->getHeight(); }

    size_t getDepth() const { return grid->getDepth(); }

    T getCellSize() const { return grid->getCellSize(); }

    size_t getNumValues() const { return field.getSize(); }

    Vector<T>& getRawValues() { return field; }

    std::shared_ptr<Grid<T>> getGrid() const { return grid; }

    Vec3<T> div(size_t x, size_t y) const {}

    Vec3<T> laplace(size_t x, size_t y, size_t z) const {
        return (getPressure(x + 1, y, z) + getPressure(x, y + 1, z) +
                getPressure(x, y, z + 1) - 6 * getPressure(x, y, z) +
                getPressure(x - 1, y, z) + getPressure(x, y - 1, z) +
                getPressure(x, y, z - 1)) /
               (grid->getCellSize() * grid->getCellSize());
    }
};

template <typename T>
class VelocityField {
    std::shared_ptr<Grid<T>> grid;
    Vector<T> field;
    BoundaryCondition bcs[Boundaries::NUM_BOUNDARIES];

  public:
    VelocityField(std::shared_ptr<Grid<T>> grid)
        : grid(grid), field(grid->getCellCount() * 3) {}

    VelocityField(std::shared_ptr<Grid<T>> grid, Vector<T> values)
        : grid(grid), field(grid->getCellCount() * 3) {
        if (values.getSize() != grid->getCellCount() * 3) {
            // TODO: Error handling
        }
        field = std::move(values);
    }

    /**
     *
     */
    T getLeftU(size_t x, size_t y, size_t z) const {
        if (!grid->inBounds(x, y, z)) {
            // TODO: BC
            return 0;
        }
        assert(grid->inBounds(x, y, z) && "Unhandled out of bounds");
        return field[grid->cellIndex(x, y, z) * 3];
    }

    T getRightU(size_t x, size_t y, size_t z) const {
        if (!grid->inBounds(x + 1, y, z)) {
            // TODO: BC
            return 0;
        }
        assert(grid->inBounds(x + 1, y, z) && "Unhandled out of bounds");
        return field[grid->cellIndex(x + 1, y, z) * 3];
    }

    T getTopV(size_t x, size_t y, size_t z) const {
        if (!grid->inBounds(x, y, z)) {
            // TODO: BC
            return 0;
        }
        assert(grid->inBounds(x, y, z) && "Unhandled out of bounds");
        return field[grid->cellIndex(x, y, z) * 3 + 1];
    }

    T getBottomV(size_t x, size_t y, size_t z) const {
        if (!grid->inBounds(x, y + 1, z)) {
            // TODO: BC
            return 0;
        }
        assert(grid->inBounds(x, y + 1, z) && "Unhandled out of bounds");
        return field[grid->cellIndex(x, y + 1, z) * 3 + 1];
    }

    T getFrontW(size_t x, size_t y, size_t z) const {
        if (!grid->inBounds(x, y, z)) {
            // TODO: BC
            return 0;
        }
        assert(grid->inBounds(x, y, z) && "Unhandled out of bounds");
        return field[grid->cellIndex(x, y, z) * 3 + 2];
    }

    T getBackW(size_t x, size_t y, size_t z) const {
        if (!grid->inBounds(x, y, z+1)) {
            // TODO: BC
            return 0;
        }
        assert(grid->inBounds(x, y, z + 1) && "Unhandled out of bounds");
        return field[grid->cellIndex(x, y, z + 1) * 3 + 2];
    }

    void setLeftU(size_t x, size_t y, size_t z, const T& value) {
        assert(grid->inBounds(x, y, z) && "Unhandled out of bounds");
        field[grid->cellIndex(x, y, z) * 3] = value;
    }

    void setRightU(size_t x, size_t y, size_t z, const T& value) {
        assert(grid->inBounds(x + 1, y, z) && "Unhandled out of bounds");
        field[grid->cellIndex(x + 1, y, z) * 3] = value;
    }

    void setTopV(size_t x, size_t y, size_t z, const T& value) {
        assert(grid->inBounds(x, y, z) && "Unhandled out of bounds");
        field[grid->cellIndex(x, y, z) * 3 + 1] = value;
    }

    void setBottomV(size_t x, size_t y, size_t z, const T& value) {
        assert(grid->inBounds(x, y + 1, z) && "Unhandled out of bounds");
        field[grid->cellIndex(x, y + 1, z) * 3 + 1] = value;
    }

    void setFrontW(size_t x, size_t y, size_t z, const T& value) {
        assert(grid->inBounds(x, y, z) && "Unhandled out of bounds");
        field[grid->cellIndex(x, y, z) * 3 + 2] = value;
    }

    void setBackW(size_t x, size_t y, size_t z, const T& value) {
        assert(grid->inBounds(x, y, z + 1) && "Unhandled out of bounds");
        field[grid->cellIndex(x, y, z + 1) * 3 + 2] = value;
    }

    size_t getWidth() const { return grid->getWidth(); }

    size_t getHeight() const { return grid->getHeight(); }

    size_t getDepth() const { return grid->getDepth(); }

    T getCellSize() const { return grid->getCellSize(); }

    size_t getNumValues() const { return field.getSize(); }

    Vector<T>& getRawValues() { return field; }

    std::shared_ptr<Grid<T>> getGrid() const { return grid; }

    /**
    * Trilinear interpolation of the velocity field at pos.
    * We assume all wall boundaries. For points outside of the domain, we use the value of the closest point within the domain.
    */
    Vec3<T> trilerp(Vec3<T> pos) const {

        auto cellSize = getCellSize();

        pos = grid->getNearestInsidePos(pos);        

        int iu = (int)(pos.x / cellSize);
        int ju = (int)(pos.y / cellSize - 0.5);
        int ku = (int)(pos.z / cellSize - 0.5);
        float ru = (pos.x / cellSize) - iu;
        float su = (pos.y / cellSize - 0.5) - ju;
        float tu = (pos.z / cellSize - 0.5) - ku;

        T u00 = (1 - ru) * getLeftU(iu, ju, ku);
        +ru* getRightU(iu, ju, ku);
        T u01 = (1 - ru) * getLeftU(iu, ju + 1, ku) +
                ru * getRightU(iu, ju + 1, ku);
        T u10 = (1 - ru) * getLeftU(iu, ju, ku + 1) +
                ru * getRightU(iu, ju, ku + 1);
        T u11 = (1 - ru) * getLeftU(iu, ju + 1, ku + 1) +
                ru * getRightU(iu, ju + 1, ku + 1);

        T u0 = (1 - su) * u00 + su * u01;
        T u1 = (1 - su) * u10 + su * u11;

        T u = (1 - tu) * u0 + tu * u1;

        int iv = (int)(pos.x / cellSize - 0.5);
        int jv = (int)(pos.y / cellSize);
        int kv = (int)(pos.z / cellSize - 0.5);
        float rv = (pos.x / cellSize - 0.5) - iv;
        float sv = (pos.y / cellSize) - jv;
        float tv = (pos.z / cellSize - 0.5) - kv;

        T v00 = (1 - rv) * getTopV(iv, jv, kv);
        +rv* getTopV(iv + 1, jv, kv);
        T v01 = (1 - rv) * getTopV(iv, jv + 1, kv) +
                rv * getTopV(iv + 1, jv + 1, kv);
        T v10 = (1 - rv) * getTopV(iv, jv, kv + 1) +
                rv * getTopV(iv + 1, jv, kv + 1);
        T v11 = (1 - rv) * getTopV(iv, jv + 1, kv + 1) +
                rv * getTopV(iv + 1, jv + 1, kv + 1);

        T v0 = (1 - sv) * v00 + sv * v01;
        T v1 = (1 - sv) * v10 + sv * v11;

        T v = (1 - tv) * v0 + tv * v1;

        int iw = (int)(pos.x / cellSize - 0.5);
        int jw = (int)(pos.y / cellSize - 0.5);
        int kw = (int)(pos.z / cellSize);
        float rw = (pos.x / cellSize - 0.5) - iw;
        float sw = (pos.y / cellSize - 0.5) - jw;
        float tw = (pos.z / cellSize) - kw;

        T w00 = (1 - rw) * getFrontW(iw, jw, kw);
        +rw* getFrontW(iw + 1, jw, kw);
        T w01 = (1 - rw) * getFrontW(iw, jw + 1, kw) +
                rw * getFrontW(iw + 1, jw + 1, kw);
        T w10 = (1 - rw) * getFrontW(iw, jw, kw + 1) +
                rw * getFrontW(iw + 1, jw, kw + 1);
        T w11 = (1 - rw) * getFrontW(iw, jw + 1, kw + 1) +
                rw * getFrontW(iw + 1, jw + 1, kw + 1);

        T w0 = (1 - sw) * w00 + sw * w01;
        T w1 = (1 - sw) * w10 + sw * w11;

        T w = (1 - tw) * w0 + tw * w1;

        return {u, v, w};
    }

    /**
     * Divergence operator using central difference.
     */
    T div(size_t x, size_t y, size_t z) const {
        // TODO: Boundary conditions
        return (getRightU(x, y, z) - getLeftU(x, y, z) + getBottomV(x, y, z) -
                getTopV(x, y, z) + getBackW(x, y, z) - getFrontW(x, y, z)) /
               grid->getCellSize();
    }

    T dudx(size_t x, size_t y, size_t z) const {
        return getLeftU(x + 1, y, z) -
               getLeftU(x - 1, y, z) / 2 * grid->getCellSize();
    }

    T dudy(size_t x, size_t y, size_t z) const {
        return getLeftU(x, y + 1, z) -
               getLeftU(x, y - 1, z) / 2 * grid->getCellSize();
    }

    T dudz(size_t x, size_t y, size_t z) const {
        return getLeftU(x, y, z + 1) -
               getLeftU(x, y, z - 1) / 2 * grid->getCellSize();
    }

    T dvdx(size_t x, size_t y, size_t z) const {
        return getTopV(x + 1, y, z) -
               getTopV(x - 1, y, z) / 2 * grid->getCellSize();
    }

    T dvdy(size_t x, size_t y, size_t z) const {
        return getTopV(x, y + 1, z) -
               getTopV(x, y - 1, z) / 2 * grid->getCellSize();
    }

    T dvdz(size_t x, size_t y, size_t z) const {
        return getTopV(x, y, z + 1) -
               getTopV(x, y, z - 1) / 2 * grid->getCellSize();
    }

    T dwdx(size_t x, size_t y, size_t z) const {
        return getFrontW(x + 1, y, z) -
               getFrontW(x - 1, y, z) / 2 * grid->getCellSize();
    }

    T dwdy(size_t x, size_t y, size_t z) const {
        return getFrontW(x, y + 1, z) -
               getFrontW(x, y - 1, z) / 2 * grid->getCellSize();
    }

    T dwdz(size_t x, size_t y, size_t z) const {
        return getFrontW(x, y, z + 1) -
               getFrontW(x, y, z - 1) / 2 * grid->getCellSize();
    }
};

#endif
