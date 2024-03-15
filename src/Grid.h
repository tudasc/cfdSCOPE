#ifndef GRID_H
#define GRID_H

#include "Vector.h"

// "indexing helper"
template <typename T>
class Grid {
    size_t width;
    size_t height;
    size_t depth;

    T resolution;

    // Vector<T> p;
    // Vector<T> v;

  public:
    Grid(size_t width, size_t height, T resolution)
        : width(width), height(height), depth(depth), resolution(resolution) {}

    size_t getWidth() const { return width; }

    size_t getHeight() const { return height; }

    size_t getDepth() const { return depth; }

    size_t cellIndex(size_t x, size_t y, size_t z) const {
        return z * width * height + y * height + x;
    }

    size_t getCellCount() const { return width * height * depth; }
};

template <typename T>
class PressureField {
    Grid<T>& grid;
    Vector<T> field;

  public:
    PressureField(Grid<T>& grid) : grid(grid), field(grid.getCellCount()) {}

    T& getPressure(size_t x, size_t y, size_t z) {
        return field[grid.cellIndex(x, y, z)];
    }

    Vec3<T> div(size_t x, size_t y) const {}
};

template <typename T>
class VelocityField {
    Grid<T>& grid;
    Vector<T> field;

  public:
    VelocityField(Grid<T>& grid) : grid(grid), field(grid.getCellCount()) {}

    /**
     *
     */
    T& getLeftU(size_t x, size_t y, size_t z) {
        return field[grid.cellIndex(x, y, z) * 3];
    }

    T& getRightU(size_t x, size_t y, size_t z) {
        assert(x + 1 < grid.getWidth() && "x out of bounds");
        return field[grid.cellIndex(x + 1, y, z) * 3];
    }

    T& getTopV(size_t x, size_t y, size_t z) {
        return field[grid.cellIndex(x, y, z) * 3];
    }

    T& getBottomV(size_t x, size_t y, size_t z) {
        assert(y + 1 < grid.getHeight() && "y out of bounds");
        return field[grid.cellIndex(x, y + 1, z) * 3];
    }

    T& getFrontW(size_t x, size_t y, size_t z) {
        return field[grid.cellIndex(x, y, z) * 3];
    }

    T& getBackW(size_t x, size_t y, size_t z) {
        assert(z + 1 < grid.getDepth() && "z out of bounds");
        return field[grid.cellIndex(x, y, z + 1) * 3];
    }

    /**
     * Divergence operator using central difference.
     */
    T div(size_t x, size_t y, size_t z) const {

        return getRightU(x, y, z) - getLeftU(x, y, z) + getBottomV(x, y, z) -
               getTopV(x, y, z) + getBackW(x, y, z) - getFrontW(x, y, z);
    }
};

#endif