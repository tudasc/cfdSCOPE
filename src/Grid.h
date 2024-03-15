#ifndef GRID_H
#define GRID_H

#include "Vector.h"
#include <memory>

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
    Grid(size_t width, size_t height, size_t depth, T resolution)
        : width(width), height(height), depth(depth), resolution(resolution) {}

    size_t getWidth() const { return width; }

    size_t getHeight() const { return height; }

    size_t getDepth() const { return depth; }

    size_t cellIndex(size_t x, size_t y, size_t z) const {
        return z * width * height + y * height + x;
    }

    T getResolution() const { return resolution; }

    size_t getCellCount() const { return width * height * depth; }
};

template <typename T>
class PressureField {
    std::shared_ptr<Grid<T>> grid;
    Vector<T> field;

  public:
    PressureField(std::shared_ptr<Grid<T>>grid) : grid(grid), field(grid->getCellCount()) {}

    const T& getPressure(size_t x, size_t y, size_t z) const {
        return field[grid->cellIndex(x, y, z)];
    }

    void setPressure(size_t x, size_t y, size_t z, const T& val) {
        field[grid->cellIndex(x, y, z)] = val;
    }

    size_t getWidth() const { return grid->getWidth(); }

    size_t getHeight() const { return grid->getHeight(); }

    size_t getDepth() const { return grid->getDepth(); }

    T getResolution() const { return grid->getResolution(); }

    Vec3<T> div(size_t x, size_t y) const {}
};

template <typename T>
class VelocityField {
    std::shared_ptr<Grid<T>> grid;
    Vector<T> field;

  public:
    VelocityField(std::shared_ptr<Grid<T>> grid) : grid(grid), field(grid->getCellCount() * 3) {}

    /**
     *
     */
    const T& getLeftU(size_t x, size_t y, size_t z) const {
        return field[grid->cellIndex(x, y, z) * 3];
    }

    const T& getRightU(size_t x, size_t y, size_t z) const {
        assert(x + 1 < grid.getWidth() && "x out of bounds");
        return field[grid->cellIndex(x + 1, y, z) * 3];
    }

    const T& getTopV(size_t x, size_t y, size_t z) const {
        return field[grid->cellIndex(x, y, z) * 3];
    }

    const T& getBottomV(size_t x, size_t y, size_t z) const {
        assert(y + 1 < grid->getHeight() && "y out of bounds");
        return field[grid->cellIndex(x, y + 1, z) * 3];
    }

    const T& getFrontW(size_t x, size_t y, size_t z) const {
        return field[grid->cellIndex(x, y, z) * 3];
    }

    const T& getBackW(size_t x, size_t y, size_t z) const {
        assert(z + 1 < grid->getDepth() && "z out of bounds");
        return field[grid->cellIndex(x, y, z + 1) * 3];
    }

    void setLeftU(size_t x, size_t y, size_t z, const T& value) {
        field[grid->cellIndex(x, y, z) * 3] = value;
    }

    void setRightU(size_t x, size_t y, size_t z, const T& value) {
        assert(x + 1 < grid->getWidth() && "x out of bounds");
        field[grid->cellIndex(x + 1, y, z) * 3] = value;
    }

    void setTopV(size_t x, size_t y, size_t z, const T& value) {
        field[grid->cellIndex(x, y, z) * 3] = value;
    }

    void setBottomV(size_t x, size_t y, size_t z, const T& value) {
        assert(y + 1 < grid->getHeight() && "y out of bounds");
        field[grid->cellIndex(x, y + 1, z) * 3] = value;
    }

    void setFrontW(size_t x, size_t y, size_t z, const T& value) {
        field[grid->cellIndex(x, y, z) * 3] = value;
    }

    void setBackW(size_t x, size_t y, size_t z, const T& value) {
        assert(z + 1 < grid->getDepth() && "z out of bounds");
        field[grid->cellIndex(x, y, z + 1) * 3] = value;
    }

    size_t getWidth() const { return grid->getWidth(); }

    size_t getHeight() const { return grid->getHeight(); }

    size_t getDepth() const { return grid->getDepth(); }

    T getResolution() const { return grid->getResolution(); }

    /**
     * Divergence operator using central difference.
     */
    T div(size_t x, size_t y, size_t z) const {

        return getRightU(x, y, z) - getLeftU(x, y, z) + getBottomV(x, y, z) -
               getTopV(x, y, z) + getBackW(x, y, z) - getFrontW(x, y, z);
    }
};

#endif