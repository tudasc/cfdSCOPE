#ifndef GRID_H
#define GRID_H

#include "Vector.h"
#include <memory>

/**
* TODO: Implement.
* For Cavity problem: Dirichlet boundary condiditions for velocity, Neumann for pressure.
*/
class BoundaryCondition {
public:
    enum Type {
        DIRICHLET, NEUMANN
    };

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

    T getCellSize() const {
        return cellSize;
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

    const T& getPressure(size_t x, size_t y, size_t z) const {
        return field[grid->cellIndex(x, y, z)];
    }

    void setPressure(size_t x, size_t y, size_t z, const T& val) {
        field[grid->cellIndex(x, y, z)] = val;
    }

    size_t getWidth() const { return grid->getWidth(); }

    size_t getHeight() const { return grid->getHeight(); }

    size_t getDepth() const { return grid->getDepth(); }

    T getCellSize() const { return grid->getCellSize(); }

    size_t getNumValues() const {
        return field.getSize();
    }

    Vector<T>& getRawValues() {
        return field;
    }

    std::shared_ptr<Grid<T>> getGrid() const {
        return grid;
    }

    Vec3<T> div(size_t x, size_t y) const {}

    Vec3<T> laplace(size_t x, size_t y, size_t z) const {
        return (getPressure(x+1, y, z) + getPressure(x, y+1, z) + getPressure(x, y, z+1) 
        - 6 * getPressure(x, y, z) 
        + getPressure(x-1, y, z) + getPressure(x, y-1, z) + getPressure(x, y, z-1)) /
         (grid->getCellSize() * grid->getCellSize());   
    }
};

template <typename T>
class VelocityField {
    std::shared_ptr<Grid<T>> grid;
    Vector<T> field;

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
    const T& getLeftU(size_t x, size_t y, size_t z) const {
        return field[grid->cellIndex(x, y, z) * 3];
    }

    const T& getRightU(size_t x, size_t y, size_t z) const {
        assert(x + 1 < grid->getWidth() && "x out of bounds");
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

    T getCellSize() const { return grid->getCellSize(); }

    size_t getNumValues() const {
        return field.getSize();
    }

    Vector<T>& getRawValues() {
        return field;
    }

    std::shared_ptr<Grid<T>> getGrid() const {
        return grid;
    }

    /**
     * Divergence operator using central difference.
     */
    T div(size_t x, size_t y, size_t z) const {
        // TODO: Boundary conditions
        return (getRightU(x, y, z) - getLeftU(x, y, z) + getBottomV(x, y, z) -
               getTopV(x, y, z) + getBackW(x, y, z) - getFrontW(x, y, z)) / grid->getCellSize();
    }

    T dudx(size_t x, size_t y, size_t z) const {
        return getLeftU(x+1,y,z) - getLeftU(x-1,y,z) / 2 * grid->getCellSize(); 
    }

    T dudy(size_t x, size_t y, size_t z) const {
        return getLeftU(x,y+1,z) - getLeftU(x,y-1,z) / 2 * grid->getCellSize(); 
    }

    T dudz(size_t x, size_t y, size_t z) const {
        return getLeftU(x,y,z+1) - getLeftU(x,y,z-1) / 2 * grid->getCellSize(); 
    }

    T dvdx(size_t x, size_t y, size_t z) const {
        return getTopV(x+1,y,z) - getTopV(x-1,y,z) / 2 * grid->getCellSize(); 
    }

    T dvdy(size_t x, size_t y, size_t z) const {
        return getTopV(x,y+1,z) - getTopV(x,y-1,z) / 2 * grid->getCellSize(); 
    }

    T dvdz(size_t x, size_t y, size_t z) const {
        return getTopV(x,y,z+1) - getTopV(x,y,z-1) / 2 * grid->getCellSize(); 
    }

    T dwdx(size_t x, size_t y, size_t z) const {
        return getFrontW(x+1,y,z) - getFrontW(x-1,y,z) / 2 * grid->getCellSize(); 
    }

    T dwdy(size_t x, size_t y, size_t z) const {
        return getFrontW(x,y+1,z) - getFrontW(x,y-1,z) / 2 * grid->getCellSize(); 
    }

    T dwdz(size_t x, size_t y, size_t z) const {
        return getFrontW(x,y,z+1) - getFrontW(x,y,z-1) / 2 * grid->getCellSize(); 
    }

};


#endif
