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

    T cellSize;

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

template<typename T, unsigned dims>
class Field {
protected:
    std::shared_ptr<Grid<T>> grid;
    Vector<T> field;

public:
    Field(std::shared_ptr<Grid<T>> grid) 
        : grid(grid), field(getVectorSize()) {}

    Field(std::shared_ptr<Grid<T>> grid, Vector<T> values) 
        : grid(grid), field(std::move(values)) {
        assert(values.getSize() == getVectorSize() &&
               "Size does not match");
    }

    const size_t getVectorSize() const {
        return grid->getCellCount() * dims;
    }

    size_t getWidth() const { return grid->getWidth(); }

    size_t getHeight() const { return grid->getHeight(); }

    size_t getDepth() const { return grid->getDepth(); }

    T getCellSize() const { return grid->getCellSize(); }

    Vector<T>& getRawValues() { return field; }

    const Vector<T>& getRawValues() const { return field; }

    std::shared_ptr<Grid<T>> getGrid() const { return grid; }

    T getValue(size_t x, size_t y, size_t z, unsigned dim = 0) const {
         if (!grid->inBounds(x, y, z)) {
            // Note: Boundary conditions need to be handled explicitly
            return 0;
        }
        return field[grid->cellIndex(x, y, z) * dims + dim];
    }

    void setValue(size_t x, size_t y, size_t z, T val, unsigned dim = 0) {
        assert(grid->inBounds(x, y, z) && "Unhandled out of bounds");
        field[grid->cellIndex(x, y, z) * dims + dim] = val;
    }
    

};


template <typename T>
class PressureField: public Field<T,1> {
public:
    PressureField(std::shared_ptr<Grid<T>> grid)
        : Field<T,1>(grid) {}

    PressureField(std::shared_ptr<Grid<T>> grid, Vector<T> values)
        : Field<T,1>(grid, values) {}

    T getPressure(size_t x, size_t y, size_t z) const {
        return this->getValue(x, y, z);
    }

    void setPressure(size_t x, size_t y, size_t z, T val) {
        this->setValue(x, y, z, val);
    }
    
};

template <typename T>
class VelocityField: public Field<T, 3> {
public:
    VelocityField(std::shared_ptr<Grid<T>> grid)
        : Field<T,3>(grid) {}

    VelocityField(std::shared_ptr<Grid<T>> grid, Vector<T> values)
        : Field<T,3>(grid, values) {}
    
    T getLeftU(size_t x, size_t y, size_t z) const {
        return this->getValue(x, y, z, 0);
    }

    T getRightU(size_t x, size_t y, size_t z) const {
        return this->getValue(x+1, y, z, 0);
    }

    T getTopV(size_t x, size_t y, size_t z) const {
        return this->getValue(x, y, z, 1);
    }

    T getBottomV(size_t x, size_t y, size_t z) const {
        return this->getValue(x, y + 1, z, 1);
    }

    T getFrontW(size_t x, size_t y, size_t z) const {
        return this->getValue(x, y, z, 2);
    }

    T getBackW(size_t x, size_t y, size_t z) const {
       return this->getValue(x, y, z+1, 2);
    }

    void setLeftU(size_t x, size_t y, size_t z, const T& value) {
        this->setValue(x, y, z, value, 0);
    }

    void setRightU(size_t x, size_t y, size_t z, const T& value) {
        this->setValue(x+1, y, z, value, 0);
    }

    void setTopV(size_t x, size_t y, size_t z, const T& value) {
        this->setValue(x, y, z, value, 1);
    }

    void setBottomV(size_t x, size_t y, size_t z, const T& value) {
        this->setValue(x, y+1, z, value, 1);
    }

    void setFrontW(size_t x, size_t y, size_t z, const T& value) {
        this->setValue(x, y, z, value, 2);
    }

    void setBackW(size_t x, size_t y, size_t z, const T& value) {
        this->setValue(x, y, z+1, value, 2);
    }

    /**
     * Trilinear interpolation of the velocity field at pos.
     * We assume all wall boundaries. For points outside of the domain, we use
     * the value of the closest point within the domain.
     */
    Vec3<T> trilerp(Vec3<T> pos) const {

        auto cellSize = this->getCellSize();

        pos = this->grid->getNearestInsidePos(pos);

        int iu = (int)(pos.x / cellSize);
        int ju = (int)(pos.y / cellSize - 0.5);
        int ku = (int)(pos.z / cellSize - 0.5);
        float ru = (pos.x / cellSize) - iu;
        float su = (pos.y / cellSize - 0.5) - ju;
        float tu = (pos.z / cellSize - 0.5) - ku;

        T u00 = (1 - ru) * getLeftU(iu, ju, ku) + ru * getRightU(iu, ju, ku);
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

        T v00 = (1 - rv) * getTopV(iv, jv, kv) + rv * getTopV(iv + 1, jv, kv);
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

        T w00 =
            (1 - rw) * getFrontW(iw, jw, kw) + rw * getFrontW(iw + 1, jw, kw);
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
        // Note: Boundary conditions not handled
        return (getRightU(x, y, z) - getLeftU(x, y, z) + getBottomV(x, y, z) -
                getTopV(x, y, z) + getBackW(x, y, z) - getFrontW(x, y, z)) /
               this->grid->getCellSize();
    }

};


#endif
