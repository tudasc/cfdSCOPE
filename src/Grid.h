#include "Vector.h"

// "indexing helper"
template <typename T>
class Grid {
    size_t width;
    size_t height;

    T resolution;

    // Vector<T> p;
    // Vector<T> v;

  public:
    Grid(size_t width, size_t height, T resolution)
        : width(width), height(height), resolution(resolution) {}

    size_t getWidth() const { return width; }

    size_t getHeight() const { return height; }

    size_t cellIndex(size_t x, size_t y) const { return 0; }
};

template <typename T>
class VelocityField {
    VelocityField(Grid<T>& grid) {}

    T& getNorthY(size_t x, size_t y) {
        //
        return 0;
    }

    Vec2<T> div(size_t x, size_t y) const {}
};
