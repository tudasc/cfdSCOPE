#ifndef VECTOR_H
#define VECTOR_H

#include <cassert>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <tuple>
#include <utility>
#include <vector>

template <typename T>
using Vec2 = std::pair<T, T>;

template <typename T>
struct Vec3 {
    T x, y, z;

    Vec3<T>(T x, T y, T z) : x(x), y(y), z(z) {}

    Vec3<T>() : Vec3<T>(0, 0, 0) {}

    Vec3<T> operator+(Vec3<T> other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    Vec3<T> operator-() const { return {-x, -y, -z}; }

    Vec3 operator-(Vec3<T> other) const { return *this + (-other); }

    Vec3<T> operator*(T scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }
};

/**
    Vector type
*/
template <typename T>
class Vector {
  public:
    Vector(size_t size) : _data(size, 0.0) {}
    Vector(std::initializer_list<T> l) : _data(l) {}

    Vector(const Vector<T>& other) : _data(other._data) {}

    T& operator[](size_t i) { return _data[i]; }
    const T& operator[](size_t i) const { return _data[i]; }

    Vector<T> operator+(const Vector<T>& other) const {
        assert(this->getSize() == other.getSize() &&
               "Vector addition is only defined for vectors of the same size.");

        Vector<T> res = *this;
        for (int i = 0; i < this->getSize(); i++) {
            res[i] += other[i];
        }
        return res;
    }

    Vector<T> operator-() const {
        Vector<T> res(*this);
        for (int i = 0; i < this->getSize(); i++) {
            res[i] = (*this)[i] * -1.0;
        }
        return res;
    }

    Vector<T> operator-(const Vector<T>& other) const {
        return *this + (-other);
    }

    Vector<T> operator*(T scalar) const {
        Vector<T> res = *this;
        for (int i = 0; i < this->getSize(); i++) {
            res[i] = (*this)[i] * scalar;
        }
        return res;
    }

    inline size_t getSize() const { return _data.size(); }

  private:
    std::vector<T> _data;
};

/**
    Vector-vector dot product
*/
template <typename T>
T dot(const Vector<T>& u, const Vector<T>& v) {
    assert(u.getSize() == v.getSize() &&
           "The dot product requires both operands to be of the same size.");

    T res = 0.0;

    for (size_t i = 0; i < u.getSize(); i++) {
        res += u[i] * v[i];
    }

    return res;
}

/**
    Euclidean vector norm
*/
template <typename T>
T norm(const Vector<T> vec) {
    T res = 0.0;

    for (size_t i = 0; i < vec.getSize(); i++) {
        res += vec[i] * vec[i];
    }

    return std::sqrt(res);
}

#endif
