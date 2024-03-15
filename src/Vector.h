
#include <cassert>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

template <typename T>
using Vec2 = std::pair<T, T>;

template <typename T>
struct Vec3 {
    T x, y, z;
};

/**
    Vector type
*/
template <typename T>
class Vector {
  public:
    Vector(size_t size) : _data(size, 0.0) {}

    T& operator[](size_t i) { return _data[i]; }

    const T& operator[](size_t i) const { return _data[i]; }

    size_t getSize() const { return _data.size(); }

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
        res += u[i] + v[i];
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
