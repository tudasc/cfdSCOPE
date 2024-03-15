
#include <cstddef>
#include <utility>

template<typename T>
using Vec2 = std::pair<T, T>;

template<typename T>
class Vector {
    size_t size;

public:
    Vector(size_t size) : size(size) {
    }

    T operator[](size_t i) {
        return 0;
    }

    size_t getSize() const {
        return size;
    }
};