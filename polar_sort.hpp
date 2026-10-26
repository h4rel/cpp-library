#pragma once

#include <cassert>
#include <type_traits>

namespace geo {

using i128 = __int128_t;

template <class T = long long>
struct Point {
    T x, y;
};

enum class Rotation {
    CCW,  // counterclockwise
    CW    // clockwise
};

enum class SameAngle {
    NearFirst,   // shorter vector first
    FarFirst,    // longer vector first
    Equivalent   // vectors of the same direction compare equivalent
};

enum class ZeroPosition {
    First,
    Last
};

// Sort points by polar angle around `origin`, starting from `ray`.
//
// Requirements:
//   * T must be an integer type (normally int or long long).
//   * `ray` must be nonzero.
//   * Differences, cross products, dot products, and squared lengths
//     must fit in signed __int128. For typical contest constraints,
//     this is sufficient.
//
// Defaults: start at +x, counterclockwise, shorter first for same
// direction, and zero vectors (points equal to origin) first.
//
// Usage:
//   using P = geo::Point<long long>;
//   std::sort(v.begin(), v.end(), geo::PolarLess<long long>{});
//   std::sort(v.begin(), v.end(),
//             geo::PolarLess<long long>{{3, 5}, {0, 1}, geo::Rotation::CW});
template <class T = long long>
struct PolarLess {
    static_assert(std::is_integral_v<T> && !std::is_same_v<T, bool>,
                  "PolarLess requires an integer coordinate type");

    using P = Point<T>;

    P origin, ray;
    Rotation rotation;
    SameAngle same_angle;
    ZeroPosition zero_position;

    PolarLess(
        P origin = {0, 0},
        P ray = {1, 0},
        Rotation rotation = Rotation::CCW,
        SameAngle same_angle = SameAngle::NearFirst,
        ZeroPosition zero_position = ZeroPosition::First
    ) : origin(origin), ray(ray), rotation(rotation),
        same_angle(same_angle), zero_position(zero_position) {
        assert(ray.x != 0 || ray.y != 0);
    }

private:
    struct Vec {
        i128 x, y;
    };

    static i128 cross(Vec a, Vec b) {
        return a.x * b.y - a.y * b.x;
    }

    static i128 dot(Vec a, Vec b) {
        return a.x * b.x + a.y * b.y;
    }

    static i128 norm2(Vec a) {
        return dot(a, a);
    }

    static bool is_zero(Vec a) {
        return a.x == 0 && a.y == 0;
    }

    Vec relative(P p) const {
        return {i128(p.x) - origin.x, i128(p.y) - origin.y};
    }

    // First half: angles in [0, pi); second half: [pi, 2*pi).
    int half(Vec v) const {
        Vec r = {i128(ray.x), i128(ray.y)};
        i128 c = cross(r, v);
        i128 d = dot(r, v);
        bool forward = (rotation == Rotation::CCW) ? (c > 0) : (c < 0);
        return (forward || (c == 0 && d > 0)) ? 0 : 1;
    }

public:
    bool operator()(const P& a, const P& b) const {
        Vec u = relative(a), v = relative(b);
        bool zu = is_zero(u), zv = is_zero(v);
        if (zu || zv) {
            if (zu && zv) return false;
            return zero_position == ZeroPosition::First ? zu : !zu;
        }

        int hu = half(u), hv = half(v);
        if (hu != hv) return hu < hv;

        i128 c = cross(u, v);
        if (c != 0) {
            return rotation == Rotation::CCW ? c > 0 : c < 0;
        }

        if (same_angle == SameAngle::Equivalent) return false;
        if (same_angle == SameAngle::NearFirst) return norm2(u) < norm2(v);
        return norm2(u) > norm2(v);
    }
};

}  // namespace geo
