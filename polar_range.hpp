#pragma once

#include <algorithm>
#include <cassert>
#include <type_traits>
#include <utility>
#include <vector>

#include "polar_sort.hpp"

namespace geo {

// Static polar-angle range counting.
//
// Boundary vectors a and b are DIRECTIONS relative to origin, not absolute
// coordinates. The interval runs from a to b in `rotation` direction.
// The reference ray is the angle-zero direction (default: positive x-axis).
//
// A point equal to origin has no angle and is ignored. Query directions and
// the reference ray must be nonzero. All intermediate integer arithmetic
// (notably cross products) must fit in signed __int128.
//
// Equal-direction endpoints describe a zero-width interval:
//   count(a, a)           = 0   ([a, a))
//   count_closed(a, a)    = number of points in direction a
//   count_open(a, a)      = 0
//   count_open_closed(a,a)= 0
// For a full circle, use size().
//
// Construction: O(N log N). Each count query: O(log N). Space: O(N).
template <class T = long long>
class PolarRange {
    static_assert(std::is_integral_v<T> && !std::is_same_v<T, bool>,
                  "PolarRange requires an integer coordinate type");

public:
    using P = Point<T>;

    explicit PolarRange(
        const std::vector<P>& points,
        P origin = {0, 0},
        P ray = {1, 0},
        Rotation rotation = Rotation::CCW
    ) : rotation_(rotation), ray_{i128(ray.x), i128(ray.y)} {
        assert(!is_zero(ray_) && "The reference ray must be nonzero");
        dirs_.reserve(points.size());
        for (const auto& p : points) {
            Vec v = {i128(p.x) - origin.x, i128(p.y) - origin.y};
            if (is_zero(v)) {
                ++ignored_zero_count_;
            } else {
                dirs_.push_back(v);
            }
        }
        std::sort(dirs_.begin(), dirs_.end(),
                  [this](Vec u, Vec v) { return less(u, v); });
    }

    // Number of nonzero points in the collection.
    int size() const { return static_cast<int>(dirs_.size()); }
    bool empty() const { return dirs_.empty(); }
    int ignored_zero_count() const { return ignored_zero_count_; }

    // [a, b): start included, end excluded (default).
    int count(P a, P b) const {
        return count_impl(a, b, true, false);
    }

    // [a, b]: both ends included.
    int count_closed(P a, P b) const {
        return count_impl(a, b, true, true);
    }

    // (a, b): both ends excluded.
    int count_open(P a, P b) const {
        return count_impl(a, b, false, false);
    }

    // (a, b]: start excluded, end included.
    int count_open_closed(P a, P b) const {
        return count_impl(a, b, false, true);
    }

    // Convenience: choose endpoint inclusion flags explicitly.
    int count(P a, P b, bool include_start, bool include_end) const {
        return count_impl(a, b, include_start, include_end);
    }

private:
    struct Vec {
        i128 x, y;
    };

    Rotation rotation_;
    Vec ray_;
    std::vector<Vec> dirs_;
    int ignored_zero_count_ = 0;

    static bool is_zero(Vec v) { return v.x == 0 && v.y == 0; }
    static i128 cross(Vec u, Vec v) { return u.x * v.y - u.y * v.x; }
    static i128 dot(Vec u, Vec v) { return u.x * v.x + u.y * v.y; }

    // Divide the circle into [0, pi) and [pi, 2pi).
    int half(Vec v) const {
        i128 c = cross(ray_, v);
        i128 d = dot(ray_, v);
        bool forward = rotation_ == Rotation::CCW ? c > 0 : c < 0;
        return (forward || (c == 0 && d > 0)) ? 0 : 1;
    }

    // Strict angle ordering; same-direction vectors compare equivalent.
    bool less(Vec u, Vec v) const {
        int hu = half(u), hv = half(v);
        if (hu != hv) return hu < hv;
        i128 c = cross(u, v);
        return rotation_ == Rotation::CCW ? c > 0 : c < 0;
    }

    int lower(Vec v) const {
        return static_cast<int>(std::lower_bound(
            dirs_.begin(), dirs_.end(), v,
            [this](Vec u, Vec q) { return less(u, q); }
        ) - dirs_.begin());
    }

    int upper(Vec v) const {
        return static_cast<int>(std::upper_bound(
            dirs_.begin(), dirs_.end(), v,
            [this](Vec q, Vec u) { return less(q, u); }
        ) - dirs_.begin());
    }

    int count_impl(P a, P b, bool include_start, bool include_end) const {
        Vec u = {i128(a.x), i128(a.y)};
        Vec v = {i128(b.x), i128(b.y)};
        assert(!is_zero(u) && !is_zero(v) &&
               "Query boundary directions must be nonzero");

        bool u_before_v = less(u, v);
        bool v_before_u = less(v, u);
        if (!u_before_v && !v_before_u) {
            // The two boundary rays coincide: no full-turn interpretation.
            return (include_start && include_end) ? upper(u) - lower(u) : 0;
        }

        int l = include_start ? lower(u) : upper(u);
        int r = include_end ? upper(v) : lower(v);

        if (u_before_v) return r - l;
        // The interval crosses the angle-zero direction.
        return size() - l + r;
    }
};

}  // namespace geo
