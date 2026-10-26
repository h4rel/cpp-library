#include <iostream>
#include <vector>
#include "polar_range.hpp"

int main() {
    using P = geo::Point<long long>;
    std::vector<P> ps = {
        {1, 0}, {2, 0}, {1, 1}, {0, 1}, {-1, 0}, {0, -1},
        {0, 0}  // ignored: the origin has no polar angle
    };

    geo::PolarRange<long long> range(ps);

    std::cout << range.count({1, 0}, {0, 1}) << '\n';
    // 3: [0 degrees, 90 degrees) -> (1,0), (2,0), (1,1)

    std::cout << range.count_closed({1, 0}, {0, 1}) << '\n';
    // 4: [0 degrees, 90 degrees], also includes (0,1)

    std::cout << range.count({0, -1}, {0, 1}) << '\n';
    // 4: [270 degrees, 90 degrees), crossing 0 degrees

    std::cout << range.count_closed({1, 0}, {2, 0}) << '\n';
    // 2: same direction, closed single-angle interval

    std::cout << range.size() << ' ' << range.ignored_zero_count() << '\n';
    // 6 1

    geo::PolarRange<long long> cw(ps, {0, 0}, {0, 1}, geo::Rotation::CW);
    std::cout << cw.count({0, 1}, {1, 0}) << '\n';
    // 2: clockwise [90 degrees, 0 degrees): (0,1), (1,1)
}
