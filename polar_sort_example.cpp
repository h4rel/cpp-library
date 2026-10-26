#include <algorithm>
#include <iostream>
#include <vector>
#include "polar_sort.hpp"

int main() {
    using P = geo::Point<long long>;
    std::vector<P> points = {
        {0, -1}, {-1, 0}, {1, 1}, {0, 1}, {1, 0}, {2, 0}
    };

    std::sort(points.begin(), points.end(), geo::PolarLess<long long>{});

    for (const auto& p : points) {
        std::cout << p.x << ' ' << p.y << '\n';
    }
}
