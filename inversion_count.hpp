#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

// C++17. An inversion is a pair (i, j) with i < j and a[j] < a[i].
// Equal values are not counted. Input vectors are never modified.
// T must be copyable and operator< must define a strict weak ordering.
// Each public function: O(N log N) time, O(N) additional space.

namespace inversion_count_detail {

class FenwickTree {
public:
    explicit FenwickTree(std::size_t n) : bit_(n + 1, 0) {}

    void add(std::size_t index) {
        for (std::size_t i = index + 1; i < bit_.size(); i += i & -i)
            ++bit_[i];
    }

    // Number of inserted elements whose ranks belong to [0, end).
    long long prefix_sum(std::size_t end) const {
        long long sum = 0;
        for (std::size_t i = end; i > 0; i -= i & -i) sum += bit_[i];
        return sum;
    }

private:
    std::vector<long long> bit_;
};

template <class T, class Callback>
void for_each_contribution(const std::vector<T>& a, Callback&& callback) {
    std::vector<T> values = a;
    std::sort(values.begin(), values.end());
    // Use comparison equivalence so T does not need operator==.
    values.erase(std::unique(values.begin(), values.end(),
                             [](const T& x, const T& y) {
                                 return !(x < y) && !(y < x);
                             }),
                 values.end());

    FenwickTree counts(values.size());
    for (std::size_t j = 0; j < a.size(); ++j) {
        const auto rank = static_cast<std::size_t>(
            std::lower_bound(values.begin(), values.end(), a[j]) - values.begin());
        // There are j earlier elements. Exclude all values <= a[j].
        const long long added = static_cast<long long>(j) - counts.prefix_sum(rank + 1);
        callback(j, added);
        counts.add(rank);
    }
}

}  // namespace inversion_count_detail

// Total inversion count. The empty vector returns 0.
template <class T>
long long inversion_count(const std::vector<T>& a) {
    long long total = 0;
    inversion_count_detail::for_each_contribution(
        a, [&](std::size_t, long long added) { total += added; });
    return total;
}

// result[k] is the inversion count in a[0..k), for 0 <= k <= N.
// Returns N+1 entries, starting with 0. The empty vector returns {0}.
template <class T>
std::vector<long long> inversion_count_prefix(const std::vector<T>& a) {
    std::vector<long long> result(a.size() + 1, 0);
    inversion_count_detail::for_each_contribution(
        a, [&](std::size_t j, long long added) { result[j + 1] = result[j] + added; });
    return result;
}

// result[j] is the number of indices i < j with a[j] < a[i].
// Returns N entries. The empty vector returns {}.
template <class T>
std::vector<long long> inversions_ending_at(const std::vector<T>& a) {
    std::vector<long long> result(a.size(), 0);
    inversion_count_detail::for_each_contribution(
        a, [&](std::size_t j, long long added) { result[j] = added; });
    return result;
}
