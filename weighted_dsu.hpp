#ifndef WEIGHTED_DSU_HPP
#define WEIGHTED_DSU_HPP 1

#include <algorithm>
#include <cassert>
#include <vector>

// T must be an additive abelian group, with T{} as its identity.
// All intermediate arithmetic must be representable in T.
template <class T = long long> struct weighted_dsu {
  public:
    weighted_dsu() : _n(0) {}
    explicit weighted_dsu(int n)
        : _n(n), parent_or_size(n, -1), diff_weight(n, T{}) {}

    // Add potential(b) - potential(a) = w (default: zero).
    // Return the leader, or -1 if the constraint contradicts existing ones.
    // On contradiction, only path compression may change internal storage.
    int merge(int a, int b, T w = T{}) {
        assert(0 <= a && a < _n);
        assert(0 <= b && b < _n);
        int x = leader(a), y = leader(b);
        if (x == y) return diff_weight[b] - diff_weight[a] == w ? x : -1;

        // potential(y) - potential(x)
        w = w + diff_weight[a] - diff_weight[b];
        if (-parent_or_size[x] < -parent_or_size[y]) {
            std::swap(x, y);
            w = -w;
        }
        parent_or_size[x] += parent_or_size[y];
        parent_or_size[y] = x;
        diff_weight[y] = w;
        return x;
    }

    bool same(int a, int b) {
        assert(0 <= a && a < _n);
        assert(0 <= b && b < _n);
        return leader(a) == leader(b);
    }

    int leader(int a) {
        assert(0 <= a && a < _n);
        if (parent_or_size[a] < 0) return a;
        int p = parent_or_size[a];
        int r = leader(p);
        diff_weight[a] = diff_weight[a] + diff_weight[p];
        return parent_or_size[a] = r;
    }

    int size(int a) {
        assert(0 <= a && a < _n);
        return -parent_or_size[leader(a)];
    }

    std::vector<std::vector<int>> groups() {
        std::vector<int> leader_buf(_n), group_size(_n);
        for (int i = 0; i < _n; i++) {
            leader_buf[i] = leader(i);
            group_size[leader_buf[i]]++;
        }
        std::vector<std::vector<int>> result(_n);
        for (int i = 0; i < _n; i++) result[i].reserve(group_size[i]);
        for (int i = 0; i < _n; i++) result[leader_buf[i]].push_back(i);
        result.erase(
            std::remove_if(result.begin(), result.end(),
                           [](const std::vector<int>& v) { return v.empty(); }),
            result.end());
        return result;
    }

    // potential(a) - potential(leader(a)); may change when the leader changes.
    T weight(int a) {
        leader(a);
        return diff_weight[a];
    }

    // potential(b) - potential(a). Requires same(a, b).
    T diff(int a, int b) {
        assert(0 <= a && a < _n);
        assert(0 <= b && b < _n);
        int x = leader(a), y = leader(b);
        assert(x == y);
        (void)x;
        (void)y;
        return diff_weight[b] - diff_weight[a];
    }

  private:
    int _n;
    // Root: negative component size; otherwise: parent vertex.
    std::vector<int> parent_or_size;
    // potential(vertex) - potential(parent); zero for roots.
    std::vector<T> diff_weight;
};

#endif  // WEIGHTED_DSU_HPP
