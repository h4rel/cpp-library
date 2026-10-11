#pragma once

#include <algorithm>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

// C++17 / 0-indexed / unweighted undirected tree.
// build(): O(N log N) time and O(N log N) space (including paths).
// All tree traversals and decomposition steps are iterative.
class CentroidDecomposition {
public:
    struct PathEntry {
        int centroid;
        int distance;  // Distance in the original tree.
    };

    int root = -1;
    std::vector<int> parent;  // Parent in the centroid tree; root has -1.
    std::vector<int> depth;   // Depth in the centroid tree; root has 0.
    // Ancestor centroids, ordered from the root to vertex v itself.
    std::vector<std::vector<PathEntry>> paths;

    explicit CentroidDecomposition(int n) {
        if (n < 0) throw std::invalid_argument("negative vertex count");
        graph_.resize(n);
        parent.assign(n, -1);
        depth.assign(n, -1);
        paths.resize(n);
    }

    // Pass a symmetric adjacency list: each edge u-v appears in both lists.
    // Lvalues are copied; std::move(tree) transfers ownership.
    // This constructor does not build the centroid decomposition.
    explicit CentroidDecomposition(std::vector<std::vector<int>> tree)
        : CentroidDecomposition(static_cast<int>(tree.size())) {
        graph_ = std::move(tree);
        std::size_t degree_sum = 0;
        for (const auto& neighbors : graph_) degree_sum += neighbors.size();
        if (degree_sum % 2 != 0)
            throw std::invalid_argument("adjacency list must be symmetric");
        edge_count_ = degree_sum / 2;
        if (size() > 0) validate_tree(0);
    }

    int size() const { return static_cast<int>(graph_.size()); }

    // All centroids of the ORIGINAL tree, in ascending vertex order.
    // O(N) time and O(N) extra space. Independent of build().
    // Returns one or two vertices for a nonempty tree, {} for an empty tree.
    std::vector<int> centroids() const {
        const int n = size();
        if (n == 0) return {};
        validate_tree(0);

        std::vector<int> dfs_parent(n, -1), subtree_size(n, 1), order{0};
        order.reserve(n);
        for (std::size_t i = 0; i < order.size(); ++i) {
            const int v = order[i];
            for (int to : graph_[v]) {
                if (to == dfs_parent[v]) continue;
                dfs_parent[to] = v;
                order.push_back(to);
            }
        }
        for (auto it = order.rbegin(); it != order.rend(); ++it) {
            const int v = *it;
            if (dfs_parent[v] != -1)
                subtree_size[dfs_parent[v]] += subtree_size[v];
        }

        std::vector<int> result;
        for (int v = 0; v < n; ++v) {
            int largest = n - subtree_size[v];
            for (int to : graph_[v]) {
                if (dfs_parent[to] == v)
                    largest = std::max(largest, subtree_size[to]);
            }
            if (largest <= n / 2) result.push_back(v);
        }
        return result;
    }

    // Add each undirected edge exactly once, before build().
    // After adding an edge, rebuild before using decomposition data.
    void add_edge(int u, int v) {
        check_vertex(u);
        check_vertex(v);
        if (u == v) throw std::invalid_argument("self-loop");
        graph_[u].push_back(v);
        graph_[v].push_back(u);
        ++edge_count_;
    }

    // start is any vertex in the input tree, not necessarily its centroid.
    // Calling build() again is supported. The empty tree has root == -1.
    void build(int start = 0) {
        const int n = size();
        if (n == 0) {
            root = -1;
            return;
        }
        check_vertex(start);
        validate_tree(start);

        root = -1;
        std::fill(parent.begin(), parent.end(), -1);
        std::fill(depth.begin(), depth.end(), -1);
        for (auto& path : paths) path.clear();

        std::vector<char> removed(n, false);
        std::vector<int> dfs_parent(n, -1), subtree_size(n);
        std::vector<int> order;
        order.reserve(n);
        std::vector<std::tuple<int, int, int>> tasks{{start, -1, 0}};
        std::vector<std::tuple<int, int, int>> stack;

        while (!tasks.empty()) {
            const auto [entry, centroid_parent, level] = tasks.back();
            tasks.pop_back();

            // Enumerate the component, ignoring previously removed centroids.
            order.clear();
            order.push_back(entry);
            dfs_parent[entry] = -1;
            for (std::size_t i = 0; i < order.size(); ++i) {
                const int v = order[i];
                for (int to : graph_[v]) {
                    if (removed[to] || to == dfs_parent[v]) continue;
                    dfs_parent[to] = v;
                    order.push_back(to);
                }
            }

            for (int v : order) subtree_size[v] = 1;
            for (auto it = order.rbegin(); it != order.rend(); ++it) {
                const int v = *it;
                if (dfs_parent[v] != -1)
                    subtree_size[dfs_parent[v]] += subtree_size[v];
            }

            const int component_size = static_cast<int>(order.size());
            int centroid = -1;
            for (int v : order) {
                int largest = component_size - subtree_size[v];
                for (int to : graph_[v]) {
                    if (!removed[to] && dfs_parent[to] == v)
                        largest = std::max(largest, subtree_size[to]);
                }
                if (largest <= component_size / 2) {
                    centroid = v;
                    break;
                }
            }

            parent[centroid] = centroid_parent;
            depth[centroid] = level;
            if (centroid_parent == -1) root = centroid;

            // Record distances before removing the new centroid.
            stack.clear();
            stack.emplace_back(centroid, -1, 0);
            while (!stack.empty()) {
                const auto [v, p, distance] = stack.back();
                stack.pop_back();
                paths[v].push_back({centroid, distance});
                for (int to : graph_[v]) {
                    if (!removed[to] && to != p)
                        stack.emplace_back(to, v, distance + 1);
                }
            }

            removed[centroid] = true;
            for (int to : graph_[centroid]) {
                if (!removed[to]) tasks.emplace_back(to, centroid, level + 1);
            }
        }
    }

private:
    std::vector<std::vector<int>> graph_;
    std::size_t edge_count_ = 0;

    void check_vertex(int v) const {
        if (v < 0 || v >= size()) throw std::out_of_range("vertex index");
    }

    void validate_tree(int start) const {
        if (edge_count_ != static_cast<std::size_t>(size() - 1))
            throw std::invalid_argument("expected N-1 edges");
        std::vector<char> seen(size(), false);
        std::vector<int> dfs_parent(size(), -1);
        std::vector<int> order{start};
        seen[start] = true;
        for (std::size_t i = 0; i < order.size(); ++i) {
            const int v = order[i];
            int parent_edges = 0;
            for (int to : graph_[v]) {
                check_vertex(to);
                if (to == dfs_parent[v]) {
                    ++parent_edges;
                    continue;
                }
                if (seen[to])
                    throw std::invalid_argument("cycle or duplicate/asymmetric edge");
                seen[to] = true;
                dfs_parent[to] = v;
                order.push_back(to);
            }
            if (dfs_parent[v] != -1 && parent_edges != 1)
                throw std::invalid_argument("adjacency list must be symmetric");
        }
        if (static_cast<int>(order.size()) != size())
            throw std::invalid_argument("tree must be connected");
    }
};
