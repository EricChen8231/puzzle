#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <thread>
#include <atomic>
#include <functional>
#include <mutex>
#include <algorithm>
#include <array>

using namespace std;

// ----------------- Basic 3D vector / frame types -----------------

struct Vec3 {
    int x, y, z;
    bool operator==(const Vec3 &o) const {
        return x == o.x && y == o.y && z == o.z;
    }
};

struct Vec3Hash {
    size_t operator()(const Vec3 &v) const noexcept {
        return (uint64_t(v.x) * 1315423911u)
             ^ (uint64_t(v.y) * 2654435761u)
             ^ (uint64_t(v.z) * 97531u);
    }
};

struct Frame {
    Vec3 u, v, n; // in-plane axes + normal
    bool operator==(const Frame &o) const {
        return u == o.u && v == o.v && n == o.n;
    }
};

inline Vec3 operator+(const Vec3 &a, const Vec3 &b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline Vec3 operator-(const Vec3 &a, const Vec3 &b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline Vec3 operator-(const Vec3 &a) {
    return {-a.x, -a.y, -a.z};
}
inline int dot(const Vec3 &a, const Vec3 &b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
inline Vec3 cross(const Vec3 &a, const Vec3 &b) {
    return {
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

// rotate v by +90 degrees around 'axis' (right-hand rule)
Vec3 rotate90_plus(const Vec3 &v, const Vec3 &axis) {
    if (dot(v, axis) != 0) return v; // parallel to axis
    return cross(axis, v);
}

// rotate v by -90 degrees around 'axis'
Vec3 rotate90_minus(const Vec3 &v, const Vec3 &axis) {
    if (dot(v, axis) != 0) return v; // parallel to axis
    return cross(v, axis);
}

// ----------------- Net / grid representation -----------------

struct Cell {
    int r, c; // 0-based row/col in the net
};

enum Dir { UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3 };

struct Edge {
    int to;
    Dir dir;
};

struct Marker {
    int r, c;    // 0-based net coordinates
    int val;     // number (not used in constraints)
    char type;   // 'n', 'c', 's'
    int cellIndex; // index in cells, or -1 if not on '#'
};

struct BoxDims {
    int minX, maxX;
    int minY, maxY;
    int minZ, maxZ;
};

// Given a surface area (number of cells) of a rectangular prism, return the
// maximum possible area of any single face. This is used to prune embeddings:
// if any outward normal already has more cells than the largest face could
// hold, the partial embedding cannot ever close into a box.
int max_face_area_from_surface(int surfaceCells) {
    // surfaceCells = 2 * (a*b + a*c + b*c)
    if (surfaceCells <= 0 || surfaceCells % 2 != 0) return surfaceCells;
    int half = surfaceCells / 2;
    int best = 0;

    for (int a = 1; a * a <= half; ++a) {
        for (int b = 1; a * b <= half; ++b) {
            int rhs = half - a * b; // a*c + b*c
            int denom = a + b;
            if (rhs <= 0 || rhs % denom != 0) continue;
            int c = rhs / denom;
            if (c <= 0) continue;
            int face1 = a * b;
            int face2 = a * c;
            int face3 = b * c;
            best = max({best, face1, face2, face3});
        }
    }

    if (best == 0) return surfaceCells; // fallback (should not happen for valid nets)
    return best;
}

// ----------------- Box and marker checks -----------------

bool compute_box_dims_and_check(const vector<Vec3> &pos3d, int numCells, BoxDims &dims) {
    if (numCells == 0) return false;

    bool first = true;
    for (int i = 0; i < numCells; ++i) {
        Vec3 p = pos3d[i];
        if (first) {
            dims.minX = dims.maxX = p.x;
            dims.minY = dims.maxY = p.y;
            dims.minZ = dims.maxZ = p.z;
            first = false;
        } else {
            dims.minX = min(dims.minX, p.x);
            dims.maxX = max(dims.maxX, p.x);
            dims.minY = min(dims.minY, p.y);
            dims.maxY = max(dims.maxY, p.y);
            dims.minZ = min(dims.minZ, p.z);
            dims.maxZ = max(dims.maxZ, p.z);
        }
    }

    auto check_face = [&](auto getterKey, auto isOnPlane) -> bool {
        vector<pair<int,int>> pts;
        pts.reserve(numCells);
        for (int i = 0; i < numCells; ++i) {
            if (!isOnPlane(i)) continue;
            pts.push_back(getterKey(i));
        }
        if (pts.empty()) return false;
        int minA = pts[0].first, maxA = pts[0].first;
        int minB = pts[0].second, maxB = pts[0].second;
        for (auto &p : pts) {
            minA = min(minA, p.first);
            maxA = max(maxA, p.first);
            minB = min(minB, p.second);
            maxB = max(maxB, p.second);
        }
        int expected = (maxA - minA + 1) * (maxB - minB + 1);
        return (int)pts.size() == expected;
    };

    bool ok = true;
    ok &= check_face(
        [&](int i) { return make_pair(pos3d[i].y, pos3d[i].z); },
        [&](int i) { return pos3d[i].x == dims.minX; }
    );
    ok &= check_face(
        [&](int i) { return make_pair(pos3d[i].y, pos3d[i].z); },
        [&](int i) { return pos3d[i].x == dims.maxX; }
    );
    ok &= check_face(
        [&](int i) { return make_pair(pos3d[i].x, pos3d[i].z); },
        [&](int i) { return pos3d[i].y == dims.minY; }
    );
    ok &= check_face(
        [&](int i) { return make_pair(pos3d[i].x, pos3d[i].z); },
        [&](int i) { return pos3d[i].y == dims.maxY; }
    );
    ok &= check_face(
        [&](int i) { return make_pair(pos3d[i].x, pos3d[i].y); },
        [&](int i) { return pos3d[i].z == dims.minZ; }
    );
    ok &= check_face(
        [&](int i) { return make_pair(pos3d[i].x, pos3d[i].y); },
        [&](int i) { return pos3d[i].z == dims.maxZ; }
    );
    if (!ok) return false;

    // Every cell must lie on one of the 6 faces
    for (int i = 0; i < numCells; ++i) {
        Vec3 p = pos3d[i];
        bool onSurface =
            (p.x == dims.minX || p.x == dims.maxX ||
             p.y == dims.minY || p.y == dims.maxY ||
             p.z == dims.minZ || p.z == dims.maxZ);
        if (!onSurface) return false;
    }

    // Check total area matches 2*(ab+bc+ac)
    int dx = dims.maxX - dims.minX + 1;
    int dy = dims.maxY - dims.minY + 1;
    int dz = dims.maxZ - dims.minZ + 1;
    int expectedCells = 2 * (dx*dy + dx*dz + dy*dz);
    if (expectedCells != numCells) return false;

    return true;
}

bool check_markers(const vector<Marker> &markers,
                   const vector<Cell> &cells,
                   const vector<Vec3> &pos3d,
                   const vector<Frame> &frame3d,
                   const BoxDims &dims) {
    int numCells = (int)cells.size();

    // (r,c) -> cell index
    unordered_map<long long,int> rcToIdx;
    rcToIdx.reserve(numCells * 2);
    for (int i = 0; i < numCells; ++i) {
        long long key =
            (static_cast<long long>(cells[i].r) << 32) |
            static_cast<unsigned long long>(cells[i].c);
        rcToIdx[key] = i;
    }

    vector<Marker> muse = markers;
    for (auto &m : muse) {
        long long key =
            (static_cast<long long>(m.r) << 32) |
            static_cast<unsigned long long>(m.c);
        auto it = rcToIdx.find(key);
        if (it == rcToIdx.end()) {
            m.cellIndex = -1;
        } else {
            m.cellIndex = it->second;
        }
    }

    // ---- Circle constraints ----
    vector<int> circleIdx;
    for (int i = 0; i < (int)muse.size(); ++i) {
        if (muse[i].type == 'c' && muse[i].cellIndex != -1)
            circleIdx.push_back(i);
    }
    vector<bool> used(circleIdx.size(), false);

    auto same3 = [](const Vec3 &a, const Vec3 &b) {
        return a.x == b.x && a.y == b.y && a.z == b.z;
    };

    for (int i = 0; i < (int)circleIdx.size(); ++i) {
        if (used[i]) continue;
        int mi = circleIdx[i];
        int ci = muse[mi].cellIndex;
        Vec3 p = pos3d[ci];
        Vec3 n = frame3d[ci].n;

        Vec3 target = p;
        if (n.x != 0) {
            target.x = (p.x == dims.minX) ? dims.maxX : dims.minX;
        } else if (n.y != 0) {
            target.y = (p.y == dims.minY) ? dims.maxY : dims.minY;
        } else if (n.z != 0) {
            target.z = (p.z == dims.minZ) ? dims.maxZ : dims.minZ;
        } else {
            return false;
        }

        bool found = false;
        for (int j = i+1; j < (int)circleIdx.size(); ++j) {
            if (used[j]) continue;
            int mj = circleIdx[j];
            int cj = muse[mj].cellIndex;
            Vec3 pj = pos3d[cj];
            if (same3(target, pj)) {
                used[i] = used[j] = true;
                found = true;
                break;
            }
        }
        if (!found) return false;
    }

    // ---- Square constraints ----
    vector<int> squareCells;
    for (auto &m : muse) {
        if (m.type == 's' && m.cellIndex != -1)
            squareCells.push_back(m.cellIndex);
    }

    auto is_on_same_face = [&](int a, int b) -> bool {
        Vec3 pa = pos3d[a], pb = pos3d[b];
        if (pa.x == pb.x && (pa.x == dims.minX || pa.x == dims.maxX)) return true;
        if (pa.y == pb.y && (pa.y == dims.minY || pa.y == dims.maxY)) return true;
        if (pa.z == pb.z && (pa.z == dims.minZ || pa.z == dims.maxZ)) return true;
        return false;
    };

    for (int idx : squareCells) {
        Vec3 p = pos3d[idx];
        bool ok = false;
        for (int jdx : squareCells) {
            if (idx == jdx) continue;
            if (!is_on_same_face(idx, jdx)) continue;
            Vec3 q = pos3d[jdx];
            int dx = abs(p.x - q.x);
            int dy = abs(p.y - q.y);
            int dz = abs(p.z - q.z);
            if (dx + dy + dz == 1) {
                ok = true;
                break;
            }
        }
        if (!ok) return false;
    }

    return true;
}

// ----------------- Helpers for DFS -----------------

Vec3 step_from_frame(const Frame &f, Dir d) {
    switch(d) {
        case UP:    return {-f.v.x, -f.v.y, -f.v.z};
        case DOWN:  return { f.v.x,  f.v.y,  f.v.z};
        case LEFT:  return {-f.u.x, -f.u.y, -f.u.z};
        case RIGHT: return { f.u.x,  f.u.y,  f.u.z};
    }
    return {0,0,0};
}

// Map a normal vector to an index 0..5 (for counting)
int normal_id(const Vec3 &n) {
    if (n.x ==  1 && n.y == 0 && n.z == 0) return 0; // +X
    if (n.x == -1 && n.y == 0 && n.z == 0) return 1; // -X
    if (n.x == 0 && n.y ==  1 && n.z == 0) return 2; // +Y
    if (n.x == 0 && n.y == -1 && n.z == 0) return 3; // -Y
    if (n.x == 0 && n.y == 0 && n.z ==  1) return 4; // +Z
    if (n.x == 0 && n.y == 0 && n.z == -1) return 5; // -Z
    return -1; // shouldn't happen if frames are orthonormal
}

struct RectangularPrism {
    vector<Vec3> pos;
    vector<Frame> frame;
    vector<bool> assigned;
    unordered_map<Vec3, int, Vec3Hash> occ;
    array<int, 6> normalCount{};
    int assignedCount = 0;
    int maxFaceArea;
    int numCells;

    RectangularPrism(int n, int maxFace)
        : pos(n), frame(n), assigned(n, false), maxFaceArea(maxFace), numCells(n) {}

    void reset() {
        fill(assigned.begin(), assigned.end(), false);
        occ.clear();
        assignedCount = 0;
        normalCount.fill(0);
    }

    bool place_root(int idx, const Vec3 &p, const Frame &f) {
        reset();
        return place(idx, p, f);
    }

    bool place(int idx, const Vec3 &p, const Frame &f) {
        if (idx < 0 || idx >= numCells) return false;
        if (assigned[idx]) return false;
        if (occ.find(p) != occ.end()) return false;

        int nid = normal_id(f.n);
        if (nid < 0) return false;
        if (normalCount[nid] + 1 > maxFaceArea) return false;

        pos[idx] = p;
        frame[idx] = f;
        assigned[idx] = true;
        occ[p] = idx;
        assignedCount++;
        normalCount[nid]++;
        return true;
    }

    void unplace(int idx) {
        if (!assigned[idx]) return;
        int nid = normal_id(frame[idx].n);
        if (nid >= 0) normalCount[nid]--;
        occ.erase(pos[idx]);
        assigned[idx] = false;
        assignedCount--;
    }

    bool finalize(BoxDims &dims) {
        if (assignedCount != numCells) return false;
        return compute_box_dims_and_check(pos, numCells, dims);
    }
};

void print_partial_grid(const vector<Cell> &cells,
                        const vector<Frame> &frame,
                        const vector<bool> &assigned,
                        int R, int C)
{
    vector<string> g(R, string(C, '.'));

    for (int i = 0; i < (int)cells.size(); ++i) {
        int r = cells[i].r;
        int c = cells[i].c;

        if (!assigned[i]) {
            g[r][c] = '?';
            continue;
        }

        Vec3 n = frame[i].n;
        int f = 0;
        if      (n.x ==  1 && n.y == 0 && n.z == 0) f = 1;
        else if (n.x == -1 && n.y == 0 && n.z == 0) f = 2;
        else if (n.x == 0 && n.y ==  1 && n.z == 0) f = 3;
        else if (n.x == 0 && n.y == -1 && n.z == 0) f = 4;
        else if (n.x == 0 && n.y == 0 && n.z ==  1) f = 5;
        else if (n.x == 0 && n.y == 0 && n.z == -1) f = 6;
        else g[r][c] = 'X'; // shouldn't happen if frames are orthonormal

        if (f != 0)
            g[r][c] = char('0' + f);
    }

    cerr << "\n--- PARTIAL EMBEDDING (by normals) ---\n";
    for (int r = 0; r < R; ++r)
        cerr << g[r] << "\n";
    cerr << "--------------------------------------\n";
}

// ----------------- Core DFS solver (fast, for parallel search) -----------------

bool solve_one_fast(const vector<string> &grid,
                    const vector<Marker> &markers,
                    int solIndex,
                    mutex *ioMutex) {
    int R = (int)grid.size();
    int C = (int)grid[0].size();

    // Build cells and index map
    vector<Cell> cells;
    vector<vector<int>> idx(R, vector<int>(C, -1));
    for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
            if (grid[r][c] == '#') {
                int id = (int)cells.size();
                cells.push_back({r,c});
                idx[r][c] = id;
            }
        }
    }
    int numCells = (int)cells.size();
    if (numCells == 0) return false;
    int maxFaceArea = max_face_area_from_surface(numCells);

    // Build adjacency
    vector<vector<Edge>> adj(numCells);
    auto add_edge = [&](int r1, int c1, int r2, int c2, Dir d) {
        int a = idx[r1][c1];
        int b = idx[r2][c2];
        if (a == -1 || b == -1) return;
        adj[a].push_back({b, d});
    };
    for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
            if (idx[r][c] == -1) continue;
            if (r > 0)   add_edge(r,c,r-1,c,UP);
            if (r+1 < R) add_edge(r,c,r+1,c,DOWN);
            if (c > 0)   add_edge(r,c,r,c-1,LEFT);
            if (c+1 < C) add_edge(r,c,r,c+1,RIGHT);
        }
    }

    // Build reverse adjacency: for each cell b, revAdj[b] lists incoming neighbors nb
    vector<vector<Edge>> revAdj(numCells);
    for (int i = 0; i < numCells; ++i) {
        for (auto &e : adj[i]) {
            revAdj[e.to].push_back({i, e.dir}); // i -> e.to with dir e.dir
        }
    }

    // Embedding structures
    RectangularPrism prism(numCells, maxFaceArea);

    // Choose a root with the highest degree to reduce branching
    int root = 0;
    for (int i = 1; i < numCells; ++i) {
        if ((int)adj[i].size() > (int)adj[root].size()) root = i;
    }

    prism.place_root(root, {0,0,0}, { {1,0,0}, {0,1,0}, {0,0,1} });

    auto compute_candidate_from_neighbor =
        [&](int a, Dir d, int mode, Vec3 &candPos, Frame &candFrame) {
            const Frame &fa = prism.frame[a];
            const Vec3 &pa = prism.pos[a];
            Vec3 step = step_from_frame(fa, d);

            // The hinge axis is the edge shared by the two cells. For moves
            // along the grid's UP/DOWN, the edge runs along "u"; for LEFT/RIGHT
            // it runs along "v". We must rotate both the frame and the offset
            // vector around this axis to properly fold the neighbor out of the
            // base plane.
            Vec3 axis = (d == UP || d == DOWN) ? fa.u : fa.v;

            if (mode == 0) {
                // same face (no fold)
                candPos = pa + step;
                candFrame = fa;
            } else {
                Vec3 rotStep = (mode == 1) ? rotate90_plus(step, axis)
                                           : rotate90_minus(step, axis);
                candPos = pa + rotStep;

                if (mode == 1) {
                    candFrame.u = rotate90_plus(fa.u, axis);
                    candFrame.v = rotate90_plus(fa.v, axis);
                    candFrame.n = rotate90_plus(fa.n, axis);
                } else {
                    candFrame.u = rotate90_minus(fa.u, axis);
                    candFrame.v = rotate90_minus(fa.v, axis);
                    candFrame.n = rotate90_minus(fa.n, axis);
                }
            }
        };

    // Fast compatibility check using only assigned neighbors of b
    auto candidate_compatible =
        [&](int b, const Vec3 &candPos, const Frame &candFrame) -> bool {
        for (auto &e : revAdj[b]) {
            int nb = e.to;
            if (!prism.assigned[nb]) continue;

            Dir dnb = e.dir;

            bool okNeighbor = false;
            for (int mode = 0; mode < 3; ++mode) {
                Vec3 p2; Frame f2;
                compute_candidate_from_neighbor(nb, dnb, mode, p2, f2);
                if (p2 == candPos && f2 == candFrame) {
                    okNeighbor = true;
                    break;
                }
            }
            if (!okNeighbor) return false;
        }
        return true;
    };

    long long dfsIters = 0;
    const long long PRINT_EVERY = 100000000; // every million dfs calls

    std::function<bool()> dfs;
    dfs = [&]() -> bool {
        ++dfsIters;
        if (dfsIters % PRINT_EVERY == 0) {
            if (ioMutex) {
                lock_guard<mutex> lock(*ioMutex);
                cerr << "Sol #" << (solIndex + 1)
                     << " DFS iters (fast): " << dfsIters << "\n";

                print_partial_grid(cells, prism.frame, prism.assigned, R, C);
            } else {
                cerr << "Sol #" << (solIndex + 1)
                     << " DFS iters (fast): " << dfsIters << "\n";
                print_partial_grid(cells, prism.frame, prism.assigned, R, C);
            }
        }

        if (prism.assignedCount == numCells) {
            BoxDims dims;
            if (!prism.finalize(dims)) return false;
            if (!check_markers(markers, cells, prism.pos, prism.frame, dims))
                return false;
            return true; // success: box + markers satisfied
        }

        int b = -1, a = -1;
        Dir dirFromAtoB = UP;
        bool foundFrontier = false;

        for (int i = 0; i < numCells && !foundFrontier; ++i) {
            if (!prism.assigned[i]) continue;
            for (auto &e : adj[i]) {
                if (!prism.assigned[e.to]) {
                    a = i;
                    b = e.to;
                    dirFromAtoB = e.dir;
                    foundFrontier = true;
                    break;
                }
            }
        }
        if (!foundFrontier) {
            return false;
        }

        for (int mode = 0; mode < 3; ++mode) { // 0: same, 1: +90, 2: -90
            Vec3 candPos;
            Frame candFrame;
            compute_candidate_from_neighbor(a, dirFromAtoB, mode, candPos, candFrame);

            if (prism.occ.find(candPos) != prism.occ.end()) continue;
            if (!candidate_compatible(b, candPos, candFrame)) continue;

            if (prism.place(b, candPos, candFrame)) {
                if (dfs()) return true;
                prism.unplace(b);
            }
        }
        return false;
    };

    bool ok = dfs();
    if (ioMutex) {
        lock_guard<mutex> lock(*ioMutex);
        cerr << "Sol #" << (solIndex + 1)
             << " DFS iters (fast): " << dfsIters << " (done)\n";
    }
    return ok;
}

// ----------------- DFS solver with capture + progress (sequential, for final solution) -----------------

bool solve_one_capture(const vector<string> &grid,
                       const vector<Marker> &markers,
                       vector<Cell> &outCells,
                       vector<Vec3> &outPos,
                       BoxDims &outDims,
                       long long &outIters) {
    int R = (int)grid.size();
    int C = (int)grid[0].size();

    vector<Cell> cells;
    vector<vector<int>> idx(R, vector<int>(C, -1));
    for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
            if (grid[r][c] == '#') {
                int id = (int)cells.size();
                cells.push_back({r,c});
                idx[r][c] = id;
            }
        }
    }
    int numCells = (int)cells.size();
    if (numCells == 0) return false;
    int maxFaceArea = max_face_area_from_surface(numCells);

    vector<vector<Edge>> adj(numCells);
    auto add_edge = [&](int r1, int c1, int r2, int c2, Dir d) {
        int a = idx[r1][c1];
        int b = idx[r2][c2];
        if (a == -1 || b == -1) return;
        adj[a].push_back({b, d});
    };
    for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
            if (idx[r][c] == -1) continue;
            if (r > 0)   add_edge(r,c,r-1,c,UP);
            if (r+1 < R) add_edge(r,c,r+1,c,DOWN);
            if (c > 0)   add_edge(r,c,r,c-1,LEFT);
            if (c+1 < C) add_edge(r,c,r,c+1,RIGHT);
        }
    }

    // Reverse adjacency for capture solver as well
    vector<vector<Edge>> revAdj(numCells);
    for (int i = 0; i < numCells; ++i) {
        for (auto &e : adj[i]) {
            revAdj[e.to].push_back({i, e.dir});
        }
    }

    RectangularPrism prism(numCells, maxFaceArea);

    int root = 0;
    for (int i = 1; i < numCells; ++i) {
        if ((int)adj[i].size() > (int)adj[root].size()) root = i;
    }

    prism.place_root(root, {0,0,0}, { {1,0,0}, {0,1,0}, {0,0,1} });

    auto compute_candidate_from_neighbor =
        [&](int a, Dir d, int mode, Vec3 &candPos, Frame &candFrame) {
            const Frame &fa = prism.frame[a];
            const Vec3 &pa = prism.pos[a];
            Vec3 step = step_from_frame(fa, d);

            Vec3 axis = (d == UP || d == DOWN) ? fa.u : fa.v;

            if (mode == 0) {
                candPos = pa + step;
                candFrame = fa;
            } else {
                Vec3 rotStep = (mode == 1) ? rotate90_plus(step, axis)
                                           : rotate90_minus(step, axis);
                candPos = pa + rotStep;

                if (mode == 1) {
                    candFrame.u = rotate90_plus(fa.u, axis);
                    candFrame.v = rotate90_plus(fa.v, axis);
                    candFrame.n = rotate90_plus(fa.n, axis);
                } else {
                    candFrame.u = rotate90_minus(fa.u, axis);
                    candFrame.v = rotate90_minus(fa.v, axis);
                    candFrame.n = rotate90_minus(fa.n, axis);
                }
            }
        };

    auto candidate_compatible =
        [&](int b, const Vec3 &candPos, const Frame &candFrame) -> bool {
        for (auto &e : revAdj[b]) {
            int nb = e.to;
            if (!prism.assigned[nb]) continue;

            Dir dnb = e.dir;

            bool okNeighbor = false;
            for (int mode = 0; mode < 3; ++mode) {
                Vec3 p2; Frame f2;
                compute_candidate_from_neighbor(nb, dnb, mode, p2, f2);
                if (p2 == candPos && f2 == candFrame) {
                    okNeighbor = true;
                    break;
                }
            }
            if (!okNeighbor) return false;
        }
        return true;
    };

    long long dfsIters = 0;
    const long long PRINT_EVERY = 100000; // print every 100k DFS calls

    std::function<bool()> dfs;
    dfs = [&]() -> bool {
        ++dfsIters;
        if (dfsIters % PRINT_EVERY == 0) {
            cerr << "DFS iterations (capture): " << dfsIters << "\r";
            cerr.flush();
        }

        if (prism.assignedCount == numCells) {
            BoxDims dims;
            if (!prism.finalize(dims)) return false;
            if (!check_markers(markers, cells, prism.pos, prism.frame, dims))
                return false;

            outCells = cells;
            outPos = prism.pos;
            outDims = dims;
            outIters = dfsIters;
            return true;
        }

        int bestA = -1, bestB = -1;
        Dir bestDir = UP;
        int bestNeighborCount = -1;

        for (int i = 0; i < numCells; ++i) {
            if (!prism.assigned[i]) continue;
            for (auto &e : adj[i]) {
                int j = e.to;
                if (prism.assigned[j]) continue;

                int cnt = 0;
                for (auto &e2 : adj[j]) {
                    if (prism.assigned[e2.to]) ++cnt;
                }
                if (cnt > bestNeighborCount) {
                    bestNeighborCount = cnt;
                    bestA = i;
                    bestB = j;
                    bestDir = e.dir;
                }
            }
        }

        if (bestB == -1) {
            return false;
        }

        int a = bestA;
        int b = bestB;
        Dir dirFromAtoB = bestDir;

        for (int mode = 0; mode < 3; ++mode) {
            Vec3 candPos;
            Frame candFrame;
            compute_candidate_from_neighbor(a, dirFromAtoB, mode, candPos, candFrame);

            if (prism.occ.find(candPos) != prism.occ.end()) continue;
            if (!candidate_compatible(b, candPos, candFrame)) continue;

            if (prism.place(b, candPos, candFrame)) {
                if (dfs()) return true;
                prism.unplace(b);
            }
        }
        return false;
    };

    bool ok = dfs();
    cerr << "DFS iterations (capture): " << dfsIters << " (done)\n";
    return ok;
}

// ----------------- Read solutions.txt -----------------

vector<vector<string>> read_solutions(const string &filename) {
    ifstream in(filename);
    if (!in) {
        cerr << "Cannot open " << filename << "\n";
        exit(1);
    }
    vector<vector<string>> sols;
    string line;
    while (getline(in, line)) {
        if (line.rfind("Solution #", 0) == 0) {
            vector<string> grid;
            for (int i = 0; i < 20; ++i) {
                string g;
                if (!getline(in, g)) {
                    cerr << "Unexpected EOF while reading grid\n";
                    exit(1);
                }
                grid.push_back(g);
            }
            sols.push_back(grid);
        }
    }
    return sols;
}

// ----------------- Face labeling helper -----------------

int face_id(const Vec3 &p, const BoxDims &d) {
    if (p.z == d.minZ) return 5; // bottom
    if (p.z == d.maxZ) return 6; // top
    if (p.y == d.minY) return 3; // front
    if (p.y == d.maxY) return 4; // back
    if (p.x == d.minX) return 1; // left
    return 2;                    // right
}

// ----------------- Main (parallel search + sequential capture) -----------------

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Optional: run the 8x8 example provided in the prompt.
    if (argc > 1 && string(argv[1]) == "--example8") {
        vector<string> grid = {
            "........",
            "..##....",
            "######..",
            "..##.###",
            "####....",
            "...#....",
            "..#####.",
            "####.#.."
        };

        vector<Marker> markers = {
            {1,3,5,'n',-1}, {2,0,2,'n',-1}, {2,1,5,'s',-1}, {2,5,4,'n',-1},
            {3,7,2,'s',-1}, {3,0,2,'n',-1}, {5,3,5,'c',-1}, {6,3,6,'n',-1},
            {6,4,6,'c',-1}, {7,1,4,'c',-1}, {7,5,4,'c',-1}
        };

        vector<Cell> cells;
        vector<Vec3> pos;
        BoxDims dims;
        long long dfsIters = 0;

        cerr << "Solving 8x8 example net...\n";
        bool ok = solve_one_capture(grid, markers, cells, pos, dims, dfsIters);
        if (!ok) {
            cout << "No valid folding found for the 8x8 example.\n";
            return 0;
        }

        cout << "Found valid folding for the 8x8 example!\n\n";
        cout << "Face-labeled grid (1–6 for faces, . for empty):\n";
        vector<string> faceGrid(grid.size(), string(grid[0].size(), '.'));
        for (size_t i = 0; i < cells.size(); ++i) {
            int r = cells[i].r;
            int c = cells[i].c;
            int f = face_id(pos[i], dims);
            faceGrid[r][c] = char('0' + f);
        }
        for (auto &row : faceGrid) cout << row << "\n";
        cout << "\nTotal DFS calls: " << dfsIters << "\n";
        return 0;
    }

    // Markers (0-based indices: row, col, value, type)
    vector<Marker> markers = {
        {1,11,4,'n',-1}, {1,15,4,'c',-1}, {2,7,5,'s',-1},
        {3,12,7,'n',-1}, {3,14,5,'c',-1}, {4,10,4,'c',-1},
        {4,13,7,'n',-1}, {4,17,4,'c',-1},
        {5,6,4,'n',-1},  {5,8,7,'s',-1},  {6,7,9,'n',-1},
        {7,17,6,'s',-1}, {8,2,7,'c',-1},  {8,11,5,'s',-1},
        {9,14,5,'n',-1}, {10,5,4,'n',-1}, {10,7,7,'n',-1},
        {10,18,3,'s',-1}, {13,3,5,'s',-1}, {13,6,6,'n',-1},
        {13,9,2,'n',-1}, {15,6,5,'n',-1}, {17,12,5,'c',-1},
        {17,13,5,'n',-1}, {18,8,4,'s',-1}
    };

    auto solutions = read_solutions("solutions.txt");
    int nSol = (int)solutions.size();

    if (nSol == 0) {
        cout << "No grids in solutions.txt\n";
        return 0;
    }

    atomic<int> nextIndex{0};
    atomic<bool> found{false};
    atomic<int> bestSol{-1};
    atomic<int> solutionsChecked{0};
    mutex ioMutex;

    int numThreads = 8; // bump this if you want more parallelism
    cerr << "Starting parallel search on " << nSol
         << " candidate grids using " << numThreads << " threads.\n";

    vector<thread> threads;
    threads.reserve(numThreads);

    for (int t = 0; t < numThreads; ++t) {
        threads.emplace_back([&]() {
            while (!found.load(memory_order_relaxed)) {
                int i = nextIndex.fetch_add(1, memory_order_relaxed);
                if (i >= nSol) break;

                {
                    lock_guard<mutex> lock(ioMutex);
                    cerr << "Starting Solution #" << (i + 1) << "\n";
                }

                bool ok = solve_one_fast(solutions[i], markers, i, &ioMutex);

                int done = solutionsChecked.fetch_add(1) + 1;
                {
                    lock_guard<mutex> lock(ioMutex);
                    double frac = (double)done / (double)nSol;
                    int barWidth = 40;
                    int filled = (int)(frac * barWidth);
                    cerr << "[";
                    for (int j = 0; j < barWidth; ++j) {
                        cerr << (j < filled ? '#' : ' ');
                    }
                    cerr << "] " << done << "/" << nSol << "\r";
                    cerr.flush();
                }

                if (ok) {
                    bool expected = false;
                    if (found.compare_exchange_strong(expected, true)) {
                        bestSol.store(i, memory_order_relaxed);
                    }
                    break;
                }
            }
        });
    }

    for (auto &th : threads) th.join();

    {
        lock_guard<mutex> lock(ioMutex);
        cerr << "\nParallel search done.\n";
    }

    int idx = bestSol.load();
    if (idx == -1) {
        cout << "No candidate grid satisfied box + circle/square constraints.\n";
        return 0;
    }

    // Sequential capture solve with DFS progress
    vector<Cell> cells;
    vector<Vec3> pos;
    BoxDims dims;
    long long dfsIters = 0;

    cerr << "Re-solving Solution #" << (idx + 1)
         << " with capture + DFS progress...\n";

    bool ok = solve_one_capture(solutions[idx], markers, cells, pos, dims, dfsIters);
    if (!ok) {
        cout << "Internal error: capture solve failed for Solution #" << (idx + 1) << "\n";
        return 0;
    }

    cout << "Found valid box net with markers: Solution #" << (idx + 1) << "\n\n";

    cout << "Original grid:\n";
    for (const auto &row : solutions[idx]) {
        cout << row << "\n";
    }
    cout << "\n";

    const int R = 20, C = 20;
    vector<string> faceGrid(R, string(C, '.'));
    for (size_t i = 0; i < cells.size(); ++i) {
        int r = cells[i].r;
        int c = cells[i].c;
        int f = face_id(pos[i], dims);
        faceGrid[r][c] = char('0' + f);
    }

    // --- face area sanity check ---
    int faceCount[7] = {0}; // 1..6 used
    for (size_t i = 0; i < cells.size(); ++i) {
        int f = face_id(pos[i], dims);
        if (f >= 1 && f <= 6) {
            faceCount[f]++;
        }
    }

    int dx = dims.maxX - dims.minX + 1;
    int dy = dims.maxY - dims.minY + 1;
    int dz = dims.maxZ - dims.minZ + 1;

    int expected_x = dy * dz; // faces 1 & 2
    int expected_y = dx * dz; // faces 3 & 4
    int expected_z = dx * dy; // faces 5 & 6

    int totalCells = (int)cells.size();
    int totalFacesArea = 0;
    for (int f = 1; f <= 6; ++f) totalFacesArea += faceCount[f];

    cout << "Face areas (cell counts):\n";
    cout << "  Face 1 (x=minX): " << faceCount[1] << "\n";
    cout << "  Face 2 (x=maxX): " << faceCount[2] << "\n";
    cout << "  Face 3 (y=minY): " << faceCount[3] << "\n";
    cout << "  Face 4 (y=maxY): " << faceCount[4] << "\n";
    cout << "  Face 5 (z=minZ): " << faceCount[5] << "\n";
    cout << "  Face 6 (z=maxZ): " << faceCount[6] << "\n\n";

    cout << "Box dims from embedding: dx=" << dx
         << ", dy=" << dy << ", dz=" << dz << "\n";
    cout << "Expected areas from box dims:\n";
    cout << "  Faces 1 & 2 (dy*dz): " << expected_x << "\n";
    cout << "  Faces 3 & 4 (dx*dz): " << expected_y << "\n";
    cout << "  Faces 5 & 6 (dx*dy): " << expected_z << "\n\n";

    cout << "Check:\n";
    cout << "  face1 == face2? " << (faceCount[1] == faceCount[2] ? "OK" : "MISMATCH") << "\n";
    cout << "  face3 == face4? " << (faceCount[3] == faceCount[4] ? "OK" : "MISMATCH") << "\n";
    cout << "  face5 == face6? " << (faceCount[5] == faceCount[6] ? "OK" : "MISMATCH") << "\n";

    cout << "  face1/2 area match dy*dz? "
         << ((faceCount[1] == expected_x && faceCount[2] == expected_x) ? "OK" : "MISMATCH") << "\n";
    cout << "  face3/4 area match dx*dz? "
         << ((faceCount[3] == expected_y && faceCount[4] == expected_y) ? "OK" : "MISMATCH") << "\n";
    cout << "  face5/6 area match dx*dy? "
         << ((faceCount[5] == expected_z && faceCount[6] == expected_z) ? "OK" : "MISMATCH") << "\n";

    cout << "  total face area = " << totalFacesArea
         << ", total cells = " << totalCells
         << (totalFacesArea == totalCells ? " (OK)\n" : " (MISMATCH)\n");
    cout << "\n";

    cout << "Face-labeled grid (1–6 for faces, . for empty):\n";
    for (int r = 0; r < R; ++r) {
        cout << faceGrid[r] << "\n";
    }
    cout << "\n";

    cout << "Total DFS calls for this embedding: " << dfsIters << "\n";

    return 0;
}
