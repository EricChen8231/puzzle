#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include <array>
#include <cstdint>
#include <queue>

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

// ----------------- Helpers for DFS / debugging -----------------

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

void print_partial_grid(const vector<Cell> &cells,
                        const vector<Vec3> &pos,
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
        else g[r][c] = 'X'; // shouldn't happen

        if (f != 0)
            g[r][c] = char('0' + f);
    }

    cerr << "\n--- PARTIAL EMBEDDING (by normals) ---\n";
    for (int r = 0; r < R; ++r)
        cerr << g[r] << "\n";
    cerr << "--------------------------------------\n";
}

// ----------------- Box + marker checks -----------------

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
    if (markers.empty()) return true;

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

// ----------------- Face labeling helper -----------------

int face_id(const Vec3 &p, const BoxDims &d) {
    if (p.z == d.minZ) return 5; // bottom
    if (p.z == d.maxZ) return 6; // top
    if (p.y == d.minY) return 3; // front
    if (p.y == d.maxY) return 4; // back
    if (p.x == d.minX) return 1; // left
    return 2;                    // right
}

// ----------------- New edge-tree BFS/DFS solver -----------------

bool solve_one_edge_tree(const vector<string> &grid,
                         const vector<Marker> &markers,
                         vector<Cell> &outCells,
                         vector<Vec3> &outPos,
                         BoxDims &outDims,
                         long long &outIters)
{
    int R = (int)grid.size();
    int C = (int)grid[0].size();

    // ----- 1. Build cells and index map -----
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
    if (numCells == 0) {
        outIters = 0;
        return false;
    }

    // ----- 2. Build adjacency lists (directed) -----
    struct EdgeLocal {
        int to;
        Dir dir;
    };

    vector<vector<EdgeLocal>> adj(numCells);
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

    // Reverse adjacency (for consistency checks)
    vector<vector<EdgeLocal>> revAdj(numCells);
    for (int i = 0; i < numCells; ++i) {
        for (auto &e : adj[i]) {
            revAdj[e.to].push_back({i, e.dir});
        }
    }

    // ----- 3. Build a BFS spanning tree -----
    vector<int> parent(numCells, -1);
    vector<Dir> parentDir(numCells, UP);  // direction from parent -> node
    vector<int> bfsOrder;
    bfsOrder.reserve(numCells);

    {
        vector<bool> vis(numCells, false);
        queue<int> q;
        int root = 0;
        vis[root] = true;
        q.push(root);
        bfsOrder.push_back(root);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (auto &e : adj[u]) {
                int v = e.to;
                if (!vis[v]) {
                    vis[v] = true;
                    parent[v] = u;
                    parentDir[v] = e.dir;
                    q.push(v);
                    bfsOrder.push_back(v);
                }
            }
        }
        if ((int)bfsOrder.size() != numCells) {
            cerr << "Warning: net graph is disconnected; cannot embed.\n";
            outIters = 0;
            return false;
        }
    }

    // ----- 4. Embedding state -----
    vector<bool> assigned(numCells, false);
    vector<Vec3> pos3d(numCells);
    vector<Frame> frame3d(numCells);
    unordered_map<Vec3,int,Vec3Hash> occ;
    occ.reserve(numCells * 4);
    vector<int> normalCount(6, 0);

    // root placement
    int root = bfsOrder[0];
    assigned[root] = true;
    pos3d[root] = {0,0,0};
    frame3d[root] = { {1,0,0}, {0,1,0}, {0,0,1} };
    occ[pos3d[root]] = root;
    int assignedCount = 1;
    {
        int nid = normal_id(frame3d[root].n);
        if (nid >= 0) normalCount[nid]++;
    }

    // helper: folding transition from parent a -> child b along dir d and mode
    auto compute_candidate_from_parent =
        [&](int a, int b, Dir d, int mode, Vec3 &candPos, Frame &candFrame) {
            const Frame &fa = frame3d[a];
            const Vec3 &pa = pos3d[a];

            Vec3 u = fa.u;
            Vec3 v = fa.v;
            Vec3 n = fa.n;

            Vec3 off = step_from_frame(fa, d);

            Vec3 u2 = u, v2 = v, n2 = n, off2 = off;

            if (mode != 0) {
                if (d == UP || d == DOWN) {
                    Vec3 axis = u;
                    if (mode == 1) {
                        v2   = rotate90_plus(v2, axis);
                        n2   = rotate90_plus(n2, axis);
                        off2 = rotate90_plus(off2, axis);
                    } else {
                        v2   = rotate90_minus(v2, axis);
                        n2   = rotate90_minus(n2, axis);
                        off2 = rotate90_minus(off2, axis);
                    }
                } else { // LEFT or RIGHT: hinge axis along v
                    Vec3 axis = v;
                    if (mode == 1) {
                        u2   = rotate90_plus(u2, axis);
                        n2   = rotate90_plus(n2, axis);
                        off2 = rotate90_plus(off2, axis);
                    } else {
                        u2   = rotate90_minus(u2, axis);
                        n2   = rotate90_minus(n2, axis);
                        off2 = rotate90_minus(off2, axis);
                    }
                }
            }

            candPos     = pa + off2;
            candFrame.u = u2;
            candFrame.v = v2;
            candFrame.n = n2;
        };

    // candidate must be consistent with *all* already-assigned neighbors via their edges
    auto candidate_compatible =
        [&](int b, const Vec3 &candPos, const Frame &candFrame) -> bool {
        for (auto &e : revAdj[b]) {
            int nb = e.to;
            if (!assigned[nb]) continue;

            Dir dnb = e.dir; // nb -> b

            bool okNeighbor = false;
            for (int mode = 0; mode < 3; ++mode) {
                Vec3 p2; Frame f2;
                compute_candidate_from_parent(nb, b, dnb, mode, p2, f2);
                if (p2.x == candPos.x && p2.y == candPos.y && p2.z == candPos.z &&
                    f2.u == candFrame.u && f2.v == candFrame.v && f2.n == candFrame.n) {
                    okNeighbor = true;
                    break;
                }
            }
            if (!okNeighbor) return false;
        }
        return true;
    };

    long long dfsIters = 0;
    const long long PRINT_EVERY = 1000000;
    static int rejectPrints = 0;
    const int MAX_REJECT_PRINTS = 10;
    int halfTotal = numCells / 2;

    // ----- 5. DFS over BFS order of cells -----
    function<bool(int)> dfs = [&](int idxInOrder) -> bool {
        ++dfsIters;
        if (dfsIters % PRINT_EVERY == 0) {
            cerr << "DFS iterations (edge-tree): " << dfsIters << "\r";
            cerr.flush();
        }

        if (idxInOrder == numCells) {
            // all cells placed, now enforce actual box + markers

            BoxDims dims;
            bool boxOK = compute_box_dims_and_check(pos3d, numCells, dims);
            bool markersOK = check_markers(markers, cells, pos3d, frame3d, dims);

            if (!boxOK || !markersOK) {
                if (rejectPrints < MAX_REJECT_PRINTS) {
                    ++rejectPrints;
                    cerr << "\n=== Full embedding reached but rejected ===\n";
                    cerr << "boxOK = " << boxOK << ", markersOK = " << markersOK << "\n";
                    cerr << "Dims: "
                         << "x:[" << dims.minX << "," << dims.maxX << "], "
                         << "y:[" << dims.minY << "," << dims.maxY << "], "
                         << "z:[" << dims.minZ << "," << dims.maxZ << "]\n";
                    print_partial_grid(cells, pos3d, frame3d, assigned, R, C);
                    cerr << "==========================================\n";
                }
                return false;
            }

            outCells = cells;
            outPos   = pos3d;
            outDims  = dims;
            outIters = dfsIters;
            return true;
        }

        int b = bfsOrder[idxInOrder];
        if (b == root) {
            // root already placed, skip to next
            return dfs(idxInOrder + 1);
        }

        int a = parent[b];
        Dir dirAB = parentDir[b];

        if (a < 0 || !assigned[a]) {
            // shouldn't happen
            return false;
        }

        // Try modes for this parent edge
        for (int mode = 0; mode < 3; ++mode) {
            Vec3 candPos;
            Frame candFrame;
            compute_candidate_from_parent(a, b, dirAB, mode, candPos, candFrame);

            if (occ.find(candPos) != occ.end()) continue;
            if (!candidate_compatible(b, candPos, candFrame)) continue;

            assigned[b] = true;
            pos3d[b] = candPos;
            frame3d[b] = candFrame;
            occ[candPos] = b;
            ++assignedCount;

            int nid = normal_id(candFrame.n);
            if (nid >= 0) normalCount[nid]++;

            bool tooBigFace = false;
            for (int k = 0; k < 6; ++k) {
                if (normalCount[k] > halfTotal) {
                    tooBigFace = true;
                    break;
                }
            }

            if (tooBigFace) {
                if (nid >= 0) normalCount[nid]--;
                --assignedCount;
                occ.erase(candPos);
                assigned[b] = false;
                continue;
            }

            if (dfs(idxInOrder + 1)) return true;

            if (nid >= 0) normalCount[nid]--;
            --assignedCount;
            occ.erase(candPos);
            assigned[b] = false;
        }

        return false;
    };

    bool ok = dfs(0);
    cerr << "\nDFS iterations (edge-tree): " << dfsIters << " (done)\n";
    outIters = dfsIters;
    return ok;
}

// ----------------- Main: 8x8 test grid + markers -----------------

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // 8x8 test grid
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

    // Markers: row, col, value, type, cellIndex=-1 (0-based)
    vector<Marker> markers = {
        {1,3,5,'n',-1},
        {2,0,2,'n',-1},
        {2,1,5,'s',-1},
        {2,5,4,'n',-1},
        {3,7,2,'s',-1},
        {3,0,2,'n',-1},
        {5,3,5,'c',-1},
        {6,3,6,'n',-1},
        {6,4,6,'c',-1},
        {7,1,4,'c',-1},
        {7,5,4,'c',-1}
    };

    vector<Cell> cells;
    vector<Vec3> pos;
    BoxDims dims;
    long long dfsIters = 0;

    cerr << "Solving 8x8 test grid WITH markers (edge-tree DFS)...\n";
    bool ok = solve_one_edge_tree(grid, markers, cells, pos, dims, dfsIters);
    if (!ok) {
        cout << "No valid box embedding for this test grid that satisfies box + marker constraints.\n";
        cout << "Total DFS calls: " << dfsIters << "\n";
        return 0;
    }

    cout << "Found valid box net WITH markers.\n";
    cout << "Total DFS calls: " << dfsIters << "\n\n";

    int R = (int)grid.size();
    int C = (int)grid[0].size();
    vector<string> faceGrid(R, string(C, '.'));
    for (size_t i = 0; i < cells.size(); ++i) {
        int r = cells[i].r;
        int c = cells[i].c;
        int f = face_id(pos[i], dims);
        faceGrid[r][c] = char('0' + f);
    }

    int faceCount[7] = {0};
    for (size_t i = 0; i < cells.size(); ++i) {
        int f = face_id(pos[i], dims);
        if (f >= 1 && f <= 6) faceCount[f]++;
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

    cout << "Face-labeled grid WITH markers (1–6 for faces, . for empty):\n";
    for (int r = 0; r < R; ++r) {
        cout << faceGrid[r] << "\n";
    }
    cout << "\n";

    return 0;
}
