#include <algorithm>
#include <array>
#include <atomic>
#include <climits>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

namespace {

constexpr bool ENFORCE_MARKER_RULES = true;
const long long REPORT_EVERY = 1000000; // DFS progress print interval

// Global stop flag and logging mutex
atomic<bool> g_stop{false};
mutex g_cerrMutex;
atomic<int> g_winnerIndex{-1};

struct Marker {
    int row = 0;
    int col = 0;
    char type = 's'; // 's' square, 'c' circle
};

struct Net {
    string header;
    vector<string> rows;
    int index = -1;
    int totalTiles = 0;
    vector<Marker> markers;
};

struct ProcessResult {
    bool valid = false;
    string header;
    int index = -1;
    int totalTiles = 0;
    int Lx = 0, Ly = 0, Lz = 0;
    int rootPatch = -1;
    vector<string> labeledGrid;
    vector<string> originalGrid;
    vector<string> mappingLines;
    string failureReason;
    long long dfsIterations = 0;
};

struct Patch {
    int axis;
    int sign;
    int i;
    int j;
    int nx;
    int ny;
    int nz;
};

struct PrismData {
    int Lx = 0;
    int Ly = 0;
    int Lz = 0;
    vector<Patch> patches;
    vector<vector<int>> adj;
    vector<vector<char>> adjMatrix;
    vector<int> oppPatch;  // index of opposite patch for each patch
};

int faceLabel(const Patch &p) {
    if (p.axis == 0) return (p.sign == -1 ? 1 : 2);
    if (p.axis == 1) return (p.sign == -1 ? 3 : 4);
    return (p.sign == -1 ? 5 : 6);
}

class FoldingSolver {
public:
    explicit FoldingSolver(vector<string> gridInput, int netIndex)
        : netIndex(netIndex) {
        if (gridInput.empty()) {
            H = 0;
            W = 0;
            return;
        }
        H = static_cast<int>(gridInput.size());
        W = 0;
        for (const auto &row : gridInput) {
            W = max(W, static_cast<int>(row.size()));
        }
        grid = std::move(gridInput);
        for (auto &row : grid) {
            if (static_cast<int>(row.size()) < W) {
                row += string(W - row.size(), '.');
            }
        }
        buildTilesAndAdj();

        tileMarker.assign(T, 'n');
        circleTiles.clear();
        squareTiles.clear();
    }

    // markers: vector of (row, col, type), type in {'c','s'}
    void setMarkers(const vector<tuple<int,int,char>> &marks) {
        tileMarker.assign(T, 'n');
        circleTiles.clear();
        squareTiles.clear();

        for (auto &mk : marks) {
            int r, c;
            char tp;
            tie(r, c, tp) = mk;

            int tid = -1;
            for (int t = 0; t < T; ++t) {
                if (tiles[t].first == r && tiles[t].second == c) {
                    tid = t;
                    break;
                }
            }
            if (tid == -1) continue; // marker on '.'

            if (tp == 'c') {
                tileMarker[tid] = 'c';
                circleTiles.push_back(tid);
            } else if (tp == 's') {
                tileMarker[tid] = 's';
                squareTiles.push_back(tid);
            }
        }
    }

    bool solveForDims(int lx, int ly, int lz,
                      vector<int> &assignmentOut,
                      int &rootPatchUsed) {
        dfsIterations = 0;
        if (T == 0) return false;
        PrismData prism = buildPrism(lx, ly, lz);
        if (static_cast<int>(prism.patches.size()) != T) {
            return false;
        }

        tileToPatch.assign(T, -1);
        patchUsed.assign(prism.patches.size(), 0);
        solutionFound = false;
        finalAssignment.clear();
        allSolutions.clear();
        currentPrism = &prism;
        dimLx = lx;
        dimLy = ly;
        dimLz = lz;

        const int P = static_cast<int>(prism.patches.size());
        canUse.assign(T, vector<char>(P, 1));
        candCount.assign(T, 0);

        for (int t = 0; t < T; ++t) {
            int tileDeg = static_cast<int>(tileAdj[t].size());
            int count = 0;
            for (int p = 0; p < P; ++p) {
                int patchDeg = static_cast<int>(prism.adj[p].size());
                if (patchDeg < tileDeg) {
                    canUse[t][p] = 0;
                } else {
                    canUse[t][p] = 1;
                    ++count;
                }
            }
            candCount[t] = count;
        }

        const int rootTile = 0;
        for (int rootPatch = 0; rootPatch < P; ++rootPatch) {
            if (g_stop.load(memory_order_relaxed)) return false;

            fill(tileToPatch.begin(), tileToPatch.end(), -1);
            fill(patchUsed.begin(), patchUsed.end(), 0);
            solutionFound = false;

            for (int t = 0; t < T; ++t) {
                int tileDeg = static_cast<int>(tileAdj[t].size());
                int count = 0;
                for (int p = 0; p < P; ++p) {
                    int patchDeg = static_cast<int>(prism.adj[p].size());
                    if (patchDeg < tileDeg) {
                        canUse[t][p] = 0;
                    } else {
                        canUse[t][p] = 1;
                        ++count;
                    }
                }
                candCount[t] = count;
            }

            tileToPatch[rootTile] = rootPatch;
            patchUsed[rootPatch] = 1;

            if (!checkCirclePartial(rootTile, rootPatch) ||
                !checkSquarePartial(rootTile, rootPatch)) {
                tileToPatch[rootTile] = -1;
                patchUsed[rootPatch] = 0;
                continue;
            }

            vector<pair<int,int>> changes;
            if (!propagateConstraints(rootTile, rootPatch, changes)) {
                tileToPatch[rootTile] = -1;
                patchUsed[rootPatch] = 0;
                continue;
            }

            if (dfsAssign(1, rootPatch)) {
                if (ENFORCE_MARKER_RULES) {
                    assignmentOut = finalAssignment;
                    rootPatchUsed = rootPatch;
                    return true;
                }
            }
        }
        return false;
    }

    vector<string> buildLabeledGrid(const vector<int> &assignment,
                                    int lx, int ly, int lz) const {
        PrismData prism = buildPrism(lx, ly, lz);
        vector<string> labeled(H, string(W, '.'));
        for (int t = 0; t < T; ++t) {
            auto [r, c] = tiles[t];
            const Patch &p = prism.patches[assignment[t]];
            labeled[r][c] = static_cast<char>('0' + faceLabel(p));
        }
        return labeled;
    }

    vector<string> describeAssignment(const vector<int> &assignment,
                                      int lx, int ly, int lz) const {
        PrismData prism = buildPrism(lx, ly, lz);
        vector<string> lines;
        lines.reserve(T);
        for (int t = 0; t < T; ++t) {
            auto [r, c] = tiles[t];
            const Patch &p = prism.patches[assignment[t]];
            ostringstream oss;
            oss << "Tile (" << r << "," << c << ") -> axis=" << p.axis
                << ", sign=" << p.sign
                << ", i=" << p.i
                << ", j=" << p.j
                << ", n=(" << p.nx << "," << p.ny << "," << p.nz << ")";
            if (tileMarker[t] == 'c') oss << " [circle]";
            if (tileMarker[t] == 's') oss << " [square]";
            lines.push_back(oss.str());
        }
        return lines;
    }

    int tileCount() const { return T; }
    int getSolutionCount() const { return static_cast<int>(allSolutions.size()); }
    long long getIterationCount() const { return dfsIterations; }

private:
    int H = 0;
    int W = 0;
    int T = 0;
    vector<string> grid;
    vector<pair<int,int>> tiles;
    vector<vector<int>> tileAdj;

    mutable vector<int> tileToPatch;
    mutable vector<char> patchUsed;
    mutable bool solutionFound = false;
    mutable vector<int> finalAssignment;
    mutable vector<pair<vector<int>, int>> allSolutions;
    mutable const PrismData *currentPrism = nullptr;

    mutable vector<vector<char>> canUse;
    mutable vector<int> candCount;

    mutable vector<char> tileMarker;     // 'n', 'c', 's'
    mutable vector<int> circleTiles;
    mutable vector<int> squareTiles;

    mutable long long dfsIterations = 0;
    int netIndex = -1;
    int dimLx = 0, dimLy = 0, dimLz = 0;

    void buildTilesAndAdj() {
        tiles.clear();
        vector<vector<int>> coord(H, vector<int>(W, -1));
        for (int r = 0; r < H; ++r) {
            for (int c = 0; c < W; ++c) {
                if (grid[r][c] == '#') {
                    coord[r][c] = static_cast<int>(tiles.size());
                    tiles.push_back({r, c});
                }
            }
        }
        T = static_cast<int>(tiles.size());
        tileAdj.assign(T, {});
        const int dr[4] = {-1, 1, 0, 0};
        const int dc[4] = {0, 0, -1, 1};
        for (int t = 0; t < T; ++t) {
            auto [r, c] = tiles[t];
            for (int k = 0; k < 4; ++k) {
                int nr = r + dr[k];
                int nc = c + dc[k];
                if (nr < 0 || nr >= H || nc < 0 || nc >= W) continue;
                int nt = coord[nr][nc];
                if (nt != -1) tileAdj[t].push_back(nt);
            }
        }
    }

    PrismData buildPrism(int lx, int ly, int lz) const {
        PrismData data;
        data.Lx = lx;
        data.Ly = ly;
        data.Lz = lz;

        auto addPatches = [&](int axis, int sign, int d1, int d2) {
            for (int i = 0; i < d1; ++i) {
                for (int j = 0; j < d2; ++j) {
                    Patch p;
                    p.axis = axis;
                    p.sign = sign;
                    p.i = i;
                    p.j = j;
                    p.nx = (axis == 0 ? (sign == -1 ? -1 : 1) : 0);
                    p.ny = (axis == 1 ? (sign == -1 ? -1 : 1) : 0);
                    p.nz = (axis == 2 ? (sign == -1 ? -1 : 1) : 0);
                    data.patches.push_back(p);
                }
            }
        };

        addPatches(0, -1, ly, lz);
        addPatches(0, +1, ly, lz);
        addPatches(1, -1, lx, lz);
        addPatches(1, +1, lx, lz);
        addPatches(2, -1, lx, ly);
        addPatches(2, +1, lx, ly);

        const int P = static_cast<int>(data.patches.size());
        data.adj.assign(P, {});

        const int vx = lx + 1;
        const int vy = ly + 1;
        const int vz = lz + 1;
        vector<int> vertexId(vx * vy * vz, -1);
        int vCount = 0;
        auto vertexIndex = [&](int x, int y, int z) -> int {
            return ((x * vy) + y) * vz + z;
        };
        auto getVertexId = [&](int x, int y, int z) -> int {
            int idx = vertexIndex(x, y, z);
            int &ref = vertexId[idx];
            if (ref == -1) ref = vCount++;
            return ref;
        };

        unordered_map<long long, vector<int>> edgeOwners;
        edgeOwners.reserve(P * 4);

        auto addEdge = [&](int patchIdx, int a, int b) {
            if (a > b) swap(a, b);
            long long key =
                (static_cast<long long>(a) << 32) |
                static_cast<unsigned int>(b);
            edgeOwners[key].push_back(patchIdx);
        };

        for (int idx = 0; idx < P; ++idx) {
            const Patch &p = data.patches[idx];
            vector<array<int,3>> corners(4);
            if (p.axis == 0) {
                int x = (p.sign == -1 ? 0 : lx);
                int y = p.i;
                int z = p.j;
                corners[0] = {x, y, z};
                corners[1] = {x, y+1, z};
                corners[2] = {x, y+1, z+1};
                corners[3] = {x, y,   z+1};
            } else if (p.axis == 1) {
                int y = (p.sign == -1 ? 0 : ly);
                int x = p.i;
                int z = p.j;
                corners[0] = {x,   y, z};
                corners[1] = {x+1, y, z};
                corners[2] = {x+1, y, z+1};
                corners[3] = {x,   y, z+1};
            } else {
                int z = (p.sign == -1 ? 0 : lz);
                int x = p.i;
                int y = p.j;
                corners[0] = {x,   y,   z};
                corners[1] = {x+1, y,   z};
                corners[2] = {x+1, y+1, z};
                corners[3] = {x,   y+1, z};
            }

            vector<int> vids(4);
            for (int k = 0; k < 4; ++k) {
                vids[k] = getVertexId(corners[k][0],
                                      corners[k][1],
                                      corners[k][2]);
            }
            for (int k = 0; k < 4; ++k) {
                addEdge(idx, vids[k], vids[(k+1)%4]);
            }
        }

        for (auto &kv : edgeOwners) {
            const vector<int> &owners = kv.second;
            for (size_t i = 0; i < owners.size(); ++i) {
                for (size_t j = i+1; j < owners.size(); ++j) {
                    int a = owners[i];
                    int b = owners[j];
                    data.adj[a].push_back(b);
                    data.adj[b].push_back(a);
                }
            }
        }

        for (auto &lst : data.adj) {
            sort(lst.begin(), lst.end());
            lst.erase(unique(lst.begin(), lst.end()), lst.end());
        }

        data.adjMatrix.assign(P, vector<char>(P, 0));
        for (int i = 0; i < P; ++i) {
            for (int j : data.adj[i]) {
                data.adjMatrix[i][j] = 1;
            }
        }

        data.oppPatch.assign(P, -1);
        for (int i = 0; i < P; ++i) {
            const Patch &A = data.patches[i];
            for (int j = 0; j < P; ++j) {
                const Patch &B = data.patches[j];
                if (A.axis == B.axis &&
                    A.i    == B.i    &&
                    A.j    == B.j    &&
                    A.sign == -B.sign) {
                    data.oppPatch[i] = j;
                    break;
                }
            }
        }

        return data;
    }

    int chooseNextTile() const {
        int best = -1;
        int bestAvail = INT_MAX;
        int bestAssignedNbrs = -1;
        int bestDeg = -1;

        const int P = static_cast<int>(currentPrism->patches.size());

        for (int t = 0; t < T; ++t) {
            if (tileToPatch[t] != -1) continue;

            int assignedNbrs = 0;
            for (int nb : tileAdj[t]) {
                if (tileToPatch[nb] != -1) ++assignedNbrs;
            }
            int deg = static_cast<int>(tileAdj[t].size());

            int avail = 0;
            for (int p = 0; p < P; ++p) {
                if (canUse[t][p] && !patchUsed[p]) ++avail;
            }

            if (avail < bestAvail ||
                (avail == bestAvail && assignedNbrs > bestAssignedNbrs) ||
                (avail == bestAvail && assignedNbrs == bestAssignedNbrs && deg > bestDeg)) {
                bestAvail = avail;
                bestAssignedNbrs = assignedNbrs;
                bestDeg = deg;
                best = t;
            }
        }
        return best;
    }

    bool checkCirclePartial(int tileIdx, int patchIdx) const {
        if (!ENFORCE_MARKER_RULES) return true;
        if (tileMarker.empty() || tileMarker[tileIdx] != 'c') return true;

        int op = currentPrism->oppPatch[patchIdx];
        if (op < 0) return false;

        for (int t = 0; t < T; ++t) {
            if (tileToPatch[t] == op) {
                if (tileMarker[t] != 'c') return false;
                return true;
            }
        }

        for (int t : circleTiles) {
            if (t == tileIdx) continue;
            if (tileToPatch[t] != -1) continue;
            if (canUse[t][op] && !patchUsed[op]) {
                return true;
            }
        }
        return false;
    }

    bool checkSquarePartial(int tileIdx, int patchIdx) const {
        if (!ENFORCE_MARKER_RULES) return true;
        if (tileMarker.empty() || tileMarker[tileIdx] != 's') return true;
        if (squareTiles.size() < 2) return false;

        const Patch &P = currentPrism->patches[patchIdx];

        for (int t : squareTiles) {
            if (t == tileIdx) continue;
            int p2 = tileToPatch[t];
            if (p2 == -1) continue;
            const Patch &Q = currentPrism->patches[p2];
            if (Q.axis == P.axis && Q.sign == P.sign &&
                currentPrism->adjMatrix[patchIdx][p2]) {
                return true;
            }
        }

        const int Pcount = static_cast<int>(currentPrism->patches.size());
        for (int t : squareTiles) {
            if (t == tileIdx) continue;
            if (tileToPatch[t] != -1) continue;

            for (int q = 0; q < Pcount; ++q) {
                if (!canUse[t][q]) continue;
                if (patchUsed[q]) continue;
                const Patch &Q = currentPrism->patches[q];
                if (Q.axis != P.axis || Q.sign != P.sign) continue;
                if (!currentPrism->adjMatrix[patchIdx][q]) continue;
                return true;
            }
        }
        return false;
    }

    bool finalCircleCheck() const {
        if (!ENFORCE_MARKER_RULES) return true;
        for (int t : circleTiles) {
            int p = tileToPatch[t];
            if (p < 0) return false;
            int op = currentPrism->oppPatch[p];
            if (op < 0) return false;
            bool ok = false;
            for (int t2 : circleTiles) {
                if (t2 == t) continue;
                if (tileToPatch[t2] == op) {
                    ok = true;
                    break;
                }
            }
            if (!ok) return false;
        }
        return true;
    }

    bool finalSquareCheck() const {
        if (!ENFORCE_MARKER_RULES) return true;
        for (int t : squareTiles) {
            int p = tileToPatch[t];
            if (p < 0) return false;
            const Patch &P = currentPrism->patches[p];
            bool ok = false;
            for (int t2 : squareTiles) {
                if (t2 == t) continue;
                int p2 = tileToPatch[t2];
                if (p2 < 0) return false;
                const Patch &Q = currentPrism->patches[p2];
                if (Q.axis == P.axis && Q.sign == P.sign &&
                    currentPrism->adjMatrix[p][p2]) {
                    ok = true;
                    break;
                }
            }
            if (!ok) return false;
        }
        return true;
    }

    bool propagateConstraints(int tileIdx, int patchIdx,
                              vector<pair<int,int>> &changes) const {
        const int P = static_cast<int>(currentPrism->patches.size());

        for (int nb : tileAdj[tileIdx]) {
            if (tileToPatch[nb] != -1) continue;

            for (int q = 0; q < P; ++q) {
                if (!canUse[nb][q]) continue;
                if (!currentPrism->adjMatrix[q][patchIdx]) {
                    canUse[nb][q] = 0;
                    --candCount[nb];
                    changes.emplace_back(nb, q);
                }
            }
            if (candCount[nb] == 0) return false;

            bool hasFree = false;
            for (int q = 0; q < P; ++q) {
                if (canUse[nb][q] && !patchUsed[q]) {
                    hasFree = true;
                    break;
                }
            }
            if (!hasFree) return false;
        }
        return true;
    }

    bool dfsAssign(int assignedCount, int rootPatch) const {
        if (g_stop.load(memory_order_relaxed)) return false;

        ++dfsIterations;
        if (dfsIterations % REPORT_EVERY == 0) {
            lock_guard<mutex> lock(g_cerrMutex);
            cerr << "\r[Net " << (netIndex + 1)
                 << " " << dimLx << "x" << dimLy << "x" << dimLz
                 << "] DFS calls: " << dfsIterations << flush;
        }

        if (assignedCount == T) {
            if (!finalCircleCheck()) return false;
            if (!finalSquareCheck()) return false;

            allSolutions.emplace_back(tileToPatch, rootPatch);

            if (ENFORCE_MARKER_RULES) {
                solutionFound = true;
                finalAssignment = tileToPatch;
                return true;
            }
            return false;
        }

        int tileIdx = chooseNextTile();
        if (tileIdx == -1) return false;

        const int patchCount = static_cast<int>(currentPrism->patches.size());

        for (int p = 0; p < patchCount; ++p) {
            if (!canUse[tileIdx][p]) continue;
            if (patchUsed[p]) continue;

            bool ok = true;
            for (int nb : tileAdj[tileIdx]) {
                int assignedPatch = tileToPatch[nb];
                if (assignedPatch == -1) continue;
                if (!currentPrism->adjMatrix[p][assignedPatch]) {
                    ok = false;
                    break;
                }
            }
            if (!ok) continue;

            if (!checkCirclePartial(tileIdx, p)) continue;
            if (!checkSquarePartial(tileIdx, p)) continue;

            tileToPatch[tileIdx] = p;
            patchUsed[p] = 1;

            vector<pair<int,int>> changes;
            if (propagateConstraints(tileIdx, p, changes)) {
                if (dfsAssign(assignedCount + 1, -1)) return true;
            }

            for (auto &chg : changes) {
                int t = chg.first;
                int q = chg.second;
                if (!canUse[t][q]) {
                    canUse[t][q] = 1;
                    ++candCount[t];
                }
            }

            tileToPatch[tileIdx] = -1;
            patchUsed[p] = 0;

            if (g_stop.load(memory_order_relaxed)) return false;
        }
        return false;
    }
};

vector<tuple<int,int,int>> enumerateTriples(int halfArea) {
    vector<tuple<int,int,int>> triples;
    for (int a = 1; a <= halfArea; ++a) {
        for (int b = a; b <= halfArea; ++b) {
            int denom = a + b;
            int numer = halfArea - a * b;
            if (denom <= 0 || numer <= 0) continue;
            if (numer % denom != 0) continue;
            int c = numer / denom;
            if (c < b) continue;
            triples.emplace_back(a, b, c);
        }
    }
    return triples;
}

vector<Net> loadNets(const string &path) {
    ifstream in(path);
    vector<Net> nets;
    if (!in) {
        lock_guard<mutex> lock(g_cerrMutex);
        cerr << "Failed to open " << path << " for reading.\n";
        return nets;
    }

    Net current;
    string line;
    while (getline(in, line)) {
        if (line.rfind("Solution #", 0) == 0) {
            if (!current.rows.empty()) {
                nets.push_back(current);
                current = Net();
            }
            current.header = line;
        } else if (line.empty()) {
            if (!current.rows.empty()) {
                nets.push_back(current);
                current = Net();
            }
        } else {
            current.rows.push_back(line);
        }
    }
    if (!current.rows.empty()) {
        nets.push_back(current);
    }

    for (size_t idx = 0; idx < nets.size(); ++idx) {
        Net &net = nets[idx];
        net.index = static_cast<int>(idx);
        int maxWidth = 0;
        net.totalTiles = 0;
        for (const auto &row : net.rows) {
            maxWidth = max(maxWidth, static_cast<int>(row.size()));
        }
        for (auto &row : net.rows) {
            if (static_cast<int>(row.size()) < maxWidth) {
                row += string(maxWidth - row.size(), '.');
            }
            net.totalTiles += static_cast<int>(count(row.begin(), row.end(), '#'));
        }
    }
    return nets;
}

// Apply the given global marker set (1-based coordinates) to all nets.
void applyMarkersToNets(vector<Net> &nets) {
    vector<tuple<int,int,char>> markerData = {
        {1, 15, 'c'},
        {2, 7,  's'},
        {3, 14, 'c'},
        {4, 10, 'c'},
        {4, 17, 'c'},
        {5, 8,  's'},
        {7, 17, 's'},
        {8, 2,  'c'},
        {8, 11, 's'},
        {10,18, 's'},
        {13,3,  's'},
        {17,12, 'c'},
        {18,8,  's'}
    };

    for (auto &net : nets) {
        net.markers.clear();
        for (auto [r1, c1, t] : markerData) {
            int r = r1 - 1;
            int c = c1 - 1;
            if (r < 0 || r >= (int)net.rows.size()) continue;
            if (c < 0 || c >= (int)net.rows[r].size()) continue;
            if (net.rows[r][c] != '#') continue;
            Marker m;
            m.row = r;
            m.col = c;
            m.type = t;
            net.markers.push_back(m);
        }
    }
}

ProcessResult processNet(const Net &net) {
    ProcessResult res;
    res.header = net.header;
    res.index = net.index;
    res.totalTiles = net.totalTiles;
    res.originalGrid = net.rows;

    if (net.rows.empty()) {
        res.failureReason = "No grid data";
        return res;
    }
    if (net.totalTiles == 0) {
        res.failureReason = "Grid contains no '#'";
        return res;
    }
    if (net.totalTiles % 2 != 0) {
        res.failureReason = "Odd number of tiles";
        return res;
    }

    const int halfArea = net.totalTiles / 2;
    vector<tuple<int,int,int>> triples = enumerateTriples(halfArea);
    if (triples.empty()) {
        res.failureReason = "No integer box dimensions for tile count";
        return res;
    }

    FoldingSolver solver(net.rows, net.index);

    if (!net.markers.empty()) {
        vector<tuple<int,int,char>> marks;
        marks.reserve(net.markers.size());
        for (const Marker &m : net.markers) {
            marks.emplace_back(m.row, m.col, m.type);
        }
        solver.setMarkers(marks);
    }

    if (solver.tileCount() != net.totalTiles) {
        res.failureReason = "Internal tile count mismatch";
        return res;
    }

    long long totalIterations = 0;

    for (auto [a, b, c] : triples) {
        if (g_stop.load(memory_order_relaxed)) break;

        vector<int> assignment;
        int rootPatch = -1;

        bool ok = solver.solveForDims(a, b, c, assignment, rootPatch);
        totalIterations += solver.getIterationCount();

        if (ok) {
            res.valid = true;
            res.Lx = a;
            res.Ly = b;
            res.Lz = c;
            res.rootPatch = rootPatch;
            res.labeledGrid = solver.buildLabeledGrid(assignment, a, b, c);
            res.mappingLines = solver.describeAssignment(assignment, a, b, c);
            res.dfsIterations = totalIterations;

            {
                lock_guard<mutex> lock(g_cerrMutex);
                cerr << "\r[Net " << (net.index + 1) << "] "
                     << "valid folding for " << a << "x" << b << "x" << c
                     << " after " << totalIterations << " DFS calls\n";
            }

            return res;
        }
    }

    res.failureReason = "No folding matched any candidate box dimensions";
    res.dfsIterations = totalIterations;

    {
        lock_guard<mutex> lock(g_cerrMutex);
        cerr << "\r[Net " << (net.index + 1) << "] "
             << "no valid folding after " << totalIterations
             << " DFS calls\n";
    }

    return res;
}

void writeWinnerResult(const vector<ProcessResult> &results,
                       const string &path,
                       int winnerIdx) {
    ofstream out(path);
    if (!out) {
        lock_guard<mutex> lock(g_cerrMutex);
        cerr << "Failed to open " << path << " for writing.\n";
        return;
    }

    const string legend = "1=-X  2=+X  3=-Y  4=+Y  5=-Z  6=+Z";

    if (winnerIdx < 0 || winnerIdx >= (int)results.size() ||
        !results[winnerIdx].valid) {
        out << "No valid embedding found.\n";
        return;
    }

    const auto &res = results[winnerIdx];

    out << res.header << "\n";
    out << "Dimensions: " << res.Lx
        << " x " << res.Ly
        << " x " << res.Lz << "\n";
    out << "Root patch index: " << res.rootPatch << "\n";
    out << legend << "\n";
    out << "DFS iterations: " << res.dfsIterations << "\n\n";

    out << "Labeled grid:\n";
    for (const auto &row : res.labeledGrid) out << row << "\n";

    out << "\nOriginal grid:\n";
    for (const auto &row : res.originalGrid) out << row << "\n";

    out << "\nTile mapping:\n";
    for (const auto &line : res.mappingLines) out << line << "\n";

    out << "\n====\n";
}

} // namespace

int main() {
    const string inputPath = "solutions.txt";
    const string outputPath = "solutionFINAL.txt";

    vector<Net> nets = loadNets(inputPath);
    if (nets.empty()) {
        lock_guard<mutex> lock(g_cerrMutex);
        cerr << "No nets available to process.\n";
        return 1;
    }

    applyMarkersToNets(nets);

    for (int i = 0; i < static_cast<int>(nets.size()); ++i) {
        nets[i].index = i;
    }

    vector<ProcessResult> results(nets.size());
    atomic<int> nextIndex{0};
    const int totalNets = static_cast<int>(nets.size());
    const int threadCount = min(8, totalNets);
    vector<thread> workers;
    workers.reserve(threadCount);

    auto worker = [&]() {
        while (true) {
            if (g_stop.load(memory_order_relaxed)) break;
            int idx = nextIndex.fetch_add(1);
            if (idx >= totalNets) break;

            const Net &net = nets[idx];
            {
                lock_guard<mutex> lock(g_cerrMutex);
                cerr << "\n[Net " << (net.index + 1) << "/" << totalNets << "] "
                     << net.header << "\n";
            }

            ProcessResult res = processNet(net);
            results[idx] = res;

            if (res.valid) {
                int expected = -1;
                if (g_winnerIndex.compare_exchange_strong(expected, idx)) {
                    g_stop.store(true, memory_order_relaxed);
                }
            }
        }
    };

    for (int i = 0; i < threadCount; ++i) {
        workers.emplace_back(worker);
    }
    for (auto &t : workers) t.join();

    int winnerIdx = g_winnerIndex.load();
    writeWinnerResult(results, outputPath, winnerIdx);

    {
        lock_guard<mutex> lock(g_cerrMutex);
        if (winnerIdx >= 0) {
            cerr << "\nFirst valid embedding found in net index "
                 << winnerIdx << " (1-based: " << (winnerIdx + 1) << ")\n";
        } else {
            cerr << "\nNo valid embedding found in any net.\n";
        }
        cerr << "Solution written to " << outputPath << "\n";
    }

    return 0;
}
