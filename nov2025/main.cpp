#include <iostream>
#include <iomanip>
#include <vector>
#include <queue>
#include <fstream>
#include <thread>
#include <atomic>
#include <mutex>
#include <array>
#include <algorithm>

using namespace std;

const int N = 20;

// Per-thread grid type: -1 = unknown, 0 = not in box, 1 = in box
using Grid = array<array<int, N>, N>;

atomic<long long> nodesVisited{0};   // total dfs calls across threads
const long long REPORT_EVERY = 1000000;

enum {UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3};
const int DR[4] = {-1, 1, 0, 0};
const int DC[4] = {0, 0, -1, 1};
const int DIR_BIT[4] = {1, 2, 4, 8}; // U,D,L,R bits

struct Arrow {
    int r, c;
    int mask;                      // which directions are shown
    vector<pair<int,int>> ray[4];  // cells in each direction, outward from arrow
};

struct NumberClue {
    int r, c;
    int val;
    vector<pair<int,int>> neigh;   // cells within one king move (including itself)
};

vector<Arrow> arrows;
vector<NumberClue> numbersClue;
vector<pair<int,int>> vars;

// solution output (shared)
ofstream out;
atomic<int> solutionsFound{0};
mutex outMutex;  // protects out and solutionsFound

bool inside(int r,int c){ return r>=0 && r<N && c>=0 && c<N; }

// ---------- Number constraints (now take Grid) ----------

// partial check: every clue still satisfiable
bool checkNumbersPartial(const Grid &state) {
    for (auto &cl : numbersClue) {
        int sum1 = 0, unknown = 0;
        for (auto [r,c] : cl.neigh) {
            if (state[r][c] == 1) sum1++;
            else if (state[r][c] == -1) unknown++;
        }
        if (sum1 > cl.val) return false;
        if (sum1 + unknown < cl.val) return false;
    }
    return true;
}

// final check: exact equality
bool checkNumbersFinal(const Grid &state) {
    for (auto &cl : numbersClue) {
        int sum1 = 0;
        for (auto [r,c] : cl.neigh)
            if (state[r][c] == 1) sum1++;
        if (sum1 != cl.val) return false;
    }
    return true;
}

// ---------- Arrow constraints ----------

// light pruning: each shown direction must still have some cell that could be 1
bool checkArrowsPartial(const Grid &state) {
    for (auto &ar : arrows) {
        for (int d = 0; d < 4; ++d) {
            if (ar.mask & DIR_BIT[d]) {
                bool possible = false;
                for (auto [r,c] : ar.ray[d]) {
                    if (state[r][c] != 0) { // 1 or unknown
                        possible = true;
                        break;
                    }
                }
                if (!possible) return false; // this direction is dead
            }
        }
    }
    return true;
}

// final arrow rule: arrows point to all directions whose nearest box cell
// is at minimum distance (in that row/column)
bool checkArrowsFinal(const Grid &state) {
    const int INF = 1e9;
    for (auto &ar : arrows) {
        int best[4];
        for (int d = 0; d < 4; ++d) {
            best[d] = INF;
            int dist = 1;
            for (auto [r,c] : ar.ray[d]) {
                if (state[r][c] == 1) {
                    best[d] = dist;
                    break;
                }
                dist++;
            }
        }
        int m = INF;
        for (int d = 0; d < 4; ++d) m = min(m, best[d]);
        if (m == INF) return false; // no box cell in any direction

        int mask = 0;
        for (int d = 0; d < 4; ++d)
            if (best[d] == m) mask |= DIR_BIT[d];

        if (mask != ar.mask) return false;
    }
    return true;
}

// ---------- connectivity + holes ----------

bool checkConnectivityAndNoHoles(const Grid &state) {
    int total1 = 0;
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            if (state[r][c] == 1) total1++;

    if (total1 == 0) return false;

    // connectivity of 1-cells
    int vis[N][N] = {};
    queue<pair<int,int>> q;

    bool started = false;
    for (int r = 0; r < N && !started; ++r)
        for (int c = 0; c < N && !started; ++c)
            if (state[r][c] == 1) {
                started = true;
                vis[r][c] = 1;
                q.push({r,c});
            }

    int seen = 0;
    while (!q.empty()) {
        auto [r,c] = q.front(); q.pop();
        seen++;
        for (int d = 0; d < 4; ++d) {
            int nr = r + DR[d], nc = c + DC[d];
            if (inside(nr,nc) && !vis[nr][nc] && state[nr][nc] == 1) {
                vis[nr][nc] = 1;
                q.push({nr,nc});
            }
        }
    }
    if (seen != total1) return false; // disconnected

    // hole check for 0-cells: every 0 must connect (4-neighbour) to boundary 0s
    int vis0[N][N] = {};
    queue<pair<int,int>> q0;

    // start BFS from all boundary 0-cells
    for (int r = 0; r < N; ++r) {
        int c = 0;
        if (state[r][c] == 0 && !vis0[r][c]) {
            vis0[r][c] = 1; q0.push({r,c});
        }
        c = N-1;
        if (state[r][c] == 0 && !vis0[r][c]) {
            vis0[r][c] = 1; q0.push({r,c});
        }
    }
    for (int c = 0; c < N; ++c) {
        int r = 0;
        if (state[r][c] == 0 && !vis0[r][c]) {
            vis0[r][c] = 1; q0.push({r,c});
        }
        r = N-1;
        if (state[r][c] == 0 && !vis0[r][c]) {
            vis0[r][c] = 1; q0.push({r,c});
        }
    }

    while (!q0.empty()) {
        auto [r,c] = q0.front(); q0.pop();
        for (int d = 0; d < 4; ++d) {
            int nr = r + DR[d], nc = c + DC[d];
            if (inside(nr,nc) && !vis0[nr][nc] && state[nr][nc] == 0) {
                vis0[nr][nc] = 1;
                q0.push({nr,nc});
            }
        }
    }

    // any 0 not reached from boundary is an interior "hole"
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            if (state[r][c] == 0 && !vis0[r][c])
                return false;

    return true;
}

// ---------- output ----------

void printSolution(const Grid &state, ostream &os) {
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c)
            os << (state[r][c] == 1 ? '#' : '.');
        os << "\n";
    }
}

void reportProgress(int depth) {
    long long nv = nodesVisited.load(memory_order_relaxed);
    if (nv % REPORT_EVERY == 0) {
        cerr << "Visited " << nv
             << " nodes, depth " << depth
             << " / " << vars.size()
             << ", solutions " << solutionsFound.load() << "\r";
    }
}

// ---------- DFS search ----------

void dfs(int idx, Grid &state) {
    ++nodesVisited;
    reportProgress(idx);

    if (idx == (int)vars.size()) {
        // final checks
        if (!checkNumbersFinal(state)) return;
        if (!checkArrowsFinal(state)) return;
        if (!checkConnectivityAndNoHoles(state)) return;

        // additional condition: even number of box cells
        int total1 = 0;
        for (int r = 0; r < N; ++r)
            for (int c = 0; c < N; ++c)
                if (state[r][c] == 1) ++total1;
        if (total1 % 2 != 0) return;

        // write solution safely
        lock_guard<mutex> lock(outMutex);
        int myIndex = ++solutionsFound;
        out << "Solution #" << myIndex
            << " (total box cells = " << total1 << ")\n";
        printSolution(state, out);
        out << "\n";
        return;
    }

    auto [r,c] = vars[idx];

    for (int v = 0; v <= 1; ++v) {
        state[r][c] = v;
        if (!checkNumbersPartial(state)) { state[r][c] = -1; continue; }
        if (!checkArrowsPartial(state))  { state[r][c] = -1; continue; }
        dfs(idx+1, state);
        state[r][c] = -1;
    }
}

// ---------- puzzle set-up ----------

void addArrow(Grid &state, int r,int c,int mask) {
    Arrow ar;
    ar.r = r; ar.c = c; ar.mask = mask;
    for (int d = 0; d < 4; ++d) {
        int rr = r + DR[d], cc = c + DC[d];
        while (inside(rr,cc)) {
            ar.ray[d].push_back({rr,cc});
            rr += DR[d]; cc += DC[d];
        }
    }
    arrows.push_back(ar);
    state[r][c] = 0; // arrows are not in the box
}

void addNumber(Grid &state, int r,int c,int val) {
    NumberClue cl;
    cl.r = r; cl.c = c; cl.val = val;
    for (int dr=-1; dr<=1; ++dr)
        for (int dc=-1; dc<=1; ++dc) {
            int rr = r+dr, cc = c+dc;
            if (inside(rr,cc))
                cl.neigh.push_back({rr,cc});
        }
    numbersClue.push_back(cl);
    state[r][c] = 1; // numbered cells are always in the box
}

struct ArrowSeed {
    int r, c;     // 0-based row, col
    int mask;     // bitmask using DIR_BIT[UP/DOWN/LEFT/RIGHT]
};

struct NumberSeed {
    int r, c;
    int val;
};

struct FixedSeed {
    int r, c;
    int val;      // 0 = forced not in box, 1 = forced in box
};

void initPuzzle(Grid &state) {
    // start all cells as unknown
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            state[r][c] = -1;

    vector<NumberSeed> numberSeeds = {
        {1,11,4}, {1,15,4}, {2,7,5}, {3,12,7}, {3,14,5}, {4,10,4}, {4,13,7}, {4,17,4},
        {5,6,4}, {5,8,7}, {6,7,9}, {7,17,6}, {8,2,7}, {8,11,5}, {9,14,5}, {10,5,4},
        {10,7,7}, {10,18,3}, {13,3,5}, {13,6,6}, {13,9,2}, {15,6,5}, {17,12,5}, {17,13,5},
        {18,8,4}
    };
    vector<ArrowSeed> arrowSeeds = {
        {0,1,DIR_BIT[RIGHT]}, {0,8,DIR_BIT[DOWN] | DIR_BIT[LEFT]}, {0,12,DIR_BIT[DOWN]},
        {0,15,DIR_BIT[DOWN]}, {0,19, DIR_BIT[LEFT] | DIR_BIT[DOWN]}, {1,5,DIR_BIT[RIGHT]},
        {1,13,DIR_BIT[LEFT] | DIR_BIT[DOWN] | DIR_BIT[RIGHT]}, {2,6,DIR_BIT[RIGHT]}, {2,9,DIR_BIT[LEFT] | DIR_BIT[DOWN] | DIR_BIT[RIGHT]},
        {2,18,DIR_BIT[UP] | DIR_BIT[LEFT] | DIR_BIT[DOWN]}, {3,1,DIR_BIT[DOWN]}, {3,3,DIR_BIT[DOWN] | DIR_BIT[RIGHT]}, {4,0,DIR_BIT[RIGHT]},
        {4,5,DIR_BIT[DOWN] | DIR_BIT[RIGHT]}, {4,15,DIR_BIT[LEFT] | DIR_BIT[UP] | DIR_BIT[RIGHT]}, {5,11,DIR_BIT[DOWN] | DIR_BIT[RIGHT]},
        {6,3,DIR_BIT[DOWN]}, {6,10,DIR_BIT[LEFT] | DIR_BIT[RIGHT]}, {6,13,DIR_BIT[UP] | DIR_BIT[LEFT] | DIR_BIT[DOWN]}, {6,16,DIR_BIT[RIGHT]},
        {6,19,DIR_BIT[LEFT] | DIR_BIT[DOWN]}, {7,1,DIR_BIT[DOWN] | DIR_BIT[RIGHT]}, {7,14,DIR_BIT[LEFT] | DIR_BIT[DOWN]}, {8,5,DIR_BIT[UP] | DIR_BIT[RIGHT] | DIR_BIT[DOWN]},
        {8,9,DIR_BIT[UP] | DIR_BIT[LEFT] | DIR_BIT[DOWN] | DIR_BIT[RIGHT]}, {9,1,DIR_BIT[UP] | DIR_BIT[RIGHT]}, {9,12,DIR_BIT[UP] | DIR_BIT[LEFT] | DIR_BIT[RIGHT]},
        {9,16,DIR_BIT[RIGHT] | DIR_BIT[DOWN]}, {9,18,DIR_BIT[UP] | DIR_BIT[LEFT] | DIR_BIT[DOWN]}, {10,1,DIR_BIT[UP]}, {10,3,DIR_BIT[UP] | DIR_BIT[DOWN]},
        {11,8,DIR_BIT[UP] | DIR_BIT[LEFT] | DIR_BIT[RIGHT]}, {11,10,DIR_BIT[RIGHT] | DIR_BIT[LEFT]}, {11,14,DIR_BIT[UP] | DIR_BIT[LEFT] | DIR_BIT[RIGHT]},
        {11,17,DIR_BIT[UP] | DIR_BIT[LEFT]}, {12,2,DIR_BIT[UP] | DIR_BIT[RIGHT]}, {12,5,DIR_BIT[UP] | DIR_BIT[LEFT] | DIR_BIT[DOWN] | DIR_BIT[RIGHT]},
        {12,18,DIR_BIT[UP]}, {13,0,DIR_BIT[RIGHT]}, {13,12,DIR_BIT[UP] | DIR_BIT[LEFT]}, {13,16,DIR_BIT[UP]}, {14,8,DIR_BIT[UP] | DIR_BIT[LEFT] | DIR_BIT[DOWN]},
        {14,11,DIR_BIT[UP] | DIR_BIT[DOWN]}, {14,13,DIR_BIT[DOWN]}, {15,2,DIR_BIT[UP] | DIR_BIT[RIGHT]}, {15,4,DIR_BIT[UP]}, {15,9,DIR_BIT[LEFT] | DIR_BIT[DOWN]},
        {15,14,DIR_BIT[LEFT] | DIR_BIT[DOWN]}, {15,19,DIR_BIT[LEFT]}, {16,5,DIR_BIT[UP]}, {16,7,DIR_BIT[UP] | DIR_BIT[RIGHT]}, {16,16,DIR_BIT[LEFT]}, {16,18,DIR_BIT[UP] | DIR_BIT[LEFT]},
        {17,1,DIR_BIT[RIGHT]}, {17,10,DIR_BIT[UP] | DIR_BIT[LEFT] | DIR_BIT[DOWN]}, {18,4,DIR_BIT[UP] | DIR_BIT[RIGHT]}, {18,6,DIR_BIT[RIGHT]}, {18,14,DIR_BIT[UP] | DIR_BIT[LEFT]},
        {19,0,DIR_BIT[RIGHT]}, {19,4,DIR_BIT[UP] | DIR_BIT[RIGHT]}, {19,7,DIR_BIT[RIGHT]}, {19,11,DIR_BIT[LEFT]}, {19,18,DIR_BIT[UP] | DIR_BIT[LEFT]}
    };

    vector<FixedSeed> fixedSeeds = {
        {0,0,0}, {0,2,0}, {0,3,0}, {0,4,0}, {0,7,1}, {0,9,0}, {0,10,0}, {0,11,0}, {0,13,0}, {0,14,0}, {0,16,0}, {0,18,1}, {1,0,0}, {1,1,0}, {1,2,0}, {1,3,0}, {1,6,0}, {1,7,1},
        {1,8,1}, {1,9,0}, {1,10,1}, {1,12,0}, {1,14,0}, {1,17,1}, {1,18,1}, {1,19,1}, {2,0,0}, {2,1,0}, {2,2,0}, {2,3,0}, {2,4,0}, {2,5,0}, {2,8,0}, {2,10,0}, {2,11,1}, {2,12,1}, {2,13,0},
        {2,14,1}, {2,15,1}, {2,17,1}, {2,19,0}, {3,0,0}, {3,2,0}, {3,4,0}, {3,5,0}, {3,6,0}, {3,7,1}, {3,8,1}, {3,9,0}, {3,11,1}, {3,13,1}, {3,15,0}, {3,16,0}, {3,17,1}, {3,18,1},
        {4,1,0}, {4,2,0}, {4,3,0}, {4,4,0}, {4,6,0}, {4,7,0}, {4,8,1}, {4,9,1}, {4,11,0}, {4,12,1}, {4,14,0}, {4,16,0}, {4,18,0}, {5,0,0}, {5,1,0}, {5,2,0}, {5,3,0}, {5,4,0}, {5,5,0}, {5,7,1},
        {5,10,0}, {5,15,0}, {5,16,0}, {5,17,1}, {5,18,0}, {5,19,0}, {6,0,0}, {6,1,0}, {6,2,0}, {6,4,0}, {6,5,0}, {6,6,1}, {6,8,1}, {6,14,0}, {6,15,0}, {6,17,1}, {6,18,1}, {7,0,0}, {7,2,1},
        {7,3,1}, {7,5,1}, {7,6,1}, {7,7,1}, {7,8,1}, {7,10,0}, {7,15,0}, {7,16,0}, {7,18,1}, {7,19,1}, {8,0,0}, {8,1,1}, {8,3,1}, {8,4,0}, {8,6,1}, {8,12,1}, {8,13,1}, {8,16,0}, {8,17,1}, {8,18,1}, {8,19,0},
        {9,0,0}, {9,2,1}, {9,3,1}, {9,5,1}, {9,11,1}, {9,13,1}, {9,15,0}, {9,17,1}, {9,19,0}, {10,0,0}, {10,2,0}, {10,4,0}, {10,8,1}, {10,10,0}, {10,11,1}, {10,12,0}, {10,14,0}, {10,16,1}, {10,17,1},
        {10,19,0}, {11,0,0}, {11,1,0}, {11,2,1}, {11,3,1}, {11,5,0}, {11,7,1}, {11,9,1}, {11,11,1}, {11,12,1}, {11,13,0}, {11,15,0}, {11,16,1}, {11,18,0}, {11,19,0}, {12,0,0},
        {12,1,0}, {12,3,1}, {12,4,0}, {12,6,0}, {12,7,1}, {12,8,0}, {12,9,0}, {12,10,0}, {12,11,1}, {12,12,1}, {12,13,0}, {12,14,0}, {12,15,0}, {12,16,0}, {12,17,0}, {12,19,0}, {13,1,0}, {13,2,0},
        {13,4,1}, {13,5,0}, {13,7,1}, {13,8,1}, {13,10,0}, {13,11,1}, {13,13,0}, {13,14,0}, {13,15,0}, {13,17,0}, {13,18,0}, {13,19,0}, {14,0,0}, {14,1,0}, {14,2,0}, {14,3,1}, {14,4,1}, {14,5,1},
        {14,6,1}, {14,7,1}, {14,9,0}, {14,10,0}, {14,12,0}, {14,14,0}, {14,15,0}, {14,16,0}, {14,17,0}, {14,18,0}, {14,19,0}, {15,0,0}, {15,1,0}, {15,3,0}, {15,5,0}, {15,7,1}, {15,8,1}, {15,10,0},
        {15,11,1}, {15,12,1}, {15,13,0}, {15,15,0}, {15,16,0}, {15,17,0}, {15,18,0}, {16,0,0}, {16,1,0}, {16,2,0}, {16,3,0}, {16,4,0}, {16,6,0}, {16,8,1}, {16,9,1}, {16,10,1}, {16,11,1}, {16,12,1},
        {16,13,0}, {16,14,0}, {16,15,0}, {16,17,0}, {16,19,0}, {17,0,0}, {17,2,0}, {17,3,0}, {17,4,0}, {17,5,0}, {17,6,0}, {17,7,0}, {17,8,0}, {17,9,1}, {17,11,0}, {17,14,1}, {17,16,0}, {17,17,0}, {17,18,0},
        {17,19,0}, {18,0,0}, {18,1,0}, {18,2,0}, {18,3,0}, {18,5,0}, {18,7,0}, {18,9,1}, {18,10,1}, {18,11,0}, {18,12,0}, {18,13,1}, {18,15,0}, {18,16,0}, {18,17,0}, {18,18,0}, {18,19,0}, {19,1,0}, {19,2,0},
        {19,3,0}, {19,5,0}, {19,6,0}, {19,8,0}, {19,9,1}, {19,10,0}, {19,12,0}, {19,13,0}, {19,14,0}, {19,15,0}, {19,16,0}, {19,17,0}, {19,19,0}
    };

    // numbered cells
    for (auto &ns : numberSeeds) {
        addNumber(state, ns.r, ns.c, ns.val);
    }

    // arrow cells
    for (auto &as : arrowSeeds) {
        addArrow(state, as.r, as.c, as.mask);
    }

    // manual fixed cells
    for (auto &fs : fixedSeeds) {
        state[fs.r][fs.c] = fs.val;
    }

    // collect variables
    vars.clear();
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            if (state[r][c] == -1)
                vars.push_back({r,c});
}

int main() {
    Grid base;
    initPuzzle(base);

    // optional sanity check
    if (!checkNumbersPartial(base)) {
        cerr << "Initial fixed cells already violate a number clue.\n";
        return 0;
    }
    if (!checkArrowsPartial(base)) {
        cerr << "Initial fixed cells already violate an arrow clue.\n";
        return 0;
    }

    cerr << "Variables to search: " << vars.size() << "\n";

    out.open("solutions.txt");
    if (!out) {
        cerr << "Failed to open solutions.txt for writing.\n";
        return 1;
    }

    // If very few vars, just run single-threaded to keep it simple
    if (vars.size() < 3) {
        Grid st = base;
        dfs(0, st);
    } else {
        // 3 variables → 2^3 = 8 threads
        auto [r0, c0] = vars[0];
        auto [r1, c1] = vars[1];
        auto [r2, c2] = vars[2];

        auto launch3 = [&](int v0, int v1, int v2) {
            return std::thread([&, v0, v1, v2] {
                Grid st = base;
                st[r0][c0] = v0;
                st[r1][c1] = v1;
                st[r2][c2] = v2;

                if (!checkNumbersPartial(st) || !checkArrowsPartial(st))
                    return;

                dfs(3, st); // already assigned 3 vars
            });
        };

        std::thread t000 = launch3(0, 0, 0);
        std::thread t001 = launch3(0, 0, 1);
        std::thread t010 = launch3(0, 1, 0);
        std::thread t011 = launch3(0, 1, 1);
        std::thread t100 = launch3(1, 0, 0);
        std::thread t101 = launch3(1, 0, 1);
        std::thread t110 = launch3(1, 1, 0);
        std::thread t111 = launch3(1, 1, 1);

        t000.join(); t001.join(); t010.join(); t011.join();
        t100.join(); t101.join(); t110.join(); t111.join();
    }

    out.close();
    cerr << "\nFinished search. Solutions found: " << solutionsFound.load() << "\n";

    if (solutionsFound == 0)
        cout << "No pattern satisfies the constraints.\n";
    else
        cout << "Solutions written to solutions.txt, count = "
             << solutionsFound.load() << "\n";

    return 0;
}
