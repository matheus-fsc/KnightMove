// knight_transfer.cpp — Broken profile DP otimizado.
//
// State bit-packed (4 bits/cell em 2 uint64 = 32 cells max).
// Vector sort+dedup em vez de unordered_map (memory-friendlier).
//
// Compilar: g++ -O3 -std=c++17 -o knight_transfer knight_transfer.cpp

#include <iostream>
#include <vector>
#include <cstdint>
#include <chrono>
#include <algorithm>
#include <cstring>

using namespace std;

int N;
int V;

const int KDR[8] = {-2, -2, -1, -1, 1, 1, 2, 2};
const int KDC[8] = {-1, 1, -2, 2, -2, 2, -1, 1};

vector<vector<int>> past_nbr;
vector<int> max_fut;
vector<vector<int>> frontier_of_v;
vector<vector<int>> profile_idx_of_v;

void precompute() {
    past_nbr.assign(V, {});
    max_fut.assign(V, -1);
    for (int v = 0; v < V; v++) {
        int r = v / N, c = v % N;
        for (int k = 0; k < 8; k++) {
            int rr = r + KDR[k], cc = c + KDC[k];
            if (0 <= rr && rr < N && 0 <= cc && cc < N) {
                int u = rr * N + cc;
                if (u < v) past_nbr[v].push_back(u);
                else if (u > v) {
                    if (u > max_fut[v]) max_fut[v] = u;
                }
            }
        }
    }
    frontier_of_v.assign(V + 1, {});
    profile_idx_of_v.assign(V + 1, vector<int>(V, -1));
    for (int v = 0; v <= V; v++) {
        for (int u = 0; u < v; u++) {
            if (max_fut[u] >= v) {
                frontier_of_v[v].push_back(u);
            }
        }
        for (size_t i = 0; i < frontier_of_v[v].size(); i++) {
            profile_idx_of_v[v][frontier_of_v[v][i]] = (int)i;
        }
    }
}

struct PState {
    uint64_t lo, hi;
    bool operator==(const PState& o) const { return lo == o.lo && hi == o.hi; }
    bool operator<(const PState& o) const {
        return hi != o.hi ? hi < o.hi : lo < o.lo;
    }
};

inline int get_cell(const PState& s, int i) {
    if (i < 16) return (int)((s.lo >> (4 * i)) & 0xFULL);
    return (int)((s.hi >> (4 * (i - 16))) & 0xFULL);
}

inline void set_cell(PState& s, int i, int val) {
    if (i < 16) {
        s.lo &= ~(0xFULL << (4 * i));
        s.lo |= (uint64_t)(val & 0xF) << (4 * i);
    } else {
        s.hi &= ~(0xFULL << (4 * (i - 16)));
        s.hi |= (uint64_t)(val & 0xF) << (4 * (i - 16));
    }
}

inline void canonicalize(PState& s, int len) {
    int next_id = 1;
    int remap[16];
    memset(remap, 0, sizeof(remap));
    for (int i = 0; i < len; i++) {
        int c = get_cell(s, i);
        if (c == 0 || c == 15) continue;
        if (remap[c] == 0) remap[c] = next_id++;
        if (remap[c] != c) set_cell(s, i, remap[c]);
    }
}

inline int find_unused_plug_bp(const PState& s, int len) {
    int used = 0;
    for (int i = 0; i < len; i++) {
        int c = get_cell(s, i);
        if (c >= 1 && c <= 14) used |= (1 << c);
    }
    for (int i = 1; i <= 14; i++) {
        if (!(used & (1 << i))) return i;
    }
    return -1;
}

inline void replace_id(PState& s, int len, int old_id, int new_id) {
    for (int i = 0; i < len; i++) {
        if (get_cell(s, i) == old_id) set_cell(s, i, new_id);
    }
}

// Para n=8, total ≈ 1.3e13, cabe em uint64.
// Para n=10+, usar __int128.
using Count = uint64_t;

string to_string128(uint64_t x) {
    return to_string(x);
}

struct Entry {
    PState st;
    Count cnt;
};

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <N>" << endl;
        return 1;
    }
    N = atoi(argv[1]);
    V = N * N;
    if (N < 1 || N > 14) {
        cerr << "N must be in [1, 14]" << endl;
        return 1;
    }

    precompute();

    int max_front = 0;
    for (int v = 0; v <= V; v++)
        max_front = max(max_front, (int)frontier_of_v[v].size());
    cerr << "N=" << N << ", V=" << V << ", max frontier=" << max_front << endl;
    if (max_front > 32) {
        cerr << "ERROR: max frontier > 32, state doesn't fit in 128 bits" << endl;
        return 1;
    }

    auto start = chrono::steady_clock::now();

    // dp_current/dp_next como vector ordenado, sem duplicatas
    vector<Entry> dp_current;
    dp_current.push_back({{0, 0}, 1});

    vector<Entry> raw_next;  // entries antes de sort/merge

    for (int v = 0; v < V; v++) {
        const auto& old_front = frontier_of_v[v];
        const auto& new_front = frontier_of_v[v + 1];
        int old_len = (int)old_front.size();
        int new_len = (int)new_front.size();
        bool v_in_new_front = (max_fut[v] > v);

        vector<int> S_choices;
        for (int u : past_nbr[v]) {
            int idx = profile_idx_of_v[v][u];
            S_choices.push_back(idx);
        }
        int K = (int)S_choices.size();

        vector<int> old_to_new(old_front.size(), -1);
        for (size_t i = 0; i < old_front.size(); i++) {
            int u = old_front[i];
            for (size_t j = 0; j < new_front.size(); j++) {
                if (new_front[j] == u) { old_to_new[i] = (int)j; break; }
            }
        }
        int v_new_idx = -1;
        if (v_in_new_front) {
            for (size_t j = 0; j < new_front.size(); j++) {
                if (new_front[j] == v) { v_new_idx = (int)j; break; }
            }
        }
        bool is_last = (v == V - 1);
        int working_len = old_len + 1;
        int v_idx_in_w = old_len;

        // Pré-alocar dp_next baseado em estimativa
        raw_next.clear();
        raw_next.reserve(dp_current.size() * 4);

        for (const Entry& e : dp_current) {
            const PState& state = e.st;
            Count cnt = e.cnt;

            for (int mask = 0; mask < (1 << K); mask++) {
                int popcnt = __builtin_popcount((unsigned)mask);
                if (popcnt > 2) continue;

                PState working = state;
                set_cell(working, v_idx_in_w, 0);

                bool valid = true;

                for (int b = 0; b < K; b++) {
                    if (!(mask & (1 << b))) continue;
                    int u_idx = S_choices[b];
                    int lu = get_cell(working, u_idx);
                    int lv = get_cell(working, v_idx_in_w);

                    if (lu == 15 || lv == 15) { valid = false; break; }

                    if (lu == 0 && lv == 0) {
                        int np = find_unused_plug_bp(working, working_len);
                        if (np < 0) { valid = false; break; }
                        set_cell(working, u_idx, np);
                        set_cell(working, v_idx_in_w, np);
                    } else if (lu == 0) {
                        set_cell(working, u_idx, lv);
                        set_cell(working, v_idx_in_w, 15);
                    } else if (lv == 0) {
                        set_cell(working, u_idx, 15);
                        set_cell(working, v_idx_in_w, lu);
                    } else if (lu == lv) {
                        if (!is_last) { valid = false; break; }
                        set_cell(working, u_idx, 15);
                        set_cell(working, v_idx_in_w, 15);
                    } else {
                        set_cell(working, u_idx, 15);
                        set_cell(working, v_idx_in_w, 15);
                        replace_id(working, working_len, lv, lu);
                    }
                }
                if (!valid) continue;

                PState new_state{0, 0};
                bool drop_valid = true;
                for (int i = 0; i < old_len; i++) {
                    int wval = get_cell(working, i);
                    if (old_to_new[i] == -1) {
                        if (wval != 15) { drop_valid = false; break; }
                    } else {
                        set_cell(new_state, old_to_new[i], wval);
                    }
                }
                if (!drop_valid) continue;

                if (v_in_new_front) {
                    set_cell(new_state, v_new_idx, get_cell(working, v_idx_in_w));
                } else {
                    if (get_cell(working, v_idx_in_w) != 15) continue;
                }

                int plug_count[16] = {0};
                for (int i = 0; i < new_len; i++) {
                    int c = get_cell(new_state, i);
                    if (c >= 1 && c <= 14) plug_count[c]++;
                }
                bool plugs_ok = true;
                for (int i = 1; i <= 14; i++) {
                    if (plug_count[i] != 0 && plug_count[i] != 2) {
                        plugs_ok = false; break;
                    }
                }
                if (!plugs_ok) continue;

                canonicalize(new_state, new_len);

                raw_next.push_back({new_state, cnt});
            }
        }

        // Sort + merge duplicates
        sort(raw_next.begin(), raw_next.end(),
             [](const Entry& a, const Entry& b) { return a.st < b.st; });

        dp_current.clear();
        if (!raw_next.empty()) {
            dp_current.reserve(raw_next.size());
            dp_current.push_back(raw_next[0]);
            for (size_t i = 1; i < raw_next.size(); i++) {
                if (raw_next[i].st == dp_current.back().st) {
                    dp_current.back().cnt += raw_next[i].cnt;
                } else {
                    dp_current.push_back(raw_next[i]);
                }
            }
        }

        auto elapsed = chrono::duration_cast<chrono::milliseconds>(
            chrono::steady_clock::now() - start).count();
        cerr << "v=" << v
             << " front=" << old_len << "→" << new_len
             << " states=" << dp_current.size()
             << " raw=" << raw_next.size()
             << " t=" << elapsed/1000.0 << "s" << endl;
    }

    Count result = 0;
    PState empty{0, 0};
    for (const Entry& e : dp_current) {
        if (e.st == empty) { result = e.cnt; break; }
    }

    auto elapsed = chrono::duration_cast<chrono::milliseconds>(
        chrono::steady_clock::now() - start).count();
    cout << "N(" << N << "x" << N << ") = " << to_string128(result) << endl;
    cout << "Time: " << elapsed/1000.0 << "s" << endl;
    cout << "N mod 8 = " << to_string128(result % 8) << endl;
    cout << "N / 8 = " << to_string128(result / 8) << endl;

    return 0;
}
