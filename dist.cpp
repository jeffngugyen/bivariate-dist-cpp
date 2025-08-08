#include <iostream>
#include <vector>
#include <cstdint>
#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

// Platform-specific 128-bit unsigned int
typedef unsigned __int128 uint128_t;

// Count number of 1-bits
int countSetBits(uint128_t n) {
    uint64_t lo = static_cast<uint64_t>(n);
    uint64_t hi = static_cast<uint64_t>(n >> 64);
    return __builtin_popcountll(lo) + __builtin_popcountll(hi);
}

// Read a binary matrix (as space-separated 0/1, each row one line) into vector<uint128_t>
vector<uint128_t> read_matrix(const string& filename, int rows, int cols) {
    vector<uint128_t> mat;
    ifstream infile(filename);
    string line;
    while (getline(infile, line)) {
        istringstream iss(line);
        uint128_t v = 0;
        for (int j = 0; j < cols; ++j) {
            int bit; iss >> bit;
            v = (v << 1) | (bit & 1);
        }
        mat.push_back(v);
    }
    return mat;
}

// Nullspace of a binary matrix over GF(2), rows as uint128_t, ncols is the number of bits
vector<uint128_t> nullspace(const vector<uint128_t>& mat, int ncols) {
    vector<uint128_t> a = mat;
    int n = a.size();
    vector<int> pivots;
    vector<bool> is_pivot(ncols, false);
    int rank = 0;
    for (int col = ncols - 1; col >= 0; --col) {
        int pivot = -1;
        for (int i = rank; i < n; ++i) {
            if ((a[i] >> col) & 1) { pivot = i; break; }
        }
        if (pivot == -1) continue;
        swap(a[rank], a[pivot]);
        for (int i = 0; i < n; ++i)
            if (i != rank && ((a[i] >> col) & 1))
                a[i] ^= a[rank];
        is_pivot[col] = true;
        pivots.push_back(col);
        rank++;
    }
    // Build basis for the nullspace (free variables)
    vector<uint128_t> basis;
    for (int i = 0; i < ncols; ++i) {
        if (!is_pivot[i]) {
            uint128_t v = uint128_t(1) << i;
            for (int j = 0; j < rank; ++j) {
                int pc = pivots[j];
                if ((a[j] >> i) & 1)
                    v ^= (uint128_t(1) << pc);
            }
            basis.push_back(v);
        }
    }
    return basis;
}

// All binary linear combinations of basis, as a set (for small basis)
set<uint128_t> span(const vector<uint128_t>& basis) {
    set<uint128_t> ans;
    int k = basis.size();
    for (uint64_t mask = 0; mask < (1ULL << k); ++mask) {
        uint128_t v = 0;
        for (int i = 0; i < k; ++i)
            if ((mask >> i) & 1)
                v ^= basis[i];
        ans.insert(v);
    }
    return ans;
}

// Minimum weight of a codeword in ns_basis not in ignore
int span_distance_ignore(const vector<uint128_t>& ns_basis, const set<uint128_t>& ignore) {
    int k = ns_basis.size();
    int d = INT32_MAX;
    for (uint64_t mask = 1; mask < (1ULL << k); ++mask) { // skip zero vector
        uint128_t v = 0;
        for (int i = 0; i < k; ++i)
            if ((mask >> i) & 1)
                v ^= ns_basis[i];
        if (!ignore.count(v))
            d = min(d, countSetBits(v));
    }
    return d;
}

int main() {
    int l = 4, m = 7;
    int nrows = l * m, ncols = 2 * l * m;

    auto H_X = read_matrix("H_X.txt", nrows, ncols);
    auto H_Z = read_matrix("H_Z.txt", nrows, ncols);

    // Find nullspaces
    auto nsX = nullspace(H_X, ncols);
    auto nsZ = nullspace(H_Z, ncols);

    // Find rowspaces
    auto rsX = span(H_X); // for small H_X
    auto rsZ = span(H_Z);

    // Compute distance (minimal X logical not in rowspace of Z-stabilizers, and vice versa)
    int dx = span_distance_ignore(nsX, rsZ);
    int dz = span_distance_ignore(nsZ, rsX);

    cout << "Code distance (X): " << dx << endl;
    cout << "Code distance (Z): " << dz << endl;
    cout << "Minimum code distance: " << min(dx, dz) << endl;
    return 0;
}

