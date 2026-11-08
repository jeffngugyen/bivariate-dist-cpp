#include "dist.h"

#include <algorithm>
#include <bitset>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr int max_columns = 128;
using Bits = std::bitset<max_columns>;
using Basis = std::vector<Bits>;

// Compare binary vectors as unsigned integers, highest bit first.
bool less(const Bits& a, const Bits& b) {
    for (int i = max_columns - 1; i >= 0; --i) {
        if (a[i] != b[i]) return !a[i];
    }
    return false;
}

Basis compress_matrix(const BinaryMatrix& matrix, int columns, const char* name) {
    Basis packed;
    packed.reserve(matrix.size());
    for (const auto& row : matrix) {
        if (row.size() != static_cast<std::size_t>(columns)) {
            throw std::invalid_argument(std::string(name) +
                                        " must have exactly 2*l*m columns in every row");
        }
        Bits bits;
        for (int j = 0; j < columns; ++j) {
            if (row[j] != 0 && row[j] != 1) {
                throw std::invalid_argument(std::string(name) + " entries must be 0 or 1");
            }
            bits[j] = row[j] != 0;
        }
        packed.push_back(bits);
    }
    return packed;
}

Basis nullspace(Basis matrix, int columns) {
    std::vector<int> pivots;
    std::vector<bool> is_pivot(columns, false);
    std::size_t rank = 0;
    for (int col = columns - 1; col >= 0; --col) {
        auto pivot = rank;
        while (pivot < matrix.size() && !matrix[pivot][col]) ++pivot;
        if (pivot == matrix.size()) continue;
        std::swap(matrix[rank], matrix[pivot]);
        for (std::size_t i = 0; i < matrix.size(); ++i) {
            if (i != rank && matrix[i][col]) matrix[i] ^= matrix[rank];
        }
        is_pivot[col] = true;
        pivots.push_back(col);
        ++rank;
    }
    Basis basis;
    for (int col = 0; col < columns; ++col) {
        if (is_pivot[col]) continue;
        Bits vector;
        vector.set(col);
        for (std::size_t j = 0; j < rank; ++j) {
            if (matrix[j][col]) vector.set(pivots[j]);
        }
        basis.push_back(vector);
    }
    return basis;
}

// Distinct leading bits, ordered from highest to lowest. XOR reduction leaves
// zero exactly for vectors in this row space; see README for the invariant.
Bits reduce(Bits vector, const Basis& basis) {
    for (const auto& row : basis) {
        if (less(vector ^ row, vector)) vector ^= row;
    }
    return vector;
}

Basis rowspace(const Basis& matrix) {
    Basis basis;
    for (const auto& row : matrix) {
        Bits remainder = reduce(row, basis);
        if (remainder.none()) continue;
        basis.push_back(remainder);
        std::sort(basis.begin(), basis.end(), [](const Bits& a, const Bits& b) {
            return less(b, a);
        });
    }
    return basis;
}

int distance_outside(const Basis& kernel_basis, const Basis& stabilizer_basis) {
    // compute_dist checks both nullities BEFORE either exponential search.
    const auto size = kernel_basis.size();
    const std::uint64_t end = std::uint64_t{1} << size;
    int best = std::numeric_limits<int>::max();
#ifdef _OPENMP
    #pragma omp parallel for reduction(min:best) schedule(static)
#endif
    for (std::uint64_t mask = 1; mask < end; ++mask) {
        Bits vector;
        for (std::size_t i = 0; i < size; ++i) {
            if ((mask >> i) & 1U) vector ^= kernel_basis[i];
        }
        const int weight = static_cast<int>(vector.count());
        // Measure the ORIGINAL vector; reduction is only a membership test.
        if (weight < best && reduce(vector, stabilizer_basis).any()) best = weight;
    }
    return best;
}
}  // namespace

std::tuple<int, int> compute_dist(const BinaryMatrix& X, const BinaryMatrix& Z,
                                  int l, int m) {
    if (l <= 0 || m <= 0 || 2LL * l * m > max_columns) {
        throw std::invalid_argument("l and m must be positive, with 2*l*m <= 128");
    }
    const int columns = 2 * l * m;
    const auto hx = compress_matrix(X, columns, "X");
    const auto hz = compress_matrix(Z, columns, "Z");
    for (const auto& x : hx) {
        for (const auto& z : hz) {
            if ((x & z).count() % 2 != 0) {
                throw std::invalid_argument("CSS checks must commute: X * Z^T = 0 over GF(2)");
            }
        }
    }
    const auto rx = rowspace(hx);
    const auto rz = rowspace(hz);
    if (rx.size() + rz.size() == static_cast<std::size_t>(columns)) {
        throw std::domain_error("code encodes no logical qubits; distance is undefined");
    }
    const auto nx = nullspace(hx, columns);
    const auto nz = nullspace(hz, columns);
    if (nx.size() >= 64 || nz.size() >= 64) {
        throw std::invalid_argument("exhaustive search requires each nullspace dimension < 64");
    }
    return {distance_outside(nx, rz), distance_outside(nz, rx)};
}
