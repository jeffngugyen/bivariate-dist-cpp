#ifndef DIST_H
#define DIST_H

#include <tuple>
#include <vector>

using BinaryMatrix = std::vector<std::vector<int>>;

// Preserve the original return order:
// (min weight in ker(X) \ row(Z), min weight in ker(Z) \ row(X)).
// With X/Z denoting X/Z stabilizer checks, this is (d_Z, d_X).
// Matrices have 2*l*m columns; redundant rows and empty check sets are allowed.
// Throws std::invalid_argument for invalid/unsupported input, and
// std::domain_error if the code encodes no logical qubits.
std::tuple<int, int> compute_dist(const BinaryMatrix& X, const BinaryMatrix& Z,
                                  int l = 4, int m = 7);

#endif
