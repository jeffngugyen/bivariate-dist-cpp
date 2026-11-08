#include "dist.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
int positive_int(const std::string& text) {
    std::size_t used = 0;
    const int value = std::stoi(text, &used);
    if (used != text.size() || value <= 0) {
        throw std::invalid_argument("l and m must be positive integers");
    }
    return value;
}

BinaryMatrix read_matrix(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::invalid_argument("cannot open matrix file: " + path);
    BinaryMatrix matrix;
    std::string line;
    while (std::getline(input, line)) {
        std::istringstream row_input(line);
        std::vector<int> row;
        std::string token;
        while (row_input >> token) {
            if (token != "0" && token != "1") {
                throw std::invalid_argument("matrix files must contain space-separated 0 or 1: " + path);
            }
            row.push_back(token == "1" ? 1 : 0);
        }
        if (!row.empty()) matrix.push_back(std::move(row));
    }
    if (input.bad()) throw std::runtime_error("error reading matrix file: " + path);
    return matrix;
}
}  // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << "Usage: dist H_X.txt H_Z.txt l m\n"
                     "Each row must have 2*l*m binary entries. Empty check files are allowed.\n";
        return 0;
    }
    if (argc != 5) {
        std::cerr << "Usage: dist H_X.txt H_Z.txt l m\nTry: make example\n";
        return 2;
    }
    try {
        const int l = positive_int(argv[3]);
        const int m = positive_int(argv[4]);
        const auto [dz, dx] = compute_dist(read_matrix(argv[1]), read_matrix(argv[2]), l, m);
        std::cout << "d_Z=" << dz << " d_X=" << dx << " distance=" << std::min(dx, dz) << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
