#include "dist.h"

using std::set;
using std::vector;
using std::ifstream;
using std::string;
using std::istringstream;
using std::cin, std::cout, std::endl;
using std::min, std::max, std::swap;

namespace py=pybind11;

typedef std::bitset<128> bitset;

/*
// std::strong_ordering operator<=> (std::bitset<N> const& a, std::bitset<N> const& b);
*/
struct Compare {
	static std::strong_ordering cmp(bitset const& a, bitset const& b) {
		for(int i = 128-1;i >= 0; --i)
		{
			std::strong_ordering x = a.test(i) <=> b.test(i);
			if(x != std::strong_ordering::equal) return x;
		}
		return std::strong_ordering::equal;
	}

	bool operator()(bitset const& a, bitset const& b) const {return cmp(a,b)<0;}
};

// Count number of 1-bits
int countSetBits(bitset n) {
	return n.count();
}

// Read a binary matrix (as space-separated 0/1, each row one line) into vector<bitset>
vector<bitset> read_matrix(const string& filename, int rows, int cols)
{
    vector<bitset> mat;
    ifstream infile(filename);
    string line;
    while (getline(infile, line)) {
        istringstream iss(line);
        bitset v;
        for (int j = 0; j < cols; ++j) {
            int bit; iss >> bit;
			if(bit & 1) v.set(j);
        }
        mat.push_back(v);
    }
    return mat;
}

vector<bitset> compress_matrix(vector<vector<int> > const& matrix)
{
    vector<bitset> mat;
	for(auto const& v: matrix) {
		bitset x;
		for(int i = 0;i < v.size(); ++i)
			if(v[i] & 1)
				x.set(i);
		mat.push_back(x);
	}
    return mat;
}

// Nullspace of a binary matrix over GF(2), rows as bitset, ncols is the number of bits
vector<bitset> nullspace(const vector<bitset>& mat, int ncols) {
    vector<bitset> a = mat;
    int n = a.size();
    vector<int> pivots;
    vector<bool> is_pivot(ncols, false);
    int rank = 0;
    for (int col = ncols - 1; col >= 0; --col) {
        int pivot = -1;
        for (int i = rank; i < n; ++i) {
            if (a[i].test(col)) { pivot = i; break; }
        }
        if (pivot == -1) continue;
        swap(a[rank], a[pivot]);
        for (int i = 0; i < n; ++i)
            if (i != rank && a[i].test(col))
                a[i] ^= a[rank];
        is_pivot[col] = true;
        pivots.push_back(col);
        rank++;
    }
    // Build basis for the nullspace (free variables)
    vector<bitset> basis;
    for (int i = 0; i < ncols; ++i) {
        if (!is_pivot[i]) {
            bitset v;
			v.set(i);
            for (int j = 0; j < rank; ++j) {
                int pc = pivots[j];
                if (a[j].test(i))
					v.flip(pc);
            }
            basis.push_back(v);
        }
    }
    return basis;
}

// All binary linear combinations of basis, as a set (for small basis)
set<bitset, Compare> span(const vector<bitset>& basis) {
    set<bitset, Compare> ans;
    int k = basis.size();
	assert(k < 63);
    for (uint64_t mask = 0; mask < (1ULL << k); ++mask) {
        bitset v;
        for (int i = 0; i < k; ++i)
            if ((mask >> i) & 1)
                v ^= basis[i];
        ans.insert(v);
    }
    return ans;
}


long long get_time_now() {
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}

// Minimum weight of a codeword in ns_basis not in ignore
int span_distance_ignore(const vector<bitset>& ns_basis, const set<bitset, Compare>& ignore) {
    int k = ns_basis.size();
    cout << "ns basis size " << k << endl;
    cout << "ignore size " << ignore.size() << endl;
    int d = INT32_MAX;

    long long curtime = get_time_now();
    for (uint64_t mask = 1; mask < (1ULL << k); ++mask) { // skip zero vector
        bitset v;
        for (int i = 0; i < k; ++i)
            if ((mask >> i) & 1)
                v ^= ns_basis[i];
        if (!ignore.count(v))
            d = min(d, countSetBits(v));
    }
    long long endtime = get_time_now();
    printf("Time: %.9fs\n", (double)(endtime - curtime) * 1e-6);
    return d;
}

int span_distance_ignore_2(const vector<bitset>& ns_basis, const vector<bitset>& ignore) {
    int k = ns_basis.size();
    cout << "ns basis size " << k << endl;
    cout << "ignore basis size " << ignore.size() << endl;

    long long curtime = get_time_now();

	int nthreads;
	int *results;
	#pragma omp parallel
	{
		#pragma omp single
		{
			nthreads = omp_get_num_threads();
			cout << "starting job with " << nthreads << " threads" << endl;
			results = new int[nthreads];
			for(int i = 0;i < nthreads; ++i)
				results[i] = INT32_MAX;
		}

		#pragma omp for
		for (uint64_t mask = 1; mask < (1ULL << k); ++mask) { // skip zero vector
			bitset v;
			for (int i = 0; i < k; ++i)
				if ((mask >> i) & 1)
					v ^= ns_basis[i];

			int val = countSetBits(v);

			int idx = omp_get_thread_num();
			if(val < results[idx]) {
				for(bitset const& x: ignore)
					if(Compare::cmp(v^x, v) == std::strong_ordering::less) v^=x;
				if(v.any()) {
					results[idx] = val;
					cout << "Current answer(" << idx << "): " << results[idx] << endl;
				}
			}
		}
	}

	int final_result = *std::min_element(results, results + nthreads);
	delete[] results;

    long long endtime = get_time_now();
    printf("Time: %.9fs\n", (double)(endtime - curtime) * 1e-6);
    return final_result;
}

vector<bitset> rowspace(vector<bitset>& v) {
    vector<bitset> ans;
    for(bitset x: v) {
        for(bitset const& b: ans)
            if(Compare::cmp(x^b, x) == std::strong_ordering::less) x^=b;
        if(x.any()) {
            ans.push_back(x);
			std::sort(ans.begin(), ans.end(), [&](bitset const& a, bitset const& b) {return Compare::cmp(a,b)==std::strong_ordering::greater;});
        }
    }
    return ans;
}

std::tuple<int, int> compute_dist(std::vector<std::vector<int> > const& X, std::vector<std::vector<int> > const& Z, int l, int m) {
	py::gil_scoped_release release;
    int nrows = l * m, ncols = 2 * l * m;

    cout << "starting to read" << endl;

    auto H_X = compress_matrix(X);
    auto H_Z = compress_matrix(Z);

    cout << "starting algorithms" << endl;

    // Find nullspaces
    auto nsX = nullspace(H_X, ncols);
    auto nsZ = nullspace(H_Z, ncols);

    cout << "starting rowspace" << endl;

    // Find rowspaces
    // auto rsX = span(H_X); // for small H_X
    // auto rsZ = span(H_Z);

    auto rsX = rowspace(H_X);
    auto rsZ = rowspace(H_Z);

    cout << "starting span dist ignore" << endl;

    // Compute distance (minimal X logical not in rowspace of Z-stabilizers, and vice versa)
    int dx = span_distance_ignore_2(nsX, rsZ);
    int dz = span_distance_ignore_2(nsZ, rsX);

    cout << "Code distance (X): " << dx << endl;
    cout << "Code distance (Z): " << dz << endl;
    cout << "Minimum code distance: " << min(dx, dz) << endl;

	return {dx, dz};
}

int main() {
    int l = 4, m = 7;
    int nrows = l * m, ncols = 2 * l * m;

    cout << "starting to read" << endl;

    auto H_X = read_matrix("H_X.txt", nrows, ncols);
    auto H_Z = read_matrix("H_Z.txt", nrows, ncols);

    cout << "starting algorithms" << endl;

    // Find nullspaces
    auto nsX = nullspace(H_X, ncols);
    auto nsZ = nullspace(H_Z, ncols);

    cout << "starting rowspace" << endl;

    // Find rowspaces
    // auto rsX = span(H_X); // for small H_X
    // auto rsZ = span(H_Z);

    auto rsX = rowspace(H_X);
    auto rsZ = rowspace(H_Z);

    cout << "starting span dist ignore" << endl;

    // Compute distance (minimal X logical not in rowspace of Z-stabilizers, and vice versa)
    int dx = span_distance_ignore_2(nsX, rsZ);
    int dz = span_distance_ignore_2(nsZ, rsX);

    cout << "Code distance (X): " << dx << endl;
    cout << "Code distance (Z): " << dz << endl;
    cout << "Minimum code distance: " << min(dx, dz) << endl;
    return 0;
}

