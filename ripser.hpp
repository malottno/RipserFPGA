 
  
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>
 
template <class Key, class T, class H, class E> using hash_map = std::unordered_map<Key, T, H, E>;
template <class Key> using hash = std::hash<Key>;

typedef float value_t;
typedef int64_t index_t;
typedef uint16_t coefficient_t;

#ifdef INDICATE_PROGRESS
static const std::chrono::milliseconds time_step(40);
#endif

static const std::string clear_line("\r\033[K");

static const size_t num_coefficient_bits = 8;

static const index_t max_simplex_index =
    (index_t(1) << (8 * sizeof(index_t) - 1 - num_coefficient_bits)) - 1;

void check_overflow(index_t i) {
	if
#ifdef USE_COEFFICIENTS
	    (i > max_simplex_index)
#else
	    (i < 0)
#endif
		;
}

typedef std::pair<value_t, index_t> diameter_index_t;
value_t get_diameter(const diameter_index_t& i) { return i.first; }
index_t get_index(const diameter_index_t& i) { return i.second; }

typedef std::pair<index_t, value_t> index_diameter_t;
index_t get_index(const index_diameter_t& i) { return i.first; }
value_t get_diameter(const index_diameter_t& i) { return i.second; }


enum compressed_matrix_layout { LOWER_TRIANGULAR, UPPER_TRIANGULAR };

template <compressed_matrix_layout Layout> struct compressed_distance_matrix {
	std::vector<value_t> distances;
	std::vector<value_t*> rows;

	compressed_distance_matrix(std::vector<value_t>&& _distances)
	    : distances(std::move(_distances)), rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
		//assert(distances.size() == size() * (size() - 1) / 2);
		init_rows();
	}

	template <typename DistanceMatrix>
	compressed_distance_matrix(const DistanceMatrix& mat)
	    : distances(mat.size() * (mat.size() - 1) / 2), rows(mat.size()) {
		init_rows();

		for (size_t i = 1; i < size(); ++i)
			for (size_t j = 0; j < i; ++j) rows[i][j] = mat(i, j);
	}

	value_t operator()(const index_t i, const index_t j) const;
	size_t size() const { return rows.size(); }
	void init_rows();
};

typedef compressed_distance_matrix<LOWER_TRIANGULAR> compressed_lower_distance_matrix;
typedef compressed_distance_matrix<UPPER_TRIANGULAR> compressed_upper_distance_matrix;

template <> void compressed_lower_distance_matrix::init_rows() {
	value_t* pointer = &distances[0];
	for (size_t i = 1; i < size(); ++i) {
		rows[i] = pointer;
		pointer += i;
	}
}

template <> void compressed_upper_distance_matrix::init_rows() {
	value_t* pointer = &distances[0] - 1;
	for (size_t i = 0; i < size() - 1; ++i) {
		rows[i] = pointer;
		pointer += size() - i - 2;
	}
}

template <>
value_t compressed_lower_distance_matrix::operator()(const index_t i, const index_t j) const {
	return i == j ? 0 : i < j ? rows[j][i] : rows[i][j];
}

template <>
value_t compressed_upper_distance_matrix::operator()(const index_t i, const index_t j) const {
	return i == j ? 0 : i > j ? rows[j][i] : rows[i][j];
}

struct sparse_distance_matrix {
	std::vector<std::vector<index_diameter_t>> neighbors;

	index_t num_edges;

	sparse_distance_matrix(std::vector<std::vector<index_diameter_t>>&& _neighbors,
	                       index_t _num_edges)
	    : neighbors(std::move(_neighbors)), num_edges(_num_edges) {}

	template <typename DistanceMatrix>
	sparse_distance_matrix(const DistanceMatrix& mat, const value_t threshold)
	    : neighbors(mat.size()), num_edges(0) {

		for (size_t i = 0; i < size(); ++i)
			for (size_t j = 0; j < size(); ++j)
				if (i != j) {
					auto d = mat(i, j);
					if (d <= threshold) {
						++num_edges;
						neighbors[i].push_back({j, d});
					}
				}
	}

	value_t operator()(const index_t i, const index_t j) const {
		auto neighbor =
		    std::lower_bound(neighbors[i].begin(), neighbors[i].end(), index_diameter_t{j, 0});
		return (neighbor != neighbors[i].end() && get_index(*neighbor) == j)
		           ? get_diameter(*neighbor)
		           : std::numeric_limits<value_t>::infinity();
	}

	size_t size() const { return neighbors.size(); }
};




class binomial_coeff_table {
	std::vector<std::vector<index_t>> B;

public:
	binomial_coeff_table(index_t n, index_t k) : B(k + 1, std::vector<index_t>(n + 1, 0)) {
		for (index_t i = 0; i <= n; ++i) {
			B[0][i] = 1;
			for (index_t j = 1; j < std::min(i, k + 1); ++j)
				B[j][i] = B[j - 1][i - 1] + B[j][i - 1];
			if (i <= k) B[i][i] = 1;
			check_overflow(B[std::min(i >> 1, k)][i]);
		}
	}

	index_t operator()(index_t n, index_t k) const {
		//assert(n < B.size() && k < B[n].size() && n >= k - 1);
		return B[k][n];
	}
};


bool is_prime(const coefficient_t n) {
	if (!(n & 1) || n < 2) return n == 2;
	for (coefficient_t p = 3; p <= n / p; p += 2)
		if (!(n % p)) return false;
	return true;
}

std::vector<coefficient_t> multiplicative_inverse_vector(const coefficient_t m) {
	std::vector<coefficient_t> inverse(m);
	inverse[1] = 1;
	// m = a * (m / a) + m % a
	// Multipying with inverse(a) * inverse(m % a):
	// 0 = inverse(m % a) * (m / a) + inverse(a)  (mod m)
	for (coefficient_t a = 2; a < m; ++a) inverse[a] = m - (inverse[m % a] * (m / a)) % m;
	return inverse;
}



class ripserFPGA {
	
	
	
	
	
};
