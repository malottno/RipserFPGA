#pragma once

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

enum file_format {
	LOWER_DISTANCE_MATRIX,
	UPPER_DISTANCE_MATRIX,
	DISTANCE_MATRIX,
	POINT_CLOUD,
	DIPHA,
	SPARSE,
	BINARY
};
 
template <class Key, class T, class H, class E> using hash_map = std::unordered_map<Key, T, H, E>;
template <class Key> using hash = std::hash<Key>;

typedef float value_t;
typedef int64_t index_t;
typedef uint16_t coefficient_t;

#ifdef INDICATE_PROGRESS
static const std::chrono::milliseconds time_step(40);
#endif
 
 
#ifdef USE_COEFFICIENTS

#ifdef _MSC_VER
#define PACK(...) __pragma(pack(push, 1)) __VA_ARGS__ __pragma(pack(pop))
#else
#define PACK(...) __attribute__((__packed__)) __VA_ARGS__
#endif

PACK(struct entry_t {
	index_t index : 8 * sizeof(index_t) - num_coefficient_bits;
	coefficient_t coefficient : num_coefficient_bits;
	entry_t(index_t _index, coefficient_t _coefficient)
	    : index(_index), coefficient(_coefficient) {}
	entry_t(index_t _index) : index(_index), coefficient(0) {}
	entry_t() : index(0), coefficient(0) {}
};)

static_assert(sizeof(entry_t) == sizeof(index_t), "size of entry_t is not the same as index_t");

extern entry_t make_entry(index_t i, coefficient_t c);
extern index_t get_index(const entry_t& e);
extern index_t get_coefficient(const entry_t& e);
extern void set_coefficient(entry_t& e, const coefficient_t c);

std::ostream& operator<<(std::ostream& stream, const entry_t& e) {
	stream << get_index(e) << ":" << get_coefficient(e);
	return stream;
}

#else

typedef index_t entry_t;
extern const index_t get_index(const entry_t& i);
extern index_t get_coefficient(const entry_t& i);
extern entry_t make_entry(index_t _index, coefficient_t _value);
extern void set_coefficient(entry_t& e, const coefficient_t c);

#endif
 
 
 


static const std::string clear_line("\r\033[K");

static const size_t num_coefficient_bits = 8;

static const index_t max_simplex_index =
    (index_t(1) << (8 * sizeof(index_t) - 1 - num_coefficient_bits)) - 1;

extern void check_overflow(index_t i);

typedef std::pair<value_t, index_t> diameter_index_t;
extern value_t get_diameter(const diameter_index_t& i);
extern index_t get_index(const diameter_index_t& i);

typedef std::pair<index_t, value_t> index_diameter_t;
extern index_t get_index(const index_diameter_t& i);
extern value_t get_diameter(const index_diameter_t& i);


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

extern template void compressed_lower_distance_matrix::init_rows();
extern template void compressed_upper_distance_matrix::init_rows(); 

extern template
value_t compressed_lower_distance_matrix::operator()(const index_t i, const index_t j) const;

extern template
value_t compressed_upper_distance_matrix::operator()(const index_t i, const index_t j) const;

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


struct euclidean_distance_matrix {
	std::vector<std::vector<value_t>> points;

	euclidean_distance_matrix(std::vector<std::vector<value_t>>&& _points)
	    : points(std::move(_points)) {
		//for (auto p : points) { assert(p.size() == points.front().size()); }
	}

	value_t operator()(const index_t i, const index_t j) const {
		//assert(i < points.size());
		//assert(j < points.size());
		return std::sqrt(std::inner_product(
		    points[i].begin(), points[i].end(), points[j].begin(), value_t(), std::plus<value_t>(),
		    [](value_t u, value_t v) { return (u - v) * (u - v); }));
	}

	size_t size() const { return points.size(); }
};

extern bool is_prime(const coefficient_t n);
extern euclidean_distance_matrix read_point_cloud(std::istream& input_stream);

extern std::vector<coefficient_t> multiplicative_inverse_vector(const coefficient_t m);


void print_usage_and_exit(int);

