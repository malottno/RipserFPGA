 
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
#include "ripser.hpp"
 
 
int main(int argc, char** argv) {
	const char* filename = nullptr;

	file_format format = DISTANCE_MATRIX;

	index_t dim_max = 1;
	value_t threshold = std::numeric_limits<value_t>::max();
	float ratio = 1;
	coefficient_t modulus = 2;

	for (index_t i = 1; i < argc; ++i) {
		const std::string arg(argv[i]);
		if (arg == "--help") {
			print_usage_and_exit(0);
		} else if (arg == "--dim") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			dim_max = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--threshold") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			threshold = std::stof(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--ratio") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			ratio = std::stof(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--format") {
			std::string parameter = std::string(argv[++i]);
			if (parameter.rfind("lower", 0) == 0)
				format = LOWER_DISTANCE_MATRIX;
			else if (parameter.rfind("upper", 0) == 0)
				format = UPPER_DISTANCE_MATRIX;
			else if (parameter.rfind("dist", 0) == 0)
				format = DISTANCE_MATRIX;
			else if (parameter.rfind("point", 0) == 0)
				format = POINT_CLOUD;
			else if (parameter == "dipha")
				format = DIPHA;
			else if (parameter == "sparse")
				format = SPARSE;
			else if (parameter == "binary")
				format = BINARY;
			else
				print_usage_and_exit(-1);
#ifdef USE_COEFFICIENTS
		} else if (arg == "--modulus") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			modulus = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size() || !is_prime(modulus)) print_usage_and_exit(-1);
#endif
		} else {
			if (filename) { print_usage_and_exit(-1); }
			filename = argv[i];
		}
	}

	std::ifstream file_stream(filename);
	if (filename && file_stream.fail()) {
		std::cerr << "couldn't open file " << filename << std::endl;
		exit(-1);
	}

	if (format == SPARSE) {
		sparse_distance_matrix dist =
		    read_sparse_distance_matrix(filename ? file_stream : std::cin);
		std::cout << "sparse distance matrix with " << dist.size() << " points and "
		          << dist.num_edges << "/" << (dist.size() * (dist.size() - 1)) / 2 << " entries"
		          << std::endl;

		ripser<sparse_distance_matrix>(std::move(dist), dim_max, threshold, ratio, modulus)
		    .compute_barcodes();
	} else if (format == POINT_CLOUD && threshold < std::numeric_limits<value_t>::max()) {
		sparse_distance_matrix dist(read_point_cloud(filename ? file_stream : std::cin), threshold);
		ripser<sparse_distance_matrix>(std::move(dist), dim_max, threshold, ratio, modulus)
				.compute_barcodes();
	} else {
		compressed_lower_distance_matrix dist =
		    read_file(filename ? file_stream : std::cin, format);

		value_t min = std::numeric_limits<value_t>::infinity(),
		        max = -std::numeric_limits<value_t>::infinity(), max_finite = max;
		int num_edges = 0;

		value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
		if (threshold == std::numeric_limits<value_t>::max()) {
			for (size_t i = 0; i < dist.size(); ++i) {
				value_t r_i = -std::numeric_limits<value_t>::infinity();
				for (size_t j = 0; j < dist.size(); ++j) r_i = std::max(r_i, dist(i, j));
				enclosing_radius = std::min(enclosing_radius, r_i);
			}
		}

		for (auto d : dist.distances) {
			min = std::min(min, d);
			max = std::max(max, d);
			if (d != std::numeric_limits<value_t>::infinity()) max_finite = std::max(max_finite, d);
			if (d <= threshold) ++num_edges;
		}
		std::cout << "value range: [" << min << "," << max_finite << "]" << std::endl;

		if (threshold == std::numeric_limits<value_t>::max()) {
			std::cout << "distance matrix with " << dist.size()
			          << " points, using threshold at enclosing radius " << enclosing_radius
			          << std::endl;
			ripser<compressed_lower_distance_matrix>(std::move(dist), dim_max, enclosing_radius,
			                                         ratio, modulus)
			    .compute_barcodes();
		} else {
			std::cout << "sparse distance matrix with " << dist.size() << " points and "
			          << num_edges << "/" << (dist.size() * dist.size() - 1) / 2 << " entries"
			          << std::endl;

			ripser<sparse_distance_matrix>(sparse_distance_matrix(std::move(dist), threshold),
			                               dim_max, threshold, ratio, modulus)
			    .compute_barcodes();
		}
		exit(0);
	}
}
