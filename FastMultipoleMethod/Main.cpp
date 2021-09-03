#include <iostream>
#include <complex>
#include <random>

#include "Quadtree.h"

#ifdef _WIN32
double my_rand()
{
	static thread_local std::mt19937 generator; // NOLINT(cert-msc51-cpp)
	const std::uniform_real_distribution<double> distribution(0.0, 1.0);
	return distribution(generator);
}
#endif

#ifdef __linux__
double my_rand(const double f_min = 0.0, const double f_max = 1.0)
{
	const double f = static_cast<double>(rand()) / RAND_MAX;
	return f_min + f * (f_max - f_min);
}
#endif

int main()
{
	constexpr size_t num_bodies = 1024;

	const auto qt = quadtree<5>();

	qt.debug_print();
	// qt.debug_print(true);

	return EXIT_SUCCESS;
}
