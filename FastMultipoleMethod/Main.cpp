#include <iostream>
#include <complex>
#include <random>

#include "Body.h"
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
	// Initialization of positions/masses
	constexpr size_t num_bodies = 1024;
	std::vector<body_ptr> bodies;

	for (size_t i = 0; i < num_bodies; ++i)
	{
		const auto& pos = std::complex<double>{my_rand(), my_rand()};
		const auto& mass = my_rand() * 1.5;

		bodies.push_back(std::make_shared<body<double>>(i, pos, mass));
	}

	auto qt = quadtree<5>();

	for (const auto& body : bodies)
	{
		qt.allocate_node_for_particle(body);
	}

	//qt.debug_print(true);

	return EXIT_SUCCESS;
}
