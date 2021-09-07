#include <iostream>
#include <complex>
#include <random>
#include <execution>

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

std::vector<std::complex<double>> compute_ground_truth(const std::vector<body_ptr>& bodies)
{
	std::vector<std::complex<double>> us(bodies.size());

	std::for_each(std::execution::par, bodies.begin(), bodies.end(), [&](const body_ptr& body_p)
	{
		std::for_each(bodies.begin(), bodies.end(), [&](const body_ptr& body_q)
		{
			if (body_p->uid == body_q->uid)
			{
				return;
			}

			us[body_p->uid] += kernel_func(body_p->pos, body_q->pos) * body_q->mass;
		});
	});


	return us;
}


int main()
{
	// Initialization of positions/masses
	constexpr size_t num_bodies = 262144; // 262144
	std::vector<body_ptr> bodies;

	for (unsigned i = 0; i < num_bodies; ++i)
	{
		const auto& mass = my_rand() * 1.5;
		//const auto& mass = 1.0;


		const auto& pos = std::complex<double>{my_rand(), my_rand()};

		//const auto width = static_cast<unsigned>(pow(2, 5));
		//const double c = 1.0 / static_cast<double>(width);
		//const auto y = i / width;
		//const auto x = i % width;

		//const auto& pos = std::complex<double>{
		//	x * c + c / 2.0,
		//	y * c + c / 2.0,
		//};


		bodies.push_back(std::make_shared<body<double>>(i, pos, mass));
	}

	// Step 0) build the quadtree 

	std::cout << "Start building the tree..." << std::endl;

	auto qt = quadtree<9>();
	//qt.debug_print();

	std::cout << "	- Inserting nodes..." << std::endl;
	std::for_each(bodies.begin(), bodies.end(), [&](const auto& body)
	{
		qt.allocate_node_for_particle(body);
	});
	std::cout << "Finished building the tree..." << std::endl;

	// Step 1) compute center of mass
	std::cout << "Starting computing COM..." << std::endl;

	qt.compute_com();

	std::cout << "Finished computing COM..." << std::endl;

	// Step 2) Compute multipoles
	std::cout << "Starting computing multipoles..." << std::endl;

	qt.compute_u();

	std::cout << "Finished computing multipoles..." << std::endl;

	// Step 3) Downward pass
	std::cout << "Starting downward pass..." << std::endl;

	qt.downward_pass();

	std::cout << "Finished downward pass..." << std::endl;

	// Step 4) Summing up with local direct N^2 neighbors.
	std::cout << "Starting summation..." << std::endl;

	qt.sum_direct_computation();

	std::cout << "Finished summation..." << std::endl;

	if constexpr (constexpr bool show_rmse = false)
	{
		const auto ground_truth = compute_ground_truth(bodies);

		double sum = 0;
		for (unsigned i = 0; i < num_bodies; ++i)
		{
			sum += pow(ground_truth[i].real() - bodies[i]->u.real(), 2);
		}

		for (int i = 0; i < 50; ++i)
		{
			std::cout << ground_truth[i].real() << " --- " << bodies[i]->u.real() << std::endl;
		}

		constexpr auto n = static_cast<double>(num_bodies);
		const auto rsme = sqrt(sum / n);

		std::cout << "RMSE = " << rsme << std::endl;
	}

	return EXIT_SUCCESS;
}
