#pragma once

#include <array>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <vector>

#include "Body.h"

namespace tree_helper
{
	constexpr std::size_t pow(const std::size_t base, const unsigned exponent)
	{
		return exponent == 0
			       ? 1
			       : exponent % 2 == 0
			       ? pow(base, exponent / 2) * pow(base, exponent / 2)
			       : base * pow(base, (exponent - 1) / 2) * pow(base, (exponent - 1) / 2);
	}

	constexpr std::size_t num_nodes(const std::size_t levels)
	{
		return levels == 0 ? 1 : pow(4, levels) + num_nodes(levels - 1);
	}

	template <std::size_t Base, const unsigned Exp>
	using pow_ = std::integral_constant<std::size_t, pow(Base, Exp)>;

	template <std::size_t Level>
	using num_nodes_ = std::integral_constant<std::size_t, num_nodes(Level - 1)>;
}

struct tree_node
{
	using index_t = unsigned;

	index_t parent{};
	std::array<index_t, 4> children{};
	std::vector<body_ptr> contents{};
	double node_mass{};

	tree_node() = default;

	tree_node(const index_t parent,
	          const index_t a,
	          const index_t b,
	          const index_t c,
	          const index_t d) : parent(parent), children{a, b, c, d}
	{
	}

	void compute_contents_com()
	{
		static auto sum = [&](const double acc, const body_ptr& ptr)
		{
			return acc + ptr->mass;
		};

		node_mass = std::accumulate(contents.begin(), contents.end(), 0.0, sum);
	}
};


template <std::size_t Level>
class quadtree
{
	static_assert(Level <= 10, "Level should be less than 10.");
	static_assert(Level > 0, "Level should be greater than 0.");

public:
	/// <summary>
	/// 
	/// </summary>
	quadtree() : data_(), levels_(Level)
	{
		//  Setup helper tables
		for (unsigned i = 0; i < Level; ++i)
		{
			num_nodes_at_level_[i] = pow(4, i);
		}

		// Setup root
		std::cout << "	- Operating on Root..." << std::endl;
		data_[0] = new tree_node(0, 1, 2, 3, 4);
		auto last_index = 1;

		// Setup tree nodes
		for (unsigned l = last_index; l < levels_ - 1; ++l)
		{
			std::cout << "	- Operating on level " << l << "..." << std::endl;
			std::cout << "		- " << num_nodes_at_level_[l] << " operations" << std::endl;

			const auto num_nodes = num_nodes_at_level_[l];
			const auto width = static_cast<int>(pow(2, l));
			const auto next_width = static_cast<int>(pow(2, l + 1));

			for (int i = 0; i < num_nodes; ++i)
			{
				auto [y, x] = std::div(i, width);

				data_[last_index + i] = new tree_node(
					// TODO: Fix this.
					0,
					//x + y * width,
					last_index + num_nodes + x * 2 + y * 2 * next_width,
					last_index + num_nodes + x * 2 + y * 2 * next_width + 1,
					last_index + num_nodes + x * 2 + (y * 2 + 1) * next_width,
					last_index + num_nodes + x * 2 + (y * 2 + 1) * next_width + 1
				);
			}

			last_index += num_nodes;
		}

		// the last layer should be leaf nodes
		std::cout << "	- Operating on leaf " << levels_ - 1 << " level..." << std::endl;
		std::cout << "		- " << data_.size() - last_index << " operations" << std::endl;

		for (unsigned i = last_index; i < data_.size(); ++i)
		{
			data_[i] = new tree_node();
		}
	}

	/// <summary>
	/// 
	/// </summary>
	/// <param name="particle"></param>
	void allocate_node_for_particle(const body_ptr& particle)
	{
		const double width = pow(2, levels_ - 1);

		const auto x = static_cast<size_t>(width * particle->x());
		const auto y = static_cast<size_t>(width * particle->y());

		const auto ptr_base = data_.size() - num_nodes_at_level_[levels_ - 1];

		const auto index = ptr_base + static_cast<size_t>(x + y * width);
		data_[index]->contents.push_back(particle);
	}

	/// <summary>
	/// 
	/// </summary>
	/// <param name="real_index"></param>
	void debug_print(const bool real_index = false) const
	{
		auto last_index = 0;


		for (auto l = 0; l < levels_; ++l)
		{
			std::cout << "Level " << l << ':';

			const auto width = static_cast<int>(pow(2, l));

			for (int i = 0; i < num_nodes_at_level_[l]; ++i)
			{
				if (i % width == 0)
				{
					std::cout << std::endl;
				}
				if (real_index)
				{
					std::cout << last_index + i << ' ';
				}
				else
				{
					std::cout << i << ' ';
				}
			}
			last_index += num_nodes_at_level_[l];

			std::cout << std::endl;
			std::cout << std::endl;
		}
	}

	/// <summary>
	/// 
	/// </summary>
	void clear() { data_.fill(nullptr); }

	void compute_com()
	{
		std::cout << "	- Operating on level " << levels_ - 1 << " (leafs)..." << std::endl;
		std::cout << "		- " << num_nodes_at_level_[levels_ - 1] << " operations" << std::endl;

		// we start from the leaf nodes
		auto n = data_.size() - 1 - num_nodes_at_level_[levels_ - 1];

		for (unsigned i = data_.size() - 1; i > n; --i)
		{
			data_[i]->compute_contents_com();
		}

		auto last_index = n;
		for (unsigned l = levels_ - 2; l > 0; --l)
		{
			std::cout << "	- Operating on level " << l << "..." << std::endl;
			std::cout << "		- " << num_nodes_at_level_[l] << " operations" << std::endl;

			n = last_index - num_nodes_at_level_[l];
			for (unsigned i = last_index; i > n; --i)
			{
				double sum = 0.0;
				for (unsigned child : data_[i]->children)
				{
					sum += data_[child]->node_mass;
				}
				data_[i]->node_mass = sum;
			}

			last_index = n;
		}

		std::cout << "	- Operating on Root..." << std::endl;
		std::cout << "		- 1 operation" << std::endl;
		double sum = 0.0;
		for (unsigned child : data_[0]->children)
		{
			sum += data_[child]->node_mass;
		}
		data_[0]->node_mass = sum;
	}

	using data_array_t = std::array<tree_node*, tree_helper::num_nodes_<Level>::value>;

protected:
	data_array_t data_;
	std::size_t levels_;

	std::vector<unsigned> get_neighbors(unsigned level, unsigned local_index)
	{
		
	}

private:
	std::array<std::size_t, Level> num_nodes_at_level_;

	unsigned local_xy_to_index(unsigned level, unsigned x, unsigned y)
	{
		auto width = static_cast<unsigned>(pow(2, level));
		return x + y * width;
	}
};

