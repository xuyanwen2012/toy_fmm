#pragma once

#include <array>
#include <iostream>
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
	using index_t = size_t;

	std::array<index_t, 4> children{};
	std::vector<body_ptr> contents{};

	tree_node() = default;

	tree_node(const index_t a, const index_t b, const index_t c, const index_t d) : children{a, b, c, d}
	{
	}
};

template <std::size_t Level>
class quadtree
{
	static_assert(Level <= 10, "Level should be less than 10.");
	static_assert(Level > 0, "Level should be greater than 0.");

public:
	quadtree() : data_(), levels_(Level)
	{
		//  Setup helper tables
		for (std::size_t i = 0; i < Level; ++i)
		{
			num_nodes_at_level_[i] = pow(4, i);
		}

		// Setup root
		data_[0] = new tree_node(1, 2, 3, 4);
		auto last_index = 1;

		// Setup tree nodes
		for (auto l = last_index; l < levels_ - 1; ++l)
		{
			const auto num_nodes = num_nodes_at_level_[l];
			const auto width = static_cast<int>(pow(2, l));
			const auto next_width = static_cast<int>(pow(2, l + 1));

			for (int i = 0; i < num_nodes; ++i)
			{
				auto [y, x] = std::div(i, width);

				data_[last_index + i] = new tree_node(
					last_index + num_nodes + x * 2 + y * 2 * next_width,
					last_index + num_nodes + x * 2 + y * 2 * next_width + 1,
					last_index + num_nodes + x * 2 + (y * 2 + 1) * next_width,
					last_index + num_nodes + x * 2 + (y * 2 + 1) * next_width + 1
				);
			}

			last_index += num_nodes;
		}

		// the last layer should be leaf nodes
		for (auto i = last_index; i < data_.size(); ++i)
		{
			data_[i] = new tree_node();
		}
	}

	void allocate_node_for_particle(const body_ptr& particle)
	{
		const double width = pow(2, levels_ - 1);

		const auto x = static_cast<size_t>(width * particle->x());
		const auto y = static_cast<size_t>(width * particle->y());

		const auto ptr_base = data_.size() - num_nodes_at_level_[levels_ - 1];

		const auto index = ptr_base + static_cast<size_t>(x + y * width);
		// std::cout<<index<<std::endl;
		data_[index]->contents.push_back(particle);
	}

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

	void clear() { data_.fill(nullptr); }

protected:
	std::array<tree_node*, tree_helper::num_nodes_<Level>::value> data_;
	std::size_t levels_;

private:
	std::array<std::size_t, Level> num_nodes_at_level_;
};
