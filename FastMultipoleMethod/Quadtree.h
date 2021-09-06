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

	unsigned uid{};
	std::array<index_t, 4> children{};
	std::vector<body_ptr> contents{};
	std::vector<index_t> interaction_list{};

	std::complex<double> u{};
	double node_mass{};
	std::complex<double> node_center{}; // could be center of mass, or center of the box?

	tree_node() = default;

	tree_node(const unsigned uid,
	          const index_t a,
	          const index_t b,
	          const index_t c,
	          const index_t d) : uid(uid), children{a, b, c, d}
	{
	}

	/// <summary>
	/// On leaf nodes we should simply compute the sum of all the masses of
	/// particles in this node. 
	/// </summary>
	void compute_contents_com()
	{
		static auto sum = [&](const double acc, const body_ptr& ptr)
		{
			return acc + ptr->mass;
		};

		node_mass = std::accumulate(contents.begin(), contents.end(), 0.0, sum);
	}

	/// <summary>
	/// Should only be called on the leaf_nodes
	/// </summary>
	void distribute_u()
	{
		std::for_each(contents.begin(), contents.end(), [&](const body_ptr& particle)
		{
			particle->u += u;
		});
	}
};


template <std::size_t Level>
class quadtree
{
	static_assert(Level <= 10, "Level should be less than 10.");
	static_assert(Level > 0, "Level should be greater than 0.");

	using index_t = unsigned;
	using data_array_t = std::array<tree_node*, tree_helper::num_nodes_<Level>::value>;

public:
	struct range
	{
		using iter = typename data_array_t::iterator;
		iter b;
		iter e;

		range(iter b, iter e) : b(b), e(e)
		{
		}

		iter begin() { return b; }
		iter end() { return e; }
	};

	constexpr quadtree() : data_(), max_levels_(Level)
	{
		setup_helpers();
		setup_children_table();
		setup_interaction_list();
	}

	range boxes_at_level(size_t l)
	{
		auto [start, n] = level_index_range_[l];
		auto base = data_.begin() + start;
		return range(base, base + n);
	}

	void debug_print(const bool real_index = false) const
	{
		auto last_index = 0;


		for (size_t l = 0; l < max_levels_; ++l)
		{
			std::cout << "Level " << l << ':';

			const auto width = static_cast<int>(pow(2, l));

			auto n = level_index_range_[l].second;
			for (size_t i = 0; i < n; ++i)
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
			last_index += n;

			std::cout << std::endl;
			std::cout << std::endl;
		}
	}

	/// <summary>
	/// 
	/// </summary>
	/// <param name="particle"></param>
	void allocate_node_for_particle(const body_ptr& particle)
	{
		const double width = pow(2, max_levels_ - 1);

		const auto x = static_cast<size_t>(width * particle->x());
		const auto y = static_cast<size_t>(width * particle->y());

		const auto ptr_base = level_index_range_[max_levels_ - 1].first;

		const auto index = ptr_base + static_cast<size_t>(x + y * width);
		data_[index]->contents.push_back(particle);
	}

	void compute_com()
	{
		std::cout << "	- Operating on level " << max_levels_ - 1 << " (leafs)..." << std::endl;
		std::cout << "		- " << level_index_range_[max_levels_ - 1].second << " operations" << std::endl;

		// we start from the leaf nodes
		for (tree_node* node : boxes_at_level(max_levels_ - 1))
		{
			node->compute_contents_com();
		}

		constexpr size_t offset = 1; // I need this offset to bypass unsigned loop
		for (size_t l = max_levels_ - 2 + offset; l >= 0 + offset; --l)
		{
			std::cout << "	- Operating on level " << l - offset << "..." << std::endl;
			std::cout << "		- " << level_index_range_[l - offset].second << " operations" << std::endl;

			for (tree_node* node : boxes_at_level(l - offset))
			{
				// dum way to compute the sum of 
				double sum = 0.0;
				for (unsigned child : node->children)
				{
					sum += data_[child]->node_mass;
				}
				node->node_mass = sum;
			}
		}

		// setup the box center(or center of mass)
		data_[0]->node_center = {0.5, 0.5};
		for (size_t l = 1; l < max_levels_; ++l)
		{
			const auto width = static_cast<unsigned>(pow(2, l));
			const double c = 1.0 / static_cast<double>(width);

			const auto [start_i, n] = level_index_range_[l];
			for (unsigned i = 0; i < n; ++i)
			{
				const auto y = i / width;
				const auto x = i % width;

				data_[start_i + i]->node_center =
				{
					x * c + c / 2.0,
					y * c + c / 2.0,
				};
			}
		}
	}

	void compute_u()
	{
		for (size_t l = 2; l < max_levels_; ++l)
		{
			for (const tree_node* box : boxes_at_level(l)) // O(N)
			{
				for (unsigned other : box->interaction_list) // O(27)
				{
					data_[other]->u += kernel_func(data_[other]->node_center, box->node_center) * box->node_mass;
				}
			}
		}
	}

	void downward_pass()
	{
		// we start from top.
		for (size_t l = 2; l < max_levels_; ++l)
		{
			for (const tree_node* box : boxes_at_level(l))
			{
				for (unsigned child : box->children)
				{
					data_[child]->u += box->u;
				}
			}
		}

		// then distribute the leaf nodes;
		for (tree_node* leaf_nodes : boxes_at_level(max_levels_ - 1))
		{
			leaf_nodes->distribute_u();
		}
	}

protected:
	void setup_helpers()
	{
		index_t last_index = 0;
		for (std::size_t l = 0; l < max_levels_; ++l)
		{
			const auto n = static_cast<index_t>(pow(4, l));
			level_index_range_[l] = std::make_pair(last_index, n);
			last_index = last_index + n;
		}
	}

	void setup_children_table()
	{
		std::cout << "	- Operating on Level 0 (Root)..." << std::endl;
		std::cout << "		- 1 operation" << std::endl;

		data_[0] = new tree_node(0, 1, 2, 3, 4);

		for (size_t l = 1; l < max_levels_ - 1; ++l)
		{
			auto [start_i, n] = level_index_range_[l];

			std::cout << "	- Operating on level " << l << "..." << std::endl;
			std::cout << "		- " << n << " operations" << std::endl;

			const auto width = static_cast<unsigned>(pow(2, l));
			const auto next_width = static_cast<unsigned>(pow(2, l + 1));

			for (unsigned i = 0; i < n; ++i)
			{
				const auto y = i / width;
				const auto x = i % width;

				data_[start_i + i] = new tree_node(
					start_i + i,
					start_i + n + x * 2 + y * 2 * next_width,
					start_i + n + x * 2 + y * 2 * next_width + 1,
					start_i + n + x * 2 + (y * 2 + 1) * next_width,
					start_i + n + x * 2 + (y * 2 + 1) * next_width + 1
				);
			}
		}

		// the last layer should be leaf nodes
		auto [start_i, n] = level_index_range_[max_levels_ - 1];
		std::cout << "	- Operating on leaf " << max_levels_ - 1 << " level..." << std::endl;
		std::cout << "		- " << n << " operations" << std::endl;

		for (unsigned i = start_i; i < data_.size(); ++i)
		{
			data_[i] = new tree_node();
			data_[i]->uid = i;
		}
	}


	void setup_interaction_list()
	{
		// No interaction lists for level 0 and 1, so we starting from level 2.
		// But we will do the for-loop on the parent level. (-1 level)
		for (size_t l = 1; l < max_levels_ - 1; ++l)
		{
			// For each parent
			const auto [start_i, n] = level_index_range_[l];
			const auto [next_level_start_i, _] = level_index_range_[l + 1];
			for (unsigned i = 0; i < n; ++i)
			{
				auto parent_neighbors_children = get_all_neighbors_children(l, start_i, i);

				for (const unsigned child_i : data_[start_i + i]->children)
				{
					std::vector<unsigned> self_neighbor = get_neighbors(l + 1, next_level_start_i,
					                                                    child_i - next_level_start_i);
					std::sort(self_neighbor.begin(), self_neighbor.end());

					std::vector<unsigned> diff;
					diff.reserve(27);

					std::set_difference(
						parent_neighbors_children.begin(), parent_neighbors_children.end(),
						self_neighbor.begin(), self_neighbor.end(),
						std::inserter(diff, diff.begin())
					);

					data_[child_i]->interaction_list = std::move(diff);
				}


				std::vector<unsigned> diff;
			}
		}
	}


private:
	data_array_t data_;
	std::size_t max_levels_;

	std::array<std::pair<index_t, index_t>, Level> level_index_range_;

public:
	/// <summary>
	/// 
	/// </summary>
	/// <param name="level"></param>
	/// <param name="start_i"></param>
	/// <param name="local_index"></param>
	/// <returns>A list of the local indices of the neighbors. </returns>
	[[nodiscard]] std::vector<index_t> get_neighbors(const size_t level, const unsigned start_i,
	                                                 const unsigned local_index) const
	{
		const auto width = static_cast<int>(pow(2, level));
		const auto y = local_index / width;
		const auto x = local_index % width;

		std::vector<index_t> neighbors;
		neighbors.reserve(8);

		for (int i = -1; i <= 1; ++i)
		{
			for (int j = -1; j <= 1; ++j)
			{
				if (i == 0 && j == 0)
				{
					continue;
				}

				const int new_x = static_cast<int>(x) + i;
				const int new_y = static_cast<int>(y) + j;

				if (new_x < 0 || new_x >= width || new_y < 0 || new_y >= width)
				{
					continue;
				}

				neighbors.push_back(start_i + new_x + new_y * width);
			}
		}

		return neighbors;
	}

	/// <summary>
	///
	/// </summary>
	/// <param name="l"></param>
	/// <param name="start_i"></param>
	/// <param name="local_index"></param>
	/// <returns> A list of the absolute indices of the neighbors' children. </returns>
	[[nodiscard]] std::vector<index_t> get_all_neighbors_children(const size_t l, const unsigned start_i,
	                                                              const unsigned local_index) const
	{
		std::vector<unsigned> potential_interaction_list;
		//std::cout << l << ": " << std::endl;
		for (const unsigned neighbor_i : get_neighbors(l, start_i, local_index))
		{
			//std::cout << neighbor_i << ' ';

			const tree_node* neighbor_ptr = data_[neighbor_i];
			for (const unsigned child : neighbor_ptr->children)
			{
				potential_interaction_list.push_back(child);
			}
		}

		//std::cout << std::endl;

		std::sort(potential_interaction_list.begin(), potential_interaction_list.end());
		return potential_interaction_list;
	}
};
