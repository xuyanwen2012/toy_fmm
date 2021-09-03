#pragma once

#include <array>
#include <iostream>
#include <type_traits>

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

	template <std::size_t Level>
	using num_nodes_ = std::integral_constant<std::size_t, num_nodes(Level)>;
}


struct tree_node
{
};

template <std::size_t Level>
class quadtree
{
public:
	constexpr quadtree() : data_(), level_(Level)
	{
	}

	void clear() { data_.fill(nullptr); }

private:
	std::array<tree_node*, tree_helper::num_nodes_<Level>::value> data_;
	std::size_t level_;  
};
