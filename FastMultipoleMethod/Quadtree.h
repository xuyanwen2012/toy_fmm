#pragma once

#include <array>
#include <iostream>
#include <type_traits>

template <class T>
constexpr T pow(const T base, const unsigned exponent)
{
	return exponent == 0
		       ? 1
		       : exponent % 2 == 0
		       ? pow(base, exponent / 2) * pow(base, exponent / 2)
		       : base * pow(base, (exponent - 1) / 2) * pow(base, (exponent - 1) / 2);
}

template <typename T, T Base, unsigned Exponent>
using pow_ = std::integral_constant<T, pow(Base, Exponent)>;

struct tree_node
{
};

template <std::size_t Level>
class quadtree
{
public:
	quadtree() = default;
private:
	std::array<double, pow_<int, 4, Level>::value> data_;
};

