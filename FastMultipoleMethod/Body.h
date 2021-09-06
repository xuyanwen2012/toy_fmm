#pragma once

#include <complex>

template <typename T,
          typename = std::enable_if_t<std::is_floating_point_v<T>, T>>
struct body
{
	body(const int uid, const std::complex<T> pos, const double mass)
		: uid(uid), pos(pos), mass(mass)
	{
	}

	int uid;
	std::complex<T> pos;
	T mass;

	T x() { return pos.real(); }
	T y() { return pos.imag(); }
};

using body_ptr = std::shared_ptr<body<double>>;

/// <summary>
/// The main kernel function used to compute the pairwise force between
///	particles. 
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="i">The first position.</param>
/// <param name="j">The second position.</param>
/// <returns>The force.</returns>
template <typename T,
          typename = std::enable_if_t<std::is_floating_point_v<T>, T>>
std::complex<T> kernel_func(const std::complex<T>& i, const std::complex<T>& j)
{
	return log(abs(i - j));
}
