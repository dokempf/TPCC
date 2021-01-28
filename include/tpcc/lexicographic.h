#ifndef TPCC_LEXICOGRAPHIC_H
#define TPCC_LEXICOGRAPHIC_H

#include <tpcc/element.h>

namespace TPCC
{

enum boundaries {both, none, periodic};
/**
 * \brief Lexicographic enumeration of the `k`-dimensional faces in a tensor product chain complex
 * of dimension `n`.
 *
 * \tparam n: The dimension of the tensor product (the order of the tensor)
 * \tparam k: The dimensions of the object considered
 * \tparam Bint: The big integer used for addressing in the whole tensor product
 * \tparam Sint: The small integer used for addressing in each fiber
 * \tparam Tint: The tiny integer with values addressing the fibers
 *
 * The numbering for `k==0` and `k==n` is lexicographic with the first coordinate
 * changing fastest and the last one slowest. For `k` between those values, the set
 * of elements consists of `n` over `k` separate sets, one for each combination of `k`
 * different coordinates.
 *
 * These sets form the slowest moving part of the enumeration, and
 * they are enumerated based such that the elements orthogonal to the first `n-k`
 * coordinates are first, those orthogonal to the last coordinates last.
 *
 * The next level in the enumaration is formed by the fact that these sets
 * can be split into 'sheets', connected sets which are not connected to each other.
 * They are characterized by their coordinates orthogonal to the chosen `k`. Thus, their
 * enumeration follows a lexicographic scheme with first fastest for these coordinates.
 *
 * The fastes level of enumeration is inside each sheet, where again the first coordinates
 * run fastest.
 */
template <int n, int k,
          boundaries bnd = both,
          typename Bint = unsigned int, typename Sint = unsigned short,
          typename Tint = unsigned char>
class Lexicographic
{
  /**
   * \brief The dimension of the fibers in each direction.
   */
  std::array<Sint, n> dimensions;
  /**
   * \brief The number of objects facing the same directions.
   */
  std::array<Bint, binomial(n, k)> block_sizes;

  std::array<Bint, binomial(n, k)> block_sizes_bnd;

public:
  /// The tensor order of the chain complex
  static constexpr Tint order() { return n; }

  /// The dimension of the elements, index in the chain complex
  static constexpr Tint cell_dimension() { return k; }

  /// Index type for addressing in the tensor product
  typedef Bint global_index_t;
  /// Index type for addressing in the fibers
  typedef Sint fiber_index_t;
  /// Index type for addressing dimensions and elements of the chain complex
  typedef Tint dimension_index_t;

  /// The type of elements of this set in the complex
  typedef Element<n, k, Sint, Tint> value_type;

  /**
   * \brief Constructor setting the dimensions of the complex.
   */
  constexpr Lexicographic(const std::array<Sint, n>& d)
    : dimensions(d)
  {
    Combinations<n, k> combinations;
    for (Tint i = 0; i < binomial(n, k); ++i)
    {
      Bint p = 1;
      auto combination = combinations[i];
      for (Tint j = 0; j < k; ++j)
        p *= dimensions[n - 1 - combination.in(j)];
      for (Tint j = 0; j < n - k; ++j)
        switch (bnd)
        {
          case both:
            p *= dimensions[n - 1 - combination.out(j)] + 1; break;
          case none:
            p *= dimensions[n - 1 - combination.out(j)] - 1; break;
          default: // perioridic
            p *= dimensions[n - 1 - combination.out(j)]; break;
        }
      block_sizes_bnd[i] = p;
    }
    for (Tint i = 0; i < binomial(n, k); ++i)
    {
      Bint p = 1;
      auto combination = combinations[i];
      for (Tint j = 0; j < k; ++j)
        p *= dimensions[n - 1 - combination.in(j)];
      for (Tint j = 0; j < n - k; ++j)
        p *= dimensions[n - 1 - combination.out(j)] + 1;
      block_sizes[i] = p;
    }
  }

  /**
   * \brief The number of elements in this set.
   */
  constexpr Bint size() const
  {
    Bint sum = 0;
    for (unsigned int i = 0; i < block_sizes.size(); ++i)
      sum += block_sizes_bnd[i];
    return sum;
  }

  /**
   * \brief The number of elements in one direction
   */
  constexpr Bint block_size(Tint block) const { return block_sizes[block]; }

  /// Dimension of the fiber with given index in the tensor product
  constexpr Sint fiber_dimension(Tint i) const { return dimensions[i]; }

  /**
   * \brief Descriptor for the element at given `index`.
   */
  value_type operator[](Bint index) const;

  Bint index_in_slice(const value_type& e) const
  {
    Bint result = 0, factor = 1;
    for (Tint i = 0; i < k; ++i)
    {
      Sint fdim = dimensions[e.along_direction(i)];
      result += e.along_coordinate(i) * factor;
      factor *= fdim;
    }
    return result;
  }


  /**
   * \brief Find index of a given element.
   */
  Bint index(const value_type& e) const;

  template<boundaries bndT = bnd>
  constexpr Lexicographic<n, k - 1, bndT, Bint, Sint, Tint> boundary() const
  {
    return Lexicographic<n, k - 1, bndT, Bint, Sint, Tint>{ dimensions };
  }
};

template <int n, int k, boundaries bnd, typename Bint, typename Sint, typename Tint>
Element<n, k, Sint, Tint> Lexicographic<n, k, bnd, Bint, Sint, Tint>::operator[](Bint index) const
{
  unsigned int b = 0;
  while (b < block_sizes.size())
  {
    if (index < block_sizes[b])
      break;
    index -= block_sizes[b];
    ++b;
  }

  if (b == block_sizes.size())
    throw(b);

  Combinations<n, k> combinations;
  auto combination = combinations[b];

  std::array<Sint, n> coordinates;
  for (int i = 0; i < k; ++i)
  {
    Sint fdim = dimensions[n - 1 - combination.in(i)];
    coordinates[n - 1 - combination.in(i)] = index % fdim;
    index /= fdim;
  }
  for (int i = 0; i < n - k; ++i)
  {
    Sint fdim = 1 + dimensions[n - 1 - combination.out(i)];
    coordinates[n - 1 - combination.out(i)] = index % fdim;
    index /= fdim;
  }
  return Element<n, k, Sint, Tint>{ combination, coordinates };
}

template <int n, int k, boundaries bnd, typename Bint, typename Sint, typename Tint>
Bint Lexicographic<n, k, bnd, Bint, Sint, Tint>::index(const value_type& e) const
{
  Bint ci = e.direction_index();
  Bint result = 0;
  for (unsigned int i = 0; i < ci; ++i)
    result += block_sizes_bnd[i];

  Bint factor = 1;
  for (Tint i = 0; i < k; ++i)
  {
    Sint fdim = dimensions[e.along_direction(i)];
    result += e.along_coordinate(i) * factor;
    factor *= fdim;
  }
  for (unsigned int i = 0; i < n - k; ++i)
  {
    Tint fdim = 1 + dimensions[e.across_direction(i)];

    if (bnd == periodic)  fdim -= 1;
    else if (bnd == none)  fdim -= 2;

    Bint cross_coord = e.across_coordinate(i);
    assert((cross_coord != 0 && cross_coord != fdim+1) || bnd != none);
    if (cross_coord == fdim && bnd == periodic)
      cross_coord = 0;
    else if (bnd == none)
      --cross_coord;

    result += cross_coord * factor;
    factor *= fdim;
  }
  return result;
}

} // namespace TPCC

#endif
