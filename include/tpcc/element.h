#ifndef TPCC_ELEMENT_H
#define TPCC_ELEMENT_H

#include <tpcc/combinations.h>

namespace TPCC
{
template <int n, int k, typename Bint, typename Sint, typename Tint>
class Slab;
/**
 * \brief Tensor coordinates for a facet of dimension `k` in the complex of dimension `n`.
 *
 * Each element is characterized by its orientation and its
 * position in the complex. The orientation is
 * described by a combination of `k` coordinates which span the
 * element. It has no extension across these coordinates.
 *
 * The coordinates of an Element are stored in the array #positions. Its length is `n`and
 * it is ordered the same way as the coordinates in the whole chain complex. Thus, the first entry
 * in there refers to the 'global' first coordinate, no matter whether this is a coordinate
 * direction in which the Element extends. These global coordinates are accessed by operator[]().
 *
 * The `k` local coordinates of the Element in the directions in which it extends are accessed
 * by along_coordinate(). Their global counterparts are numbered by along_direction().
 *
 * The `n-k` coordinates orthogonal to the element are accessed by across_coordinate(). The
 * corresponding global directions are numbered by across_direction().
 *
 * An imbortant feature of Element is that you can obtain an Element object for each of its boundary
 * elements through facet().
 *
 * \tparam n: The dimension of the tensor chain complex
 * \tparam k: The dimension of this element
 * \tparam Sint: The integer type used for indexing in a single fiber
 * \tparam Tint: The integer type used for indexing fibers, directions, etc.
 */
template <int n, int k, typename Sint = unsigned short, typename Tint = unsigned char>
class Element
{
  /**
   * \brief A Combination enumerating the coordinate directions along which the element is aligned.
   *
   * Since combinations are generated in descending order, coordinates usually are referenced in
   * ascending order, the directions are actually not the value in the combination or its
   * complement. The index `i` returned by the combination is immediately replaced by `n-1-i`.
   */
  Combination<n, k> orientation;
  /// The integer coordinates within the `n`-dimensional complex
  std::array<Sint, n> positions;

public:
  /// Constructor with both data elements
  template <typename T>
  Element(const Combination<n, k>& combi, const std::array<T, n>& pos)
    : orientation(combi)
    , positions(pos)
  {
  }

  static constexpr unsigned int n_val = n;
  static constexpr unsigned int k_val = k;

  /// The number of facets in the boundary of this object
  static constexpr Tint n_facets() { return 2 * k; }

  /// The index of the combination enumerating directions
  constexpr Tint direction_index() const { return Combinations<n, k>::index(orientation); }

  /**
   * \brief The coordinates in the `n`-dimensional chain complex
   *
   * \param index between 0 and `n`, the direction in which the coordinate is returned.
   */
  constexpr Sint operator[](Tint index) const { return positions[index]; }

  /**
   * \brief Mapping from local to global coordinate directions.
   *
   * The Element extends along `k` integer coordinates, which are mapped through #orientation
   * to the `n` integer coordinates of the complex.
   * This function maps the local coordinate direction of the Element to the global direction in the
   * chain complex.
   */
  constexpr Tint along_direction(Tint index) const { return n - 1 - orientation.in(index); }

  /// The coordinates in the `k`-dimensional element
  constexpr Sint along_coordinate(Tint index) const { return positions[along_direction(index)]; }

  /**
   * \brief Mapping of orthogonal coordinate directions to global.
   */
  constexpr Tint across_direction(Tint index) const { return n - 1 - orientation.out(index); }

  /// The position of the element orthogonal to its extension
  constexpr Sint across_coordinate(Tint index) const { return positions[across_direction(index)]; }

  /**
   * \brief Function for printing the data stored in the element.
   *
   * The format used is within enclosing parentheses as follows:
   * <ol>
   * <li> Each `along_direction` in order</li>
   * <li> A colon</li>
   * <li> Each `across_direction` in order </li>
   * <li> A space</li>
   * <li> Coordinates from 0 to `n-1` separated by commas</li>
   * </ol>
   * \note This function is also used in the consistency tests,
   * therefore, the comparison data must be changed if the output
   * format is changed.
   */
  void print_debug(std::ostream& os) const
  {
    os << " (";
    for (Tint i = 0; i < k; ++i)
      os << (unsigned int)along_direction(i);
    os << ':';
    for (Tint i = 0; i < n - k; ++i)
      os << (unsigned int)across_direction(i);
    os << ' ';
    for (unsigned int i = 0; i < n - 1; ++i)
      os << positions[i] << ',';
    os << positions[n - 1] << ")";
  }

  /**
   * \brief Enumeration of the boundary facets of the element.
   *
   * Such a boundary facet shares all directions into which it
   * extends, except one. Conversely, it is orthogonal to the same
   * coordinate directions plus one. It has the same one-dimensional
   * coordinates for all shared directions. For the one direction
   * which was eliminated, the coordinate is either the same as for
   * the element for the lower boundary, or plus one for the upper.
   */
  constexpr Element<n, k - 1, Sint, Tint> facet(Tint index) const
  {
    Tint i2 = index / 2;           // The direction index out of k
    Tint im = index % 2;           // Lower or upper boundary in this direction?
    Tint gi = along_direction(i2); // The global direction out of n belonging to index
    Combination<n, k - 1> combi = orientation.eliminate(i2);
    std::array<Sint, n> new_positions = positions;
    if (im == 1)
      ++new_positions[gi];

    return Element<n, k - 1, Sint, Tint>{ combi, new_positions };
  }

  template <int, int, typename, typename, typename>
  friend class Slab;
};
} // namespace TPCC

#endif
