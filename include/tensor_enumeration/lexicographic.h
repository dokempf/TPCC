#ifndef TENSOR_ENUMERATION_LEXICOGRAPHIC_H
#define TENSOR_ENUMERATION_LEXICOGRAPHIC_H

#include <tensor_enumeration/combinations.h>

namespace TensorEnumeration
{
  /**
   * \brief The `k`-dimensional faces in a tensor product grid of dimension `n`.
   *
   * \tparam n: The dimension of the tensor product (the order of the tensor)
   * \tparam k: The dimensions of the object considered
   * \tparam Bint: The big integer used for addressing in the whole tensor product
   * \tparam Sint: The small integer used for addressing in each component
   * \tparam Tint: The tiny integer with values addressing the components
   */
  template <int n, int k, typename Bint = unsigned int, typename Sint = unsigned short, typename Tint = unsigned char>
  class Lexicographic
  {
    /**
     * \brief The dimension of the fibers in each direction.
     */
    std::array<Sint, n> dimensions;
    /**
     * \brief The number of objects facing the same directions.
     */
    std::array<Bint, binomial(n,k)> block_sizes;
  public:
    /// The tensor order of the chain complex
    static constexpr Tint order()
    { return n; }

    /// The dimension of the elements
    static constexpr Tint cell_dimension()
    { return k; }

    typedef Bint global_index_t;
    typedef Sint fiber_index_t;
    typedef Tint dimension_index_t;
    
    /**
     * \brief Descriptor for a facet of dimension `k` in the complex
     */
    struct Element
    {
      Combination<n,k> directions;
      std::array<Sint, k> position_along;
      std::array<Sint,n-k> position_across;

      void print_debug(std::ostream& os) const
      {
	directions.print_debug(os);
	os << " (";
	std::array<Sint, n> data;
	for (unsigned int i=0;i<k;++i)
	  data[directions.in(i)] = position_along[i];
	for (unsigned int i=0;i<n-k;++i)
	  data[directions.out(i)] = position_across[i];
	
	for (unsigned int i=0;i<n-1;++i)
	  os << data[i] << ',';
	os << data[n-1] << ")";
      }
    };

    typedef Element value_type;
    
    /**
     * \brief Constructor setting the dimensions of the complex.
     */
    constexpr Lexicographic(const std::array<Sint, n>& d)
      : dimensions(d)
    {
      Combinations<n,k> combinations;
      for (Tint i=0;i<binomial(n,k);++i)
	{
	  Bint p = 1;
	  auto combination = combinations[i];
	  for (Tint j=0;j<k;++j)
	    p *= dimensions[combination.in(j)];
	  for (Tint j=0;j<n-k;++j)
	    p *= 1 + dimensions[combination.out(j)];
	  block_sizes[i] = p;
	}
    }
    
    /**
     * \brief The number of elements in this set.
     */
    constexpr Bint size () const
    {
      Bint sum = 0;
      for (unsigned int i=0;i<block_sizes.size();++i)
	sum += block_sizes[i];
      return sum;
    }

    /**
     * \brief The number of elements in one direction
     */
    constexpr Bint block_size (Tint block) const
    {
      return block_sizes[block];
    }

    /**
     * \brief Descriptor for the element at given `index`.
     */
    Element operator[] (Bint index) const;

    /**
     * \brief Find index of a given element.
     */
    Bint index (const Element& e) const;
  };    

template <int n, int k, typename Bint, typename Sint, typename Tint>
typename Lexicographic<n,k,Bint,Sint,Tint>::Element
Lexicographic<n,k,Bint,Sint,Tint>::operator[] (Bint index) const
{
  unsigned int b=0;
  while (b<block_sizes.size())
    {
      if (index<block_sizes[b])
	break;
      index -= block_sizes[b];
      ++b;
    }

  if (b == block_sizes.size())
    throw (b);
  
  Combinations<n,k> combinations;
  auto combination = combinations[b];

  std::array<Sint, k> along;
  for (unsigned int i=0;i<k;++i)
    {
      Tint fdim = dimensions[combination.in(i)];
      along[i] = index % fdim;
      index /= fdim;
    }
  std::array<Sint, n-k> across;
  for (unsigned int i=0;i<n-k;++i)
    {
      Tint fdim = 1 + dimensions[combination.out(i)];
      across[i] = index % fdim;
      index /= fdim;
    }
  return Element{combination, along, across};
}

template <int n, int k, typename Bint, typename Sint, typename Tint>
Bint
Lexicographic<n,k,Bint,Sint,Tint>::index (const Element& e) const
{
  Bint ci = Combinations<n,k>::index(e.directions);
  Bint result = 0;
  for (unsigned int i=0;i<ci;++i)
    result += block_sizes[i];  

  Bint factor = 1;
  for (unsigned int i=0;i<k;++i)
    {
      Tint fdim = dimensions[e.directions.in(i)];
      result += e.position_along[i] *factor;
      factor *= fdim;
    } 
  for (unsigned int i=0;i<n-k;++i)
    {
      Tint fdim = 1 + dimensions[e.directions.out(i)];
      result += e.position_across[i] *factor;
      factor *= fdim;
    }
  return result;
}
  
}


#endif
