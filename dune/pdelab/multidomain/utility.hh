#ifndef DUNE_MULTIDOMAIN_UTILITY_HH
#define DUNE_MULTIDOMAIN_UTILITY_HH

namespace Dune {

namespace PDELab {

namespace MultiDomain {

template<int Value, int... Pack>
struct pack_contains;

template<int Value, int First, int... Pack>
struct pack_contains<Value,First,Pack...>
{
  static const bool value = Value == First || pack_contains<Value,Pack...>::value;
};

template<int Value>
struct pack_contains<Value>
{
  static const bool value = false;
};

template<int N, int... ChildIndices>
struct check_indices;

template<int N, int FirstIndex, int... ChildIndices>
struct check_indices<N,FirstIndex,ChildIndices...>
{
  static const bool value = check_indices<N, ChildIndices...>::value && FirstIndex < N && !pack_contains<FirstIndex,ChildIndices...>::value;
};

template<int N>
struct check_indices<N>
{
  static const bool value = true;
};


/*
 * tags for testing whether a given type is a subproblem / coupling
 */

template<typename T>
struct is_subproblem
{
  static const bool value = false;
};

template<typename T>
struct is_coupling
{
  static const bool value = false;
};

} // namespace MultiDomain

} // namespace PDELab

} // namespace Dune

#endif // DUNE_MULTIDOMAIN_UTILITY_HH
