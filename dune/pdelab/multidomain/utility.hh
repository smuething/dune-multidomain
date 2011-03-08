#ifndef DUNE_MULTIDOMAIN_UTILITY_HH
#define DUNE_MULTIDOMAIN_UTILITY_HH

#include <type_traits> // for enable_if

namespace Dune {

namespace PDELab {

namespace MultiDomain {

/**
 * TMP for deriving storage types in VariadicCompositeNode
 */
template<typename... OArgs>
struct replace;

template<typename OHead, typename... OTail>
struct replace<OHead,OTail...> {
  template<typename Transform, typename... SArgs>
  struct with {
    typedef typename replace<OTail...>::template with<Transform,SArgs...,typename Transform::template transform<OHead>::type>::type type;
  };
};

/**
 * End of recursion - export the transformed argument list as a tuple
 */
template<>
struct replace<> {
  template<typename Transform, typename... SArgs>
  struct with {
    typedef typename Transform::template container<SArgs...>::type type;
  };
};

/**
 * wrapper to simplify calling the TMP
 */
template<typename Transform, typename... Args>
struct transform {
  typedef typename replace<Args...>::template with<Transform>::type type;
};

/**
 * TMP for deriving storage types in VariadicCompositeNode
 */
template<std::size_t i, typename... OArgs>
struct indexed_replace;

template<std::size_t i, typename OHead, typename... OTail>
struct indexed_replace<i,OHead,OTail...> {
  template<typename Transform, typename... SArgs>
  struct with {
    typedef typename indexed_replace<i+1,OTail...>::template with<Transform,SArgs...,typename Transform::template transform<i,OHead>::type>::type type;
  };
};

/**
 * End of recursion - export the transformed argument list as a tuple
 */
template<std::size_t i>
struct indexed_replace<i> {
  template<typename Transform, typename... SArgs>
  struct with {
    typedef typename Transform::template container<SArgs...>::type type;
  };
};

/**
 * wrapper to simplify calling the TMP
 */
template<typename Transform, typename... Args>
struct indexed_transform {
  typedef typename indexed_replace<0,Args...>::template with<Transform>::type type;
};


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


/*
 * helper functions for the convenience versions of constraints() and interpolate().
 * These functions leave any non-subproblem parameter alone,
 * but replace a subproblem by the correct local function space
 * (either trial or test).
 */

template<typename MultiLFS, typename NoSubProblem>
typename std::enable_if<!is_subproblem<NoSubProblem>::value,const NoSubProblem&>::type
extract_trial_lfs(const MultiLFS& multiLFS, const NoSubProblem& noSubProblem)
{
  return noSubProblem;
}

template<typename MultiLFS, typename SubProblem>
typename std::enable_if<is_subproblem<SubProblem>::value,typename SubProblem::Traits::LocalTrialFunctionSpace>::type
extract_trial_lfs(const MultiLFS& multiLFS, const SubProblem& subProblem)
{
  return typename SubProblem::Traits::LocalTrialFunctionSpace(multiLFS,subProblem,subProblem.trialGridFunctionSpaceConstraints());
}


template<typename MultiLFS, typename NoSubProblem>
typename std::enable_if<!is_subproblem<NoSubProblem>::value,const NoSubProblem&>::type
extract_test_lfs(const MultiLFS& multiLFS, const NoSubProblem& noSubProblem)
{
  return noSubProblem;
}

template<typename MultiLFS, typename SubProblem>
typename std::enable_if<is_subproblem<SubProblem>::value,typename SubProblem::Traits::LocalTestFunctionSpace>::type
extract_test_lfs(const MultiLFS& multiLFS, const SubProblem& subProblem)
{
  return typename SubProblem::Traits::LocalTestFunctionSpace(multiLFS,subProblem,subProblem.testGridFunctionSpaceConstraints());
}

} // namespace MultiDomain

} // namespace PDELab

} // namespace Dune

#endif // DUNE_MULTIDOMAIN_UTILITY_HH
