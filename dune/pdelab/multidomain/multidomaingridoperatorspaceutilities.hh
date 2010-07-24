// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_MULTIDOMAIN_MULTIDOMAINGRIDOPERATORSPACEUTILITIES_HH
#define DUNE_MULTIDOMAIN_MULTIDOMAINGRIDOPERATORSPACEUTILITIES_HH

#include<map>
#include<tuple>

#include<dune/common/exceptions.hh>
#include <dune/pdelab/multidomain/typemap.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {

template<typename GFSU, typename GFSV, typename B>

struct MultiDomainGridOperatorSpaceTraits
{
  typedef GFSU TrialGridFunctionSpace;

  typedef GFSV TestGridFunctionSpace;

  //! \brief the grid view where grid function is defined upon
  typedef typename GFSU::Traits::GridType GridType;
  typedef typename GridType::LeafGridView GridViewType;

  //! \brief vector backend
  typedef B BackendType;

  //! \brief short cut for size type exported by Backend
  typedef typename B::size_type SizeType;
};

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

template<typename T, std::size_t i, std::size_t j>
struct SubProblemEntry
{
  typedef T type;
  typedef T key;
  static const std::size_t globalPos = i;
  static const std::size_t localPos = j;
};

template<std::size_t i, std::size_t j, typename... ProblemsAndCouplings>
struct extract_problems_helper;

template<std::size_t i, std::size_t j, bool is_subproblem_, typename T, typename... ProblemsAndCouplings>
struct extract_problem_switch;

template<std::size_t i, std::size_t j, typename T, typename... ProblemsAndCouplings>
struct extract_problem_switch<i,j,false,T,ProblemsAndCouplings...>
{
  template<typename... Problems>
  struct result
  {
    typedef typename extract_problems_helper<i+1,j,ProblemsAndCouplings...>::template result<Problems...>::type type;
    typedef typename extract_problems_helper<i+1,j,ProblemsAndCouplings...>::template result<Problems...>::map_type map_type;
  };
};

template<std::size_t i, std::size_t j, typename T, typename... ProblemsAndCouplings>
struct extract_problem_switch<i,j,true,T,ProblemsAndCouplings...>
{
  template<typename... Problems>
  struct result
  {
    typedef typename extract_problems_helper<i+1,j+1,ProblemsAndCouplings...>::template result<Problems...,SubProblemEntry<T,i,j> >::type type;
    typedef typename extract_problems_helper<i+1,j+1,ProblemsAndCouplings...>::template result<Problems...,SubProblemEntry<T,i,j> >::map_type map_type;
  };
};

template<std::size_t i, std::size_t j, typename T, typename... ProblemsAndCouplings>
struct extract_problems_helper<i,j,T,ProblemsAndCouplings...>
{
  template<typename... Problems>
  struct result
  {
    typedef typename extract_problem_switch<i,j,is_subproblem<T>::value,T,ProblemsAndCouplings...>::template result<Problems...>::type type;
    typedef typename extract_problem_switch<i,j,is_subproblem<T>::value,T,ProblemsAndCouplings...>::template result<Problems...>::map_type map_type;
  };
};

template<std::size_t i, std::size_t j>
struct extract_problems_helper<i,j>
{
  template<typename... Problems>
  struct result
  {
    typedef std::tuple<Problems...> type;
    typedef embedded_key_map<Problems...> map_type;
  };
};

template<typename... ProblemsAndCouplings>
struct extract_problems
{
  typedef typename extract_problems_helper<0,0,ProblemsAndCouplings...>::template result<>::type type;
  typedef typename extract_problems_helper<0,0,ProblemsAndCouplings...>::template result<>::map_type map_type;
};

template<typename T1, typename T2>
struct lookup_subproblem
{
  typedef T1 type;
};

template<typename T, std::size_t i, std::size_t j, typename SubProblems>
struct CouplingEntry
{
  typedef T type;
  static const std::size_t globalPos = i;
  static const std::size_t localPos = j;
  typedef typename get_map_entry<typename T::Traits::LocalSubProblem,SubProblems>::type LocalSubProblem;
  typedef typename get_map_entry<typename T::Traits::RemoteSubProblem,SubProblems>::type RemoteSubProblem;

};

template<typename SubProblems, std::size_t i, std::size_t j, typename... ProblemsAndCouplings>
struct extract_couplings_helper;

template<typename SubProblems, std::size_t i, std::size_t j, bool is_coupling_, typename T, typename... ProblemsAndCouplings>
struct extract_coupling_switch;

template<typename SubProblems, std::size_t i, std::size_t j, typename T, typename... ProblemsAndCouplings>
struct extract_coupling_switch<SubProblems,i,j,false,T,ProblemsAndCouplings...>
{
  template<typename... Couplings>
  struct result
  {
    typedef typename extract_couplings_helper<SubProblems,i+1,j,ProblemsAndCouplings...>::template result<Couplings...>::type type;
  };
};

template<typename SubProblems, std::size_t i, std::size_t j, typename T, typename... ProblemsAndCouplings>
struct extract_coupling_switch<SubProblems,i,j,true,T,ProblemsAndCouplings...>
{
  template<typename... Couplings>
  struct result
  {
    typedef typename extract_couplings_helper<SubProblems,i+1,j+1,ProblemsAndCouplings...>::template result<Couplings...,CouplingEntry<T,i,j,SubProblems> >::type type;
  };
};

template<typename SubProblems, std::size_t i, std::size_t j, typename T, typename... ProblemsAndCouplings>
struct extract_couplings_helper<SubProblems,i,j,T,ProblemsAndCouplings...>
{
  template<typename... Couplings>
  struct result
  {
    typedef typename extract_coupling_switch<SubProblems,i,j,is_coupling<T>::value,T,ProblemsAndCouplings...>::template result<Couplings...>::type type;
  };
};

template<typename SubProblems,std::size_t i, std::size_t j>
struct extract_couplings_helper<SubProblems,i,j>
{
  template<typename... Couplings>
  struct result
  {
    typedef std::tuple<Couplings...> type;
  };
};

template<typename... ProblemsAndCouplings>
struct extract_couplings
{
  typedef typename extract_couplings_helper<typename extract_problems<ProblemsAndCouplings...>::map_type,
                                            0,
                                            0,
                                            ProblemsAndCouplings...>::template result<>::type type;
};


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_MULTIDOMAINGRIDOPERATORSPACEUTILITIES_HH
