#ifndef DUNE_MULTIDOMAIN_SUBPROBLEM_HH
#define DUNE_MULTIDOMAIN_SUBPROBLEM_HH

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/multidomain/subdomainset.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


template<
  typename SubProblem,
  typename GFSU,
  typename GFSV,
  typename LOP,
  typename Condition_,
  std::size_t... Indices
  >
struct SubProblemTraits
{

  typedef GFSU MultiDomainTrialGridFunctionSpace;
  typedef Dune::PDELab::LocalFunctionSpace<GFSU> MultiDomainTrialLocalFunctionSpace;
  typedef SubProblemLocalFunctionSpace<MultiDomainTrialLocalFunctionSpace,SubProblem,Indices...> LocalTrialFunctionSpace;
  typedef LocalTrialFunctionSpace TrialLocalFunctionSpace;
  typedef GFSV MultiDomainTestGridFunctionSpace;
  typedef Dune::PDELab::LocalFunctionSpace<GFSV> MultiDomainTestLocalFunctionSpace;
  typedef SubProblemLocalFunctionSpace<MultiDomainTestLocalFunctionSpace,SubProblem,Indices...> LocalTestFunctionSpace;
  typedef LocalTestFunctionSpace TestLocalFunctionSpace;
  typedef Condition_ Condition;
  typedef LOP LocalOperator;

};


struct SubProblemTag;

template<
  typename GFSU,
  typename GFSV,
  typename LocalOperator,
  typename Condition,
  std::size_t... Indices
  >
class SubProblem
  : public Dune::PDELab::TypeTree::LeafNode
{

  dune_static_assert(sizeof...(Indices) >= 1,"You need to provide at least one index");

public:

  typedef SubProblemTag MultiDomainComponentTag;

  typedef SubProblemTraits<SubProblem,GFSU,GFSV,LocalOperator,Condition,Indices...> Traits;

  SubProblem(LocalOperator& lop, const Condition& condition) :
    _lop(lop),
    _condition(condition)
  {}

  LocalOperator& localOperator() {
    return _lop;
  }

  const LocalOperator& localOperator() const {
    return _lop;
  }

  typename Traits::TrialLocalFunctionSpace&& localTrialFunctionSpace(typename Traits::MultiDomainTrialLocalFunctionSpace& mdlfs) {
    return typename Traits::LocalTrialFunctionSpace(mdlfs,_condition);
  }

  typename Traits::TestLocalFunctionSpace&& localTestFunctionSpace(typename Traits::MultiDomainTestLocalFunctionSpace& mdlfs) {
    return typename Traits::LocalTrialFunctionSpace(mdlfs,_condition);
  }

  template<typename EG>
  bool appliesTo(const EG& eg) const {
    return _condition(eg.subDomains());
  }

  const Condition& condition() const {
    return _condition;
  }

  template<typename TReal>
  void setTime(const TReal& time)
  {
    _lop.setTime(time);
  }

protected:
  LocalOperator& _lop;
  const Condition _condition;

};

template<
  typename GFSU,
  typename GFSV,
  typename LocalOperator,
  typename Condition,
  std::size_t... Indices
  >
struct is_subproblem<SubProblem<GFSU,GFSV,LocalOperator,Condition,Indices...> >
{
  static const bool value = true;
};


template<
  typename GFSU,
  typename GFSV,
  typename LocalOperator,
  typename Condition,
  typename... GridFunctionSpaces
  >
class TypeBasedSubProblem
  : public SubProblem<GFSU,
                      GFSV,
                      LocalOperator,
                      Condition,
                      get_subproblem_index<GFSU,GridFunctionSpaces>::value...>
{

  typedef SubProblem<GFSU,
                     GFSV,
                     LocalOperator,
                     Condition,
                     get_subproblem_index<GFSU,GridFunctionSpaces>::value...
                     > BaseT;

public:

  typedef typename BaseT::Traits Traits;

  TypeBasedSubProblem(LocalOperator& lop, const Condition& condition) :
    BaseT(lop,condition)
  {}

};

template<
  typename GFSU,
  typename GFSV,
  typename LocalOperator,
  typename Condition,
  typename... GridFunctionSpaces
  >
struct is_subproblem<TypeBasedSubProblem<GFSU,GFSV,LocalOperator,Condition,GridFunctionSpaces...> >
{
  static const bool value = true;
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune


#endif // DUNE_MULTIDOMAIN_SUBPROBLEM_HH
