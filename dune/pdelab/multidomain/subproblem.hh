#ifndef DUNE_MULTIDOMAIN_SUBPROBLEM_HH
#define DUNE_MULTIDOMAIN_SUBPROBLEM_HH

namespace Dune {

namespace PDELab {

namespace MultiDomain {


template<
  typename SubProblem,
  typename GFSU,
  typename CONU,
  typename GFSV,
  typename CONV,
  typename LOP,
  typename Condition_,
  std::size_t... Indices
  >
struct SubProblemTraits
{

  typedef GFSU MultiDomainTrialGridFunctionSpace;
  typedef typename GFSU::LocalFunctionSpace MultiDomainLocalTrialFunctionSpace;
  typedef SubProblemLocalFunctionSpace<typename GFSU::LocalFunctionSpace,SubProblem,CONU,Indices...> LocalTrialFunctionSpace;
  typedef LocalTrialFunctionSpace TrialLocalFunctionSpace;
  typedef CONU TrialGridFunctionSpaceConstraints;
  typedef GFSV MultiDomainTestGridFunctionSpace;
  typedef typename GFSV::LocalFunctionSpace MultiDomainLocalTestFunctionSpace;
  typedef SubProblemLocalFunctionSpace<typename GFSV::LocalFunctionSpace,SubProblem,CONV,Indices...> LocalTestFunctionSpace;
  typedef LocalTestFunctionSpace TestLocalFunctionSpace;
  typedef CONV TestGridFunctionSpaceConstraints;
  typedef Condition_ Condition;
  typedef LOP LocalOperator;

};


template<
  typename GFSU,
  typename CONU,
  typename GFSV,
  typename CONV,
  typename LocalOperator,
  typename Condition,
  std::size_t... Indices
  >
class SubProblem
{

  dune_static_assert(sizeof...(Indices) >= 1,"You need to provide at least one index");

public:

  typedef SubProblemTraits<SubProblem,GFSU,CONU,GFSV,CONV,LocalOperator,Condition,Indices...> Traits;

  SubProblem(const CONU& conu, const CONV& conv, LocalOperator& lop, const Condition& condition) :
    _conu(conu),
    _conv(conv),
    _lop(lop),
    _condition(condition)
  {}

  LocalOperator& localOperator() {
    return _lop;
  }

  const LocalOperator& localOperator() const {
    return _lop;
  }

  typename Traits::LocalTrialFunctionSpace&& localTrialFunctionSpace(typename Traits::MultiDomainLocalTrialFunctionSpace& mdlfs) {
    return typename Traits::LocalTrialFunctionSpace(mdlfs,_condition);
  }

  typename Traits::LocalTestFunctionSpace&& localTestFunctionSpace(typename Traits::MultiDomainLocalTestFunctionSpace& mdlfs) {
    return typename Traits::LocalTrialFunctionSpace(mdlfs,_condition);
  }

  template<typename SDS>
  bool appliesTo(const SDS& sds) const {
    return _condition(sds);
  }

  const Condition& condition() const {
    return _condition;
  }

  typename Traits::TrialGridFunctionSpaceConstraints& trialGridFunctionSpaceConstraints() {
    return _conu;
  }

  typename Traits::TestGridFunctionSpaceConstraints& testGridFunctionSpaceConstraints() {
    return _conv;
  }

  const typename Traits::TrialGridFunctionSpaceConstraints& trialGridFunctionSpaceConstraints() const {
    return _conu;
  }

  const typename Traits::TestGridFunctionSpaceConstraints& testGridFunctionSpaceConstraints() const {
    return _conv;
  }


protected:
  CONU _conu;
  CONV _conv;
  LocalOperator& _lop;
  const Condition _condition;

};

template<
  typename GFSU,
  typename CONU,
  typename GFSV,
  typename CONV,
  typename LocalOperator,
  typename Condition,
  std::size_t... Indices
  >
struct tag<SubProblem<GFSU,CONU,GFSV,CONV,LocalOperator,Condition,Indices...> >
{
  static const bool value = true;
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune


#endif // DUNE_MULTIDOMAIN_SUBPROBLEM_HH
