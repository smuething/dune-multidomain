#ifndef DUNE_MULTIDOMAIN_INSTATIONARYSUBPROBLEM_HH
#define DUNE_MULTIDOMAIN_INSTATIONARYSUBPROBLEM_HH

#include <dune/pdelab/multidomain/subproblem.hh>

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
  typename TOP,
  typename Condition_,
  std::size_t... Indices
  >
struct InstationarySubProblemTraits :
    public SubProblemTraits<SubProblem,GFSU,CONU,GFSV,CONV,LOP,Condition_,Indices...>
{

  typedef TOP TemporalOperator;

};


template<
  typename GFSU,
  typename CONU,
  typename GFSV,
  typename CONV,
  typename LocalOperator,
  typename TemporalOperator,
  typename Condition,
  std::size_t... Indices
  >
class InstationarySubProblem :
    public SubProblem<GFSU,CONU,GFSV,CONV,LocalOperator,Condition,Indices...>
{

  typedef SubProblem<GFSU,CONU,GFSV,CONV,LocalOperator,Condition,Indices...> BaseT;

public:

  typedef InstationarySubProblemTraits<InstationarySubProblem,GFSU,CONU,GFSV,CONV,LocalOperator,TemporalOperator,Condition,Indices...> Traits;

  InstationarySubProblem(const CONU& conu, const CONV& conv, const LocalOperator& lop, const TemporalOperator& top, const Condition& condition) :
    BaseT(conu,conv,lop,condition),
    _top(top)
  {}

  const TemporalOperator& temporalOperator() const {
    return _top;
  }

private:
  const TemporalOperator& _top;

};

template<
  typename GFSU,
  typename CONU,
  typename GFSV,
  typename CONV,
  typename LocalOperator,
  typename TemporalOperator,
  typename Condition,
  std::size_t... Indices
  >
struct tag<InstationarySubProblem<GFSU,CONU,GFSV,CONV,LocalOperator,TemporalOperator,Condition,Indices...> >
{
  static const bool value = true;
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune


#endif // DUNE_MULTIDOMAIN_INSTATIONARYSUBPROBLEM_HH
