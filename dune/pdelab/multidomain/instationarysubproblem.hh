#ifndef DUNE_MULTIDOMAIN_INSTATIONARYSUBPROBLEM_HH
#define DUNE_MULTIDOMAIN_INSTATIONARYSUBPROBLEM_HH

#include <dune/pdelab/multidomain/subproblem.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


template<
  typename SubProblem,
  typename TReal,
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
  typedef TReal TimeType;
  typedef unsigned int StepType;
  typedef unsigned int StageType;

};


template<
  typename TReal,
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

  typedef InstationarySubProblemTraits<InstationarySubProblem,TReal,GFSU,CONU,GFSV,CONV,LocalOperator,TemporalOperator,Condition,Indices...> Traits;

  InstationarySubProblem(const CONU& conu, const CONV& conv, LocalOperator& lop, TemporalOperator& top, const Condition& condition) :
    BaseT(conu,conv,lop,condition),
    _top(top)
  {}

  const TemporalOperator& temporalOperator() const {
    return _top;
  }

  TemporalOperator& temporalOperator() {
    return _top;
  }

  void preStep(typename Traits::TimeType time, typename Traits::TimeType dt, typename Traits::StepType step) {
    this->localOperator().preStep(time,dt,step);
    this->temporalOperator().preStep(time,dt,step);
  }

  void postStep() {
    this->localOperator().postStep();
    this->temporalOperator().postStep();
  }

  void preStage(TReal time, typename Traits::StageType stage) {
    this->localOperator().preStage(time,stage);
    this->temporalOperator().preStage(time,stage);
  }

  void postStage() {
    this->localOperator().postStage();
    this->temporalOperator().postStage();
  }

  void setTime(TReal t) const {
    this->_lop.setTime(t);
    this->_top.setTime(t);
  }

  typename Traits::TimeType suggestTimestep(typename Traits::TimeType dt) const {
    return this->localOperator().suggestTimestep(dt);
  }

private:
  TemporalOperator& _top;

};

template<
  typename TReal,
  typename GFSU,
  typename CONU,
  typename GFSV,
  typename CONV,
  typename LocalOperator,
  typename TemporalOperator,
  typename Condition,
  std::size_t... Indices
  >
struct is_subproblem<InstationarySubProblem<TReal,GFSU,CONU,GFSV,CONV,LocalOperator,TemporalOperator,Condition,Indices...> >
{
  static const bool value = true;
};


template<
  typename TReal,
  typename GFSU,
  typename CONU,
  typename GFSV,
  typename CONV,
  typename LocalOperator,
  typename TemporalOperator,
  typename Condition,
  typename... GridFunctionSpaces
  >
class TypeBasedInstationarySubProblem
  : public InstationarySubProblem<TReal,
                                  GFSU,
                                  CONU,
                                  GFSV,
                                  CONV,
                                  LocalOperator,
                                  TemporalOperator,
                                  Condition,
                                  get_subproblem_index<GFSU,GridFunctionSpaces>::value...>
{

  typedef InstationarySubProblem<TReal,
                                 GFSU,
                                 CONU,
                                 GFSV,
                                 CONV,
                                 LocalOperator,
                                 TemporalOperator,
                                 Condition,
                                 get_subproblem_index<GFSU,GridFunctionSpaces>::value...
                                 > BaseT;

public:

  typedef typename BaseT::Traits Traits;

  TypeBasedInstationarySubProblem(const CONU& conu, const CONV& conv, LocalOperator& lop, TemporalOperator& top, const Condition& condition) :
    BaseT(conu,conv,lop,top,condition)
  {}

};

template<
  typename TReal,
  typename GFSU,
  typename CONU,
  typename GFSV,
  typename CONV,
  typename LocalOperator,
  typename TemporalOperator,
  typename Condition,
  typename... GridFunctionSpaces
  >
struct is_subproblem<TypeBasedInstationarySubProblem<TReal,GFSU,CONU,GFSV,CONV,LocalOperator,TemporalOperator,Condition,GridFunctionSpaces...> >
{
  static const bool value = true;
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune


#endif // DUNE_MULTIDOMAIN_INSTATIONARYSUBPROBLEM_HH
