#ifndef DUNE_MULTIDOMAIN_INSTATIONARYCOUPLING_HH
#define DUNE_MULTIDOMAIN_INSTATIONARYCOUPLING_HH

#include <dune/pdelab/multidomain/coupling.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {

template<
  typename TReal,
  typename LocalSubProblem_,
  typename RemoteSubProblem_,
  typename CouplingOperator_
  >
struct InstationaryCouplingTraits :
    public CouplingTraits<LocalSubProblem_,RemoteSubProblem_,CouplingOperator_>
{

  typedef TReal TimeType;
  typedef unsigned int StepType;
  typedef unsigned int StageType;

};

template<
  typename TReal,
  typename LocalSubProblem,
  typename RemoteSubProblem,
  typename CouplingOperator
  >
class InstationaryCoupling :
    public Coupling<LocalSubProblem,RemoteSubProblem,CouplingOperator>
{

  typedef Coupling<LocalSubProblem,RemoteSubProblem,CouplingOperator> BaseT;

public:

  typedef InstationaryCouplingTraits<TReal,LocalSubProblem,RemoteSubProblem,CouplingOperator> Traits;

  InstationaryCoupling(const LocalSubProblem& localSubProblem,
           const RemoteSubProblem& remoteSubProblem,
           CouplingOperator& couplingOperator)
    : BaseT(localSubProblem,remoteSubProblem,couplingOperator)
  {}

  void preStep(typename Traits::TimeType time, typename Traits::TimeType dt, typename Traits::StepType step) {
    this->_operator.preStep(time,dt,step);
  }

  void postStep() {
    this->_operator.postStep();
  }

  void preStage(TReal time, typename Traits::StageType stage) {
    this->_operator.preStage(time,stage);
  }

  void postStage() {
    this->_operator.postStage();
  }

  void setTime(TReal t) const {
    this->_operator.setTime(t);
  }

  typename Traits::TimeType suggestTimestep(typename Traits::TimeType dt) const {
    return this->couplingOperator().suggestTimestep(dt);
  }

};

template<
  typename TReal,
  typename LocalSubProblem,
  typename RemoteSubProblem,
  typename CouplingOperator
  >
struct is_coupling<InstationaryCoupling<TReal,LocalSubProblem,RemoteSubProblem,CouplingOperator> >
{
  static const bool value = true;
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_INSTATIONARYCOUPLING_HH
