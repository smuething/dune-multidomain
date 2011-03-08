// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_LOCALASSEMBLER_HH
#define DUNE_PDELAB_MULTIDOMAIN_LOCALASSEMBLER_HH

#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/common/typetree/filteredcompositenode.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>

#include <dune/pdelab/multidomain/multidomaingridoperatorspaceutilities.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/multidomain/visitor.hh>
#include <dune/pdelab/multidomain/datawrappers.hh>
#include <dune/pdelab/multidomain/operatorflagtests.hh>
#include <dune/pdelab/multidomain/residualassemblerengine.hh>
#include <dune/pdelab/multidomain/jacobianassemblerengine.hh>
#include <dune/pdelab/multidomain/patternassemblerengine.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


struct SubProblemFilter
  : public Dune::PDELab::TypeTree::SimpleFilter
{
  template<typename T, std::size_t new_k, std::size_t old_k>
  struct apply
  {
    static const bool value = is_same<typename T::MultiDomainComponentTag,SubProblemTag>::value;
  };
};

struct CouplingFilter
  : public Dune::PDELab::TypeTree::SimpleFilter
{
  template<typename T, std::size_t new_k, std::size_t old_k>
  struct apply
  {
    static const bool value =
      is_same<typename T::MultiDomainComponentTag,CouplingTag>::value ||
      is_same<typename T::MultiDomainComponentTag,EnrichedCouplingTag>::value;
  };
};

namespace functors {

  template<typename Data>
  struct set_time
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Participant>
    void operator()(Participant& participant)
    {
      participant.localOperator().setTime(data().time());
    }

  };


  template<typename Data>
  struct pre_step
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Participant>
    void operator()(Participant& participant)
    {
      participant.localOperator().preStep(data().time(),data().dt(),data().stages());
    }

  };

  template<typename Data>
  struct post_step
  {

    template<typename Participant>
    void operator()(Participant& participant)
    {
      participant.localOperator().postStep();
    }

  };

  template<typename Data>
  struct pre_stage
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Participant>
    void operator()(Participant& participant)
    {
      participant.localOperator().preStage(data().time(),data().stage());
    }

  };

  template<typename Data>
  struct post_stage
    : public data_accessor<Data>
  {

    template<typename Participant>
    void operator()(Participant& participant)
    {
      participant.localOperator().postStage();
    }

  };

  template<typename Data>
  struct suggest_time_step
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Participant>
    void operator()(const Participant& participant)
    {
      data().dt() = participant.localOperator().suggestTimestep(data().dt());
    }

  };

}


template<typename GridOperator, typename... AssemblyParticipants>
class LocalAssembler
  : public Dune::PDELab::TypeTree::VariadicCompositeNode<AssemblyParticipants...>
  , public Dune::PDELab::LocalAssemblerBase<typename GridOperator::Traits::MatrixBackend,
                                            typename GridOperator::Traits::TrialGridFunctionSpaceConstraints,
                                            typename GridOperator::Traits::TestGridFunctionSpaceConstraints>
{

  template<typename>
  friend class ResidualAssemblerEngine;

  template<typename>
  friend class JacobianAssemblerEngine;

  template<typename>
  friend class PatternAssemblerEngine;


  typedef Dune::PDELab::TypeTree::VariadicCompositeNode<AssemblyParticipants...> NodeT;
  typedef Dune::PDELab::LocalAssemblerBase<
    typename GridOperator::Traits::MatrixBackend,
    typename GridOperator::Traits::TrialGridFunctionSpaceConstraints,
    typename GridOperator::Traits::TestGridFunctionSpaceConstraints
    > BaseT;

  typedef Dune::PDELab::TypeTree::FilteredCompositeNode<LocalAssembler,SubProblemFilter> SubProblems;
  typedef Dune::PDELab::TypeTree::FilteredCompositeNode<LocalAssembler,CouplingFilter> Couplings;

public:

  typedef typename GridOperator::Traits::Domain Domain;
  typedef typename GridOperator::Traits::Range Range;
  typedef typename GridOperator::Traits::Jacobian Jacobian;
  typedef typename GridOperator::Traits::MatrixBackend::Pattern Pattern;

  typedef Dune::PDELab::MultiDomain::JacobianAssemblerEngine<LocalAssembler> JacobianAssemblerEngine;
  typedef Dune::PDELab::MultiDomain::ResidualAssemblerEngine<LocalAssembler> ResidualAssemblerEngine;
  typedef Dune::PDELab::MultiDomain::PatternAssemblerEngine<LocalAssembler> PatternAssemblerEngine;

  template<typename TReal>
  void setTime(TReal time)
  {
    applyToParticipants(visitor<functors::set_time>::add_data(store_time(time)));
  }

  void setWeight(double weight)
  {
    _weight = weight;
  }

  template<typename TReal>
  void preStep (TReal time, TReal dt, std::size_t stages)
  {
    applyToParticipants(visitor<functors::pre_step>::add_data(store_time(time),
                                                              store_dt(dt),
                                                              store_stages(stages)));
  }

  void postStep ()
  {
    applyToParticipants(visitor<functors::post_step>::add_data());
  }

  template<typename TReal>
  void preStage (TReal time, std::size_t stage)
  {
    applyToParticipants(visitor<functors::pre_stage>::add_data(store_time(time),
                                                               store_stage(stage)));
  }

  void postStage ()
  {
    applyToParticipants(visitor<functors::post_stage>::add_data());
  }


  template<typename TReal>
  TReal suggestTimestep (TReal dt) const
  {
    return applyToParticipants(visitor<functors::suggest_time_step>::add_data(store_dt(dt))).dt();
  }

  JacobianAssemblerEngine& jacobianAssemblerEngine(const Domain& x, Jacobian& a)
  {
    _jacobianAssemblerEngine.setSolution(x);
    _jacobianAssemblerEngine.setJacobian(a);
    return _jacobianAssemblerEngine;
  }

  ResidualAssemblerEngine& residualAssemblerEngine(const Domain& x, Range& r)
  {
    _residualAssemblerEngine.setSolution(x);
    _residualAssemblerEngine.setResidual(r);
    return _residualAssemblerEngine;
  }

  PatternAssemblerEngine& patternAssemblerEngine(Pattern& pattern)
  {
    _patternAssemblerEngine.setPattern(pattern);
    return _patternAssemblerEngine;
  }

  template<typename Visitor>
  const Visitor& applyToSubProblems(Visitor&& v)
  {
    Dune::PDELab::TypeTree::applyToTree(_subProblems,v);
    return v;
  }

  template<typename Visitor>
  const Visitor& applyToCouplings(Visitor&& v)
  {
    Dune::PDELab::TypeTree::applyToTree(_couplings,v);
    return v;
  }

  template<typename Visitor>
  const Visitor& applyToParticipants(Visitor&& v)
  {
    Dune::PDELab::TypeTree::applyToTree(*this,v);
    return v;
  }

  template<typename Visitor>
  const Visitor& applyToSubProblems(Visitor&& v) const
  {
    Dune::PDELab::TypeTree::applyToTree(_subProblems,v);
    return v;
  }

  template<typename Visitor>
  const Visitor& applyToCouplings(Visitor&& v) const
  {
    Dune::PDELab::TypeTree::applyToTree(_couplings,v);
    return v;
  }

  template<typename Visitor>
  const Visitor& applyToParticipants(Visitor&& v) const
  {
    Dune::PDELab::TypeTree::applyToTree(*this,v);
    return v;
  }

  double weight() const
  {
    return _weight;
  }

  LocalAssembler(AssemblyParticipants&... participants)
    : NodeT(stackobject_to_shared_ptr(participants)...)
    , _subProblems(*this)
    , _couplings(*this)
    , _weight(1.0)
    , _jacobianAssemblerEngine(*this)
    , _residualAssemblerEngine(*this)
    , _patternAssemblerEngine(*this)
  {}

  LocalAssembler(const typename GridOperator::Traits::TrialGridFunctionSpaceConstraints& cu,
                 const typename GridOperator::Traits::TestGridFunctionSpaceConstraints& cv,
                 AssemblyParticipants&... participants)
    : NodeT(stackobject_to_shared_ptr(participants)...)
    , BaseT(cu,cv)
    , _subProblems(*this)
    , _couplings(*this)
    , _weight(1.0)
    , _jacobianAssemblerEngine(*this)
    , _residualAssemblerEngine(*this)
    , _patternAssemblerEngine(*this)
  {}


private:

  SubProblems _subProblems;
  Couplings _couplings;
  double _weight;

  JacobianAssemblerEngine _jacobianAssemblerEngine;
  ResidualAssemblerEngine _residualAssemblerEngine;
  PatternAssemblerEngine _patternAssemblerEngine;

};


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_LOCALASSEMBLER_HH
