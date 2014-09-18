// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_LOCALASSEMBLER_HH
#define DUNE_PDELAB_MULTIDOMAIN_LOCALASSEMBLER_HH

#include <dune/typetree/typetree.hh>
#include <dune/typetree/filteredcompositenode.hh>

#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>

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
  : public TypeTree::SimpleFilter
{
  template<typename T, std::size_t new_k, std::size_t old_k>
  struct apply
    : public is_same<typename T::MultiDomainComponentTag,SubProblemTag>
  {};
};

struct CouplingFilter
  : public TypeTree::SimpleFilter
{
  template<typename T, std::size_t new_k, std::size_t old_k>
  struct apply
    : public std::integral_constant<bool,
                                    is_same<typename T::MultiDomainComponentTag,CouplingTag>::value ||
                                    is_same<typename T::MultiDomainComponentTag,EnrichedCouplingTag>::value
                                    >::type
  {};
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


template<typename GO>
struct LocalAssemblerTraits
  : public Dune::PDELab::LocalAssemblerTraits<GO>
{

  struct Spaces
  {

    typedef typename GO::Traits::TrialGridFunctionSpace GFSU;
    typedef typename GO::Traits::TestGridFunctionSpace GFSV;

    typedef LocalFunctionSpace<GFSU> LFSU;
    typedef LocalFunctionSpace<GFSV> LFSV;

    typedef CouplingLocalFunctionSpace<GFSU> LFSU_C;
    typedef CouplingLocalFunctionSpace<GFSV> LFSV_C;

    typedef LFSIndexCache<LFSU,typename LocalAssemblerTraits::TrialGridFunctionSpaceConstraints> LFSU_Cache;
    typedef LFSIndexCache<LFSV,typename LocalAssemblerTraits::TestGridFunctionSpaceConstraints> LFSV_Cache;

    typedef LFSIndexCache<LFSU_C,typename LocalAssemblerTraits::TrialGridFunctionSpaceConstraints> LFSU_C_Cache;
    typedef LFSIndexCache<LFSV_C,typename LocalAssemblerTraits::TestGridFunctionSpaceConstraints> LFSV_C_Cache;

  };

};


template<typename GO, typename... AssemblyParticipants>
class LocalAssembler
  : public TypeTree::CompositeNode<AssemblyParticipants...>
  , public Dune::PDELab::LocalAssemblerBase<typename GO::Traits::MatrixBackend,
                                            typename GO::Traits::TrialGridFunctionSpaceConstraints,
                                            typename GO::Traits::TestGridFunctionSpaceConstraints>
{

  template<typename>
  friend class ResidualAssemblerEngine;

  template<typename>
  friend class JacobianAssemblerEngine;

  template<typename>
  friend class PatternAssemblerEngine;

  template<typename, typename...>
  friend class LocalAssembler;

public:

  typedef LocalAssemblerTraits<GO> Traits;

private:

  typedef TypeTree::CompositeNode<AssemblyParticipants...> NodeT;
  typedef Dune::PDELab::LocalAssemblerBase<
    typename GO::Traits::MatrixBackend,
    typename GO::Traits::TrialGridFunctionSpaceConstraints,
    typename GO::Traits::TestGridFunctionSpaceConstraints
    > BaseT;

  typedef TypeTree::FilteredCompositeNode<LocalAssembler,SubProblemFilter> SubProblems;
  typedef TypeTree::FilteredCompositeNode<LocalAssembler,CouplingFilter> Couplings;

public:

  typedef typename Traits::Solution Domain;
  typedef typename Traits::Residual Range;
  typedef typename Traits::Jacobian Jacobian;
  typedef typename Traits::MatrixPattern Pattern;

  typedef Dune::PDELab::MultiDomain::JacobianAssemblerEngine<LocalAssembler> LocalJacobianAssemblerEngine;
  typedef Dune::PDELab::MultiDomain::ResidualAssemblerEngine<LocalAssembler> LocalResidualAssemblerEngine;
  typedef Dune::PDELab::MultiDomain::PatternAssemblerEngine<LocalAssembler> LocalPatternAssemblerEngine;

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

  LocalJacobianAssemblerEngine& localJacobianAssemblerEngine(Jacobian& a, const Domain& x)
  {
    _jacobianAssemblerEngine.setSolution(x);
    _jacobianAssemblerEngine.setJacobian(a);
    return _jacobianAssemblerEngine;
  }

  LocalResidualAssemblerEngine& localResidualAssemblerEngine(Range& r, const Domain& x)
  {
    _residualAssemblerEngine.setSolution(x);
    _residualAssemblerEngine.setResidual(r);
    return _residualAssemblerEngine;
  }

  LocalPatternAssemblerEngine& localPatternAssemblerEngine(Pattern& pattern)
  {
    _patternAssemblerEngine.setPattern(pattern);
    return _patternAssemblerEngine;
  }

  template<typename Visitor>
  const Visitor& applyToSubProblems(Visitor&& v)
  {
    TypeTree::applyToTree(_subProblems,v);
    return v;
  }

  template<typename Visitor>
  const Visitor& applyToCouplings(Visitor&& v)
  {
    TypeTree::applyToTree(_couplings,v);
    return v;
  }

  template<typename Visitor>
  const Visitor& applyToParticipants(Visitor&& v)
  {
    TypeTree::applyToTree(*this,v);
    return v;
  }

  template<typename Visitor>
  const Visitor& applyToSubProblems(Visitor&& v) const
  {
    TypeTree::applyToTree(_subProblems,v);
    return v;
  }

  template<typename Visitor>
  const Visitor& applyToCouplings(Visitor&& v) const
  {
    TypeTree::applyToTree(_couplings,v);
    return v;
  }

  template<typename Visitor>
  const Visitor& applyToParticipants(Visitor&& v) const
  {
    TypeTree::applyToTree(*this,v);
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

  LocalAssembler(const typename Traits::TrialGridFunctionSpaceConstraints& cu,
                 const typename Traits::TestGridFunctionSpaceConstraints& cv,
                 AssemblyParticipants&... participants)
    : NodeT(stackobject_to_shared_ptr(participants)...)
    , BaseT(cu,cv)
    , _subProblems(*this)
    , _couplings(*this)
    , _weight(1.0)
    , _readData(true)
    , _writeData(true)
    , _jacobianAssemblerEngine(*this)
    , _residualAssemblerEngine(*this)
    , _patternAssemblerEngine(*this)
  {}

  bool readData() const
  {
    return _readData;
  }

  bool writeData() const
  {
    return _writeData;
  }

  template<typename LA>
  void shareData(LA& la)
  {
    la._writeData = false;
    _readData = false;
    _jacobianAssemblerEngine.shareData(la._jacobianAssemblerEngine);
    _patternAssemblerEngine.shareData(la._patternAssemblerEngine);
    _residualAssemblerEngine.shareData(la._residualAssemblerEngine);
  }

private:

  SubProblems _subProblems;
  Couplings _couplings;
  double _weight;
  bool _readData;
  bool _writeData;

  LocalJacobianAssemblerEngine _jacobianAssemblerEngine;
  LocalResidualAssemblerEngine _residualAssemblerEngine;
  LocalPatternAssemblerEngine _patternAssemblerEngine;

};


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_LOCALASSEMBLER_HH
