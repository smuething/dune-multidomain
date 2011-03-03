// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_LOCALASSEMBLER_HH
#define DUNE_PDELAB_MULTIDOMAIN_LOCALASSEMBLER_HH

#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/common/typetree/filteredcompositenode.hh>

#include <dune/pdelab/multidomain/multidomaingridoperatorspaceutilities.hh>
#include <dune/pdelab/multidomain/visitor.hh>
#include <dune/pdelab/multidomain/datawrappers.hh>
#include <dune/pdelab/multidomain/operatorflagtests.hh>


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


template<typename... AssemblyParticipants>
class LocalAssembler
  : public Dune::PDELab::TypeTree::VariadicCompositeNode<AssemblyParticipants...>
{

  typedef Dune::PDELab::TypeTree::VariadicCompositeNode<AssemblyParticipants...> NodeT;
  typedef Dune::PDELab::TypeTree::FilteredCompositeNode<LocalAssembler,SubProblemFilter> SubProblems;
  typedef Dune::PDELab::TypeTree::FilteredCompositeNode<LocalAssembler,CouplingFilter> Couplings;

public:

  bool requireIntersections() const
  {
    return requireUVSkeleton() || requireVSkeleton() ||
      requireUVEnrichedCoupling() || requireVEnrichedCoupling() ||
      requireUVBoundary() || requireVBoundary();
  }

  bool requireIntersectionsTwoSided() const
  {
    return requireUVBoundary() || requireVBoundary() ||
      any_child<SubProblems,do_skeleton_two_sided<> >::value;
  }

  bool requireUVVolume() const
  {
    return any_child<SubProblems,do_alpha_volume<> >::value;
  }

  bool requireVVolume() const
  {
    return any_child<SubProblems,do_lambda_volume<> >::value;
  }

  bool requireUVSkeleton() const
  {
    return any_child<SubProblems,do_alpha_skeleton<> >::value ||
      any_child<SubProblems,do_alpha_boundary<> >::value ||
      any_child<Couplings,do_alpha_coupling<> >::value ||
      any_child<Couplings,do_alpha_enriched_coupling<> >::value;
  }

  bool requireVSkeleton() const
  {
    return any_child<SubProblems,do_lambda_skeleton<> >::value ||
      any_child<SubProblems,do_lambda_boundary<> >::value ||
      any_child<Couplings,do_lambda_coupling<> >::value ||
      any_child<Couplings,do_lambda_enriched_coupling<> >::value;
  }

  bool requireUVBoundary() const
  {
    return any_child<SubProblems,do_alpha_boundary<> >::value;
  }

  bool requireVBoundary() const
  {
    return any_child<SubProblems,do_lambda_boundary<> >::value;
  }

  bool requireUVEnrichedCoupling() const
  {
    return any_child<Couplings,do_alpha_enriched_coupling<> >::value;
  }

  bool requireVEnrichedCoupling() const
  {
    return any_child<Couplings,do_lambda_enriched_coupling<> >::value;
  }

  bool requireUVVolumePostSkeleton() const
  {
    return any_child<SubProblems,do_alpha_volume_post_skeleton<> >::value;
  }

  bool requireVVolumePostSkeleton() const
  {
    return any_child<SubProblems,do_alpha_volume_post_skeleton<> >::value;
  }


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

  template<typename Visitor>
  const Visitor& applyToSubProblems(Visitor&& v)
  {
    Dune::PDELab::TypeTree::applyToTree(_subProblems,v);
    return v;
  }

  template<typename Visitor>
  const Visitor& applyToSubCouplings(Visitor&& v)
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
  const Visitor& applyToSubCouplings(Visitor&& v) const
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

  LocalAssembler(AssemblyParticipants&... participants)
    : NodeT(stackobject_to_shared_ptr(participants)...)
    , _subProblems(*this)
    , _couplings(*this)
    , _weight(1.0)
  {}

private:

  SubProblems _subProblems;
  Couplings _couplings;
  double _weight;

};


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_LOCALASSEMBLER_HH
