// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_LOCALASSEMBLER_HH
#define DUNE_PDELAB_MULTIDOMAIN_LOCALASSEMBLER_HH

#include<map>
#include<tuple>

#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>

#include <dune/pdelab/multidomain/multidomaingridoperatorspaceutilities.hh>
#include <dune/pdelab/multidomain/visitor.hh>
#include <dune/pdelab/multidomain/datawrappers.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


struct SubProblemFilter
{
  template<typename T, std::size_t new_k, std::size_t old_k>
  struct apply
  {
    static const bool value = is_same<typename T::MultiDomainComponentTag,SubProblemTag>::value;
  };
};

struct CouplingFilter
{
  template<typename T, std::size_t new_k, std::size_t old_k>
  struct apply
  {
    static const bool value = is_same<typename T::MultiDomainComponentTag,CouplingTag>::value;
  };
};

namespace functors {

  template<typename Data>
  struct set_time
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Participant>
    void operator()(const Participant& participant)
    {
      participant.setTime(data().time());
    }

  };


  template<typename Data>
  struct set_weight
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Participant>
    void operator()(const Participant& participant)
    {
      participant.setWeight(data().weight());
    }

  };

  template<typename Data>
  struct pre_step
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Participant>
    void operator()(const Participant& participant)
    {
      participant.preStep(data().time(),data().dt(),data().stages());
    }

  };

  template<typename Data>
  struct post_step
  {

    template<typename Participant>
    void operator()(const Participant& participant)
    {
      participant.postStep();
    }

  };

  template<typename Data>
  struct pre_stage
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Participant>
    void operator()(const Participant& participant)
    {
      participant.preStage(data().time(),data().stage());
    }

  };

  template<typename Data>
  struct post_stage
    : public data_accessor<Data>
  {

    template<typename Participant>
    void operator()(const Participant& participant)
    {
      participant.postStage();
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
      participant.postStage();
    }

  };

}


template<...>
class MultiDomainLocalAssembler
{

  bool requireIntersections() const
  {
    return requireUVSkeleton() || requireVSkeleton() ||
      requireUVEnrichedCoupling() || requireVEnrichedCoupling();
  }

  bool requireIntersectionsTwoSided() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,require_intersections_two_sided<> >::value;
  }

  bool requireUVVolume() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_volume<> >::value;
  }

  bool requireVVolume() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_lambda_volume<> >::value;
  }

  bool requireUVSkeleton() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_skeleton<> >::value ||
      any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_boundary<> >::value ||
      any_child<MultiDomainLocalAssembler,Couplings,do_alpha_coupling<> >::value ||
      any_child<MultiDomainLocalAssembler,Couplings,do_alpha_enriched_coupling<> >::value;
  }

  bool requireVSkeleton() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_lambda_skeleton<> >::value ||
      any_child<MultiDomainLocalAssembler,SubProblems,do_lambda_boundary<> >::value ||
      any_child<MultiDomainLocalAssembler,Couplings,do_lambda_coupling<> >::value ||
      any_child<MultiDomainLocalAssembler,Couplings,do_lambda_enriched_coupling<> >::value;
  }

  bool requireUVBoundary() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_boundary<> >::value;
  }

  bool requireVBoundary() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_lambda_boundary<> >::value;
  }

  bool requireUVEnrichedCoupling() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_enriched_coupling<> >::value;
  }

  bool requireVEnrichedCoupling() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_lambda_enriched_coupling<> >::value;
  }

  bool requireUVVolumePostSkeleton() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_volume_post_skeleton<> >::value;
  }

  bool requireVVolumePostSkeleton() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_volume_post_skeleton<> >::value;
  }


  template<typename TReal>
  void setTime(TReal time)
  {
    applyToParticipants(visitor<functors::set_time>::add_data(wrap_time(time)));
  }

  template<typename RF>
  void setWeight(RF weight)
  {
    applyToParticipants(visitor<functors::set_weight>::add_data(wrap_weight(weight)));
  }

  template<typename TReal>
  void preStep (TReal time, TReal dt, std::size_t stages)
  {
    applyToParticipants(visitor<functors::pre_step>::add_data(wrap_time(time),
                                                              wrap_dt(dt),
                                                              wrap_stages(stages)));
  }

  void postStep ()
  {
    applyToParticipants(visitor<functors::post_step>::add_data());
  }

  template<typename TReal>
  void preStage (TReal time, std::size_t stage)
  {
    applyToParticipants(visitor<functors::pre_stage>::add_data(wrap_time(time),
                                                               wrap_stage(stage)));
  }

  void postStage ()
  {
    applyToParticipants(visitor<functors::post_stage>::add_data());
  }


  template<typename TReal>
  TReal suggestTimestep (TReal dt) const
  {
    return applyToParticipants(visitor<functors::suggest_time_step>::add_data(wrap_dt(dt))).dt();
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
    Dune::PDELab::TypeTree::applyToTree(_assemblyParticipants,v);
    return v;
  }

  shared_ptr<AssemblyParticipants> _assemblyParticipants;
  Dune::PDELab::TypeTree::FilteredCompositeNode<AssemblyParticipants,SubProblemFilter> _subProblems;
  Dune::PDELab::TypeTree::FilteredCompositeNode<AssemblyParticipants,CouplingFilter> _couplings;

};


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_LOCALASSEMBLER_HH
