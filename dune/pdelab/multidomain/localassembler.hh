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

}


template<...>
class MultiDomainLocalAssembler
{

  bool requireIntersections() const
  {
    return requireAlphaSkeleton() || requireLambdaSkeleton() ||
      requireAlphaEnrichedCoupling() || requireLambdaEnrichedCoupling();
  }

  bool requireIntersectionsTwoSided() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,require_intersections_two_sided<> >::value;
  }

  bool requireAlphaVolume() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_volume<> >::value;
  }

  bool requireLambdaVolume() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_lambda_volume<> >::value;
  }

  bool requireAlphaSkeleton() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_skeleton<> >::value ||
      any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_boundary<> >::value ||
      any_child<MultiDomainLocalAssembler,Couplings,do_alpha_coupling<> >::value ||
      any_child<MultiDomainLocalAssembler,Couplings,do_alpha_enriched_coupling<> >::value;
  }

  bool requireLambdaSkeleton() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_lambda_skeleton<> >::value ||
      any_child<MultiDomainLocalAssembler,SubProblems,do_lambda_boundary<> >::value ||
      any_child<MultiDomainLocalAssembler,Couplings,do_lambda_coupling<> >::value ||
      any_child<MultiDomainLocalAssembler,Couplings,do_lambda_enriched_coupling<> >::value;
  }

  bool requireAlphaBoundary() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_boundary<> >::value;
  }

  bool requireLambdaBoundary() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_lambda_boundary<> >::value;
  }

  bool requireAlphaEnrichedCoupling() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_enriched_coupling<> >::value;
  }

  bool requireLambdaEnrichedCoupling() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_lambda_enriched_coupling<> >::value;
  }

  bool requireAlphaVolumePostSkeleton() const
  {
    return any_child<MultiDomainLocalAssembler,SubProblems,do_alpha_volume_post_skeleton<> >::value;
  }

  bool requireAlphaVolumePostSkeleton() const
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
