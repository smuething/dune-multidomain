// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_RESIDUALASSEMBLERENGINE_HH
#define DUNE_PDELAB_MULTIDOMAIN_RESIDUALASSEMBLERENGINE_HH

#include <algorithm>

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>
#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>

#include <dune/pdelab/multidomain/policy.hh>
#include <dune/pdelab/multidomain/datawrappers.hh>
#include <dune/pdelab/multidomain/residualassemblerfunctors.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


namespace {

  template<typename Domain_View, typename Range_View,
           typename Domain_C_View, typename Range_C_View>
  struct ResidualAssemblerEngineData
  {

    Domain_View domain_s_view;
    Domain_View domain_n_view;
    Domain_C_View domain_c_view;

    Range_View range_s_view;
    Range_View range_n_view;
    Range_C_View range_c_view;

    typedef LocalVector<typename Domain_View::ElementType, TrialSpaceTag> SolutionVector;
    typedef LocalVector<typename Range_View::ElementType, TestSpaceTag> ResidualVector;

    SolutionVector x_s;
    ResidualVector r_s;
    SolutionVector x_n;
    ResidualVector r_n;
    SolutionVector x_c;
    ResidualVector r_c;

    typename ResidualVector::WeightedAccumulationView r_s_view;
    typename ResidualVector::WeightedAccumulationView r_n_view;
    typename ResidualVector::WeightedAccumulationView r_c_view;

    ResidualAssemblerEngineData()
      : r_s_view(r_s.weightedAccumulationView(1.0))
      , r_n_view(r_n.weightedAccumulationView(1.0))
      , r_c_view(r_c.weightedAccumulationView(1.0))
    {}

  };

} // anonymous namespace


template<typename LA>
class ResidualAssemblerEngine
  : public LocalAssemblerEngineBase
{

  template<typename>
  friend class ResidualAssemblerEngine;

public:

  static const bool needs_constraints_caching = false;

  typedef LA LocalAssembler;
  typedef typename LA::Range Range;
  typedef typename LA::Domain Domain;

  struct Traits
  {

    typedef typename LA::Traits::TrialGridFunctionSpaceConstraints TrialGridFunctionSpaceConstraints;
    typedef typename LA::Traits::TestGridFunctionSpaceConstraints TestGridFunctionSpaceConstraints;

    typedef NoConstraintsCachingPolicy CachePolicy;
    typedef typename LA::Traits::template Spaces<CachePolicy> Spaces;

  };

  typedef typename Traits::Spaces::GFSU GFSU;
  typedef typename Traits::Spaces::GFSV GFSV;

  typedef typename Traits::Spaces::LFSU LFSU;
  typedef typename Traits::Spaces::LFSV LFSV;
  typedef typename Traits::Spaces::LFSU_C LFSU_C;
  typedef typename Traits::Spaces::LFSV_C LFSV_C;

  typedef typename Traits::Spaces::LFSU_Cache LFSU_Cache;
  typedef typename Traits::Spaces::LFSV_Cache LFSV_Cache;
  typedef typename Traits::Spaces::LFSU_C_Cache LFSU_C_Cache;
  typedef typename Traits::Spaces::LFSV_C_Cache LFSV_C_Cache;

  typedef typename Domain::template ConstLocalView<LFSU_Cache> Domain_View;
  typedef typename Range::template LocalView<LFSV_Cache> Range_View;

  typedef typename Domain::template ConstLocalView<LFSU_C_Cache> Domain_C_View;
  typedef typename Range::template LocalView<LFSV_C_Cache> Range_C_View;

  typedef ResidualAssemblerEngineData<
    Domain_View,Range_View,
    Domain_C_View,Range_C_View
    > EngineData;

  bool requireSkeleton() const
  {
    return
      requireUVSkeleton() || requireVSkeleton() ||
      requireUVEnrichedCoupling() || requireVEnrichedCoupling() ||
      requireUVBoundary() || requireVBoundary();
  }

  bool requireSkeletonTwoSided() const
  {
    return
      requireUVBoundary() || requireVBoundary() ||
      requireUVEnrichedCoupling() || requireVEnrichedCoupling() ||
      any_child<typename LocalAssembler::SubProblems,do_skeleton_two_sided<> >::value ||
      any_child<typename LocalAssembler::Couplings,do_alpha_coupling<> >::value;;
  }

  bool requireUVVolume() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_alpha_volume<> >::value;
  }

  bool requireVVolume() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_lambda_volume<> >::value;
  }

  bool requireUVSkeleton() const
  {
    return
      any_child<typename LocalAssembler::SubProblems,do_alpha_skeleton<> >::value ||
      any_child<typename LocalAssembler::SubProblems,do_alpha_boundary<> >::value ||
      any_child<typename LocalAssembler::Couplings,do_alpha_coupling<> >::value;
  }

  bool requireVSkeleton() const
  {
    return
      any_child<typename LocalAssembler::SubProblems,do_lambda_skeleton<> >::value ||
      any_child<typename LocalAssembler::SubProblems,do_lambda_boundary<> >::value ||
      any_child<typename LocalAssembler::Couplings,do_lambda_coupling<> >::value;
  }

  bool requireUVBoundary() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_alpha_boundary<> >::value;
  }

  bool requireVBoundary() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_lambda_boundary<> >::value;
  }

  bool requireUVEnrichedCoupling() const
  {
    return any_child<typename LocalAssembler::Couplings,do_alpha_enriched_coupling<> >::value;
  }

  bool requireVEnrichedCoupling() const
  {
    return any_child<typename LocalAssembler::Couplings,do_lambda_enriched_coupling<> >::value;
  }

  bool requireUVVolumePostSkeleton() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_alpha_volume_post_skeleton<> >::value;
  }

  bool requireVVolumePostSkeleton() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_alpha_volume_post_skeleton<> >::value;
  }

  template<typename EG>
  void onBindLFSUV(const EG& eg, const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache)
  {
    if (localAssembler().readData())
      {
        data().x_s.resize(lfsu_s_cache.size());
      }
  }

  template<typename EG>
  void onBindLFSV(const EG& eg, const LFSV_Cache& lfsv_s_cache)
  {
    // clear local residual
    if (localAssembler().readData())
      {
        data().r_s.resize(lfsv_s_cache.size());
        std::fill(data().r_s.base().begin(),data().r_s.base().end(),0.0);
      }
  }


  template<typename IG>
  void onBindLFSUVOutside(const IG& ig,
                          const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache,
                          const LFSU_Cache& lfsu_n_cache, const LFSV_Cache& lfsv_n_cache)
  {
    if (localAssembler().readData())
      {
        data().x_n.resize(lfsu_n_cache.size());
      }
  }

  template<typename IG>
  void onBindLFSVOutside(const IG& ig, const LFSV_Cache& lfsv_s_cache, const LFSV_Cache& lfsv_n_cache)
  {
    // clear local residual
    if (localAssembler().readData())
      {
        data().r_n.resize(lfsv_n_cache.size());
        std::fill(data().r_n.base().begin(),data().r_n.base().end(),0.0);
      }
  }


  template<typename IG>
  void onBindLFSUVCoupling(const IG& ig,
                           const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache,
                           const LFSU_Cache& lfsu_n_cache, const LFSV_Cache& lfsv_n_cache,
                           const LFSU_C_Cache& lfsu_c_cache, const LFSV_C_Cache& lfsv_c_cache)
  {
    if (localAssembler().readData())
      {
        data().x_c.resize(lfsu_c_cache.size());
      }
  }

  template<typename IG>
  void onBindLFSVCoupling(const IG& ig,
                          const LFSV_Cache& lfsv_s_cache,
                          const LFSV_Cache& lfsv_n_cache,
                          const LFSV_C_Cache& lfsv_c_cache)
  {
    // clear local residual
    if (localAssembler().readData())
      {
        data().r_c.resize(lfsv_c_cache.size());
        std::fill(data().r_c.base().begin(),data().r_c.base().end(),0.0);
      }
  }


  template<typename EG>
  void onUnbindLFSV(const EG& eg, const LFSV_Cache & lfsv_s_cache)
  {
    // accumulate local residual into global residual
    if (localAssembler().writeData())
      {
        if (data().r_s_view.modified())
          {
            data().range_s_view.bind(lfsv_s_cache);
            data().range_s_view.add(data().r_s);
            data().range_s_view.commit();
          }
        data().r_s_view.resetModified();
      }
  }

  template<typename IG>
  void onUnbindLFSVOutside(const IG& ig,
                           const LFSV_Cache & lfsv_s_cache,
                           const LFSV_Cache & lfsv_n_cache)
  {
    // accumulate local residual into global residual
    if (localAssembler().writeData())
      {
        if (data().r_n_view.modified())
          {
            data().range_n_view.bind(lfsv_n_cache);
            data().range_n_view.add(data().r_n);
            data().range_n_view.commit();
          }
        data().r_n_view.resetModified();
      }
  }

  template<typename IG>
  void onUnbindLFSVCoupling(const IG& ig,
                            const LFSV_Cache & lfsv_s_cache,
                            const LFSV_Cache & lfsv_n_cache,
                            const LFSV_C_Cache & lfsv_c_cache)
  {
    // accumulate local residual into global residual
    if (localAssembler().writeData())
      {
        if (data().r_c_view.modified())
          {
            data().range_c_view.bind(lfsv_c_cache);
            data().range_c_view.add(data().r_c);
            data().range_c_view.commit();
          }
        data().r_c_view.resetModified();
      }
  }


  void loadCoefficientsLFSUInside(const LFSU_Cache& lfsu_s_cache)
  {
    // read local data
    if (localAssembler().readData())
      {
        data().domain_s_view.bind(lfsu_s_cache);
        data().domain_s_view.read(data().x_s);
      }
  }

  void loadCoefficientsLFSUOutside(const LFSU_Cache& lfsu_n_cache)
  {
    // read local data
    if (localAssembler().readData())
      {
        data().domain_n_view.bind(lfsu_n_cache);
        data().domain_n_view.read(data().x_n);
      }
  }

  void loadCoefficientsLFSUCoupling(const LFSU_C_Cache& lfsu_c_cache)
  {
    // read local data
    if (localAssembler().readData())
      {
        data().domain_c_view.bind(lfsu_c_cache);
        data().domain_c_view.read(data().x_c);
      }
  }


  template<typename EG>
  void assembleUVVolume(const EG& eg, const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache)
  {
    typedef visitor<functors::alpha_volume,do_alpha_volume<> > Visitor;
    data().r_s_view.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_eg(eg),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_x(data().x_s),wrap_r(data().r_s_view)
      )
    );
  }

  template<typename EG>
  void assembleVVolume(const EG& eg, const LFSV_Cache& lfsv_s_cache)
  {
    typedef visitor<functors::lambda_volume,do_lambda_volume<> > Visitor;
    data().r_s_view.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_eg(eg),
        wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_r(data().r_s_view)
      )
    );
  }


  template<typename IG>
  void assembleUVSkeleton(const IG& ig,
                          const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache,
                          const LFSU_Cache& lfsu_n_cache, const LFSV_Cache& lfsv_n_cache)
  {
    typedef visitor<functors::alpha_skeleton_or_boundary,do_alpha_skeleton_or_boundary<> > SubProblemVisitor;
    data().r_s_view.setWeight(localAssembler().weight());
    data().r_n_view.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(
      SubProblemVisitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_x_s(data().x_s),wrap_r_s(data().r_s_view),
        wrap_lfsu_n_cache(lfsu_n_cache),wrap_lfsv_n_cache(lfsv_n_cache),
        wrap_x_n(data().x_n),wrap_r_n(data().r_n_view)
      )
    );

    typedef visitor<functors::alpha_coupling,do_alpha_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(
      CouplingVisitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_x_s(data().x_s),wrap_r_s(data().r_s_view),
        wrap_lfsu_n_cache(lfsu_n_cache),wrap_lfsv_n_cache(lfsv_n_cache),
        wrap_x_n(data().x_n),wrap_r_n(data().r_n_view)
      )
    );
  }

  template<typename IG>
  void assembleVSkeleton(const IG& ig,
                         const LFSV_Cache& lfsv_s_cache,
                         const LFSV_Cache& lfsv_n_cache)
  {
    typedef visitor<functors::lambda_skeleton_or_boundary,do_lambda_skeleton_or_boundary<> > Visitor;
    data().r_s_view.setWeight(localAssembler().weight());
    data().r_n_view.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_r_s(data().r_s_view),
        wrap_lfsv_n_cache(lfsv_n_cache),
        wrap_r_n(data().r_n_view)
      )
    );

    typedef visitor<functors::lambda_coupling,do_lambda_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(
       CouplingVisitor::add_data(
          wrap_operator_type(DefaultOperator()),
          wrap_ig(ig),
          wrap_lfsv_s_cache(lfsv_s_cache),
          wrap_r_s(data().r_s_view),
          wrap_lfsv_n_cache(lfsv_n_cache),
          wrap_r_n(data().r_n_view)
       )
    );
  }


  template<typename IG>
  void assembleUVBoundary(const IG& ig, const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache)
  {
    typedef visitor<functors::alpha_boundary,do_alpha_boundary<> > Visitor;
    data().r_s_view.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_x(data().x_s),wrap_r(data().r_s_view)
      )
    );
  }

  template<typename IG>
  void assembleVBoundary(const IG& ig, const LFSV_Cache& lfsv_s_cache)
  {
    typedef visitor<functors::lambda_boundary,do_lambda_boundary<> > Visitor;
    data().r_s_view.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_r(data().r_s_view)
      )
    );
  }


  template<typename IG>
  void assembleUVEnrichedCoupling(const IG& ig,
                                  const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache,
                                  const LFSU_Cache& lfsu_n_cache, const LFSV_Cache& lfsv_n_cache,
                                  const LFSU_C_Cache& lfsu_c_cache, const LFSV_C_Cache& lfsv_c_cache)
  {
    typedef visitor<functors::alpha_enriched_coupling,do_alpha_enriched_coupling<> > CouplingVisitor;
    data().r_s_view.setWeight(localAssembler().weight());
    data().r_n_view.setWeight(localAssembler().weight());
    data().r_c_view.setWeight(localAssembler().weight());
    localAssembler().applyToCouplings(
      CouplingVisitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_x_s(data().x_s),wrap_r_s(data().r_s_view),
        wrap_lfsu_n_cache(lfsu_n_cache),wrap_lfsv_n_cache(lfsv_n_cache),
        wrap_x_n(data().x_n),wrap_r_n(data().r_n_view),
        wrap_lfsu_c_cache(lfsu_c_cache),wrap_lfsv_c_cache(lfsv_c_cache),
        wrap_x_c(data().x_c),wrap_r_c(data().r_c_view)
      )
    );
  }

  template<typename IG>
  void assembleVEnrichedCoupling(const IG& ig,
                                 const LFSV_Cache& lfsv_s_cache,
                                 const LFSV_Cache& lfsv_n_cache,
                                 const LFSV_C_Cache& lfsv_c_cache)
  {
    typedef visitor<functors::lambda_enriched_coupling,do_lambda_enriched_coupling<> > CouplingVisitor;
    data().r_s_view.setWeight(localAssembler().weight());
    data().r_n_view.setWeight(localAssembler().weight());
    data().r_c_view.setWeight(localAssembler().weight());
    localAssembler().applyToCouplings(
      CouplingVisitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_r_s(data().r_s_view),
        wrap_lfsv_n_cache(lfsv_n_cache),
        wrap_r_n(data().r_n_view),
        wrap_lfsv_c_cache(lfsv_c_cache),
        wrap_r_c(data().r_c_view)
      )
    );
  }


  template<typename EG>
  void assembleUVVolumePostSkeleton(const EG& eg, const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache)
  {
    typedef visitor<functors::alpha_volume_post_skeleton,do_alpha_volume_post_skeleton<> > Visitor;
    data().r_s_view.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_eg(eg),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_x(data().x_s),wrap_r(data().r_s_view)
      )
    );
  }

  template<typename EG>
  void assembleVVolumePostSkeleton(const EG& eg, const LFSV_Cache& lfsv_s_cache)
  {
    typedef visitor<functors::lambda_volume_post_skeleton,do_lambda_volume_post_skeleton<> > Visitor;
    data().r_s_view.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_eg(eg),
        wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_r(data().r_s_view)
      )
    );
  }


  void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
  {
    if (localAssembler().writeData())
      {
        Dune::PDELab::constrain_residual(localAssembler().testConstraints(),data().range_s_view.container());

        // detach container views
        data().domain_s_view.detach();
        data().domain_n_view.detach();
        data().domain_c_view.detach();
        data().range_s_view.detach();
        data().range_n_view.detach();
        data().range_c_view.detach();
      }
  }

  void setSolution(const Domain& x)
  {
    if (localAssembler().readData())
      {
        data().domain_s_view.attach(x);
        data().domain_n_view.attach(x);
        data().domain_c_view.attach(x);
      }
  }

  void setResidual(Range& r)
  {
    if (localAssembler().readData())
      {
        data().range_s_view.attach(r);
        data().range_n_view.attach(r);
        data().range_c_view.attach(r);
      }
  }

  const LocalAssembler& localAssembler() const
  {
    return _localAssembler;
  }

  const typename Traits::TrialGridFunctionSpaceConstraints& trialConstraints() const
  {
    return localAssembler().trialConstraints();
  }

  const typename Traits::TestGridFunctionSpaceConstraints& testConstraints() const
  {
    return localAssembler().testConstraints();
  }

  ResidualAssemblerEngine(const LocalAssembler& localAssembler)
    : _localAssembler(localAssembler)
    , _engineData(make_shared<EngineData>())
  {}

  template<typename Engine>
  void shareData(Engine& engine)
  {
    _engineData = engine._engineData;
  }

private:

  EngineData& data()
  {
    return *_engineData;
  }

  const LocalAssembler& _localAssembler;

  shared_ptr<EngineData> _engineData;

};



} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_RESIDUALASSEMBLERENGINE_HH
