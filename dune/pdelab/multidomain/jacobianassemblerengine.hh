// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_JACOBIANASSEMBLERENGINE_HH
#define DUNE_PDELAB_MULTIDOMAIN_JACOBIANASSEMBLERENGINE_HH

#include <algorithm>

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>
#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>

#include <dune/pdelab/multidomain/datawrappers.hh>
#include <dune/pdelab/multidomain/jacobianassemblerfunctors.hh>
#include <dune/pdelab/multidomain/operatorflagtests.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {

namespace {

  template<typename Domain_View, typename Domain_C_View,
           typename Jacobian_SS_View, typename Jacobian_SC_View,
           typename Jacobian_CS_View, typename Jacobian_CC_View>
  struct JacobianAssemblerEngineData
  {

    Domain_View domain_s_view;
    Domain_View domain_n_view;
    Domain_C_View domain_c_view;

    typedef LocalVector<typename Domain_View::ElementType, TrialSpaceTag> SolutionVector;
    typedef LocalMatrix<typename Jacobian_SS_View::ElementType> LocalJacobian;

    SolutionVector x_s;
    SolutionVector x_n;
    SolutionVector x_c;

    Jacobian_SS_View jacobian_ss_view;
    Jacobian_SS_View jacobian_sn_view;
    Jacobian_SC_View jacobian_sc_view;

    Jacobian_SS_View jacobian_nn_view;
    Jacobian_SS_View jacobian_ns_view;
    Jacobian_SC_View jacobian_nc_view;

    Jacobian_CC_View jacobian_cc_view;
    Jacobian_CS_View jacobian_cs_view;
    Jacobian_CS_View jacobian_cn_view;

    LocalJacobian a_ss;
    LocalJacobian a_sn;
    LocalJacobian a_sc;

    LocalJacobian a_nn;
    LocalJacobian a_ns;
    LocalJacobian a_nc;

    LocalJacobian a_cc;
    LocalJacobian a_cs;
    LocalJacobian a_cn;

    typename LocalJacobian::WeightedAccumulationView a_ss_view;
    typename LocalJacobian::WeightedAccumulationView a_sn_view;
    typename LocalJacobian::WeightedAccumulationView a_sc_view;

    typename LocalJacobian::WeightedAccumulationView a_nn_view;
    typename LocalJacobian::WeightedAccumulationView a_ns_view;
    typename LocalJacobian::WeightedAccumulationView a_nc_view;

    typename LocalJacobian::WeightedAccumulationView a_cc_view;
    typename LocalJacobian::WeightedAccumulationView a_cs_view;
    typename LocalJacobian::WeightedAccumulationView a_cn_view;

    JacobianAssemblerEngineData()
      : a_ss_view(a_ss,1.0)
      , a_sn_view(a_sn,1.0)
      , a_sc_view(a_sc,1.0)
      , a_nn_view(a_nn,1.0)
      , a_ns_view(a_ns,1.0)
      , a_nc_view(a_nc,1.0)
      , a_cc_view(a_cc,1.0)
      , a_cs_view(a_cs,1.0)
      , a_cn_view(a_cn,1.0)
    {}

  };

} // anonymous namespace

template<typename LA>
class JacobianAssemblerEngine
  : public LocalAssemblerEngineBase
{

  template<typename>
  friend class JacobianAssemblerEngine;

public:

  static const bool needs_constraints_caching = true;

  typedef LA LocalAssembler;
  typedef typename LA::Domain Domain;
  typedef typename LA::Jacobian Jacobian;

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
  typedef typename Domain::template ConstLocalView<LFSU_C_Cache> Domain_C_View;

  typedef typename Jacobian::template LocalView<LFSV_Cache,LFSU_Cache> Jacobian_SS_View;
  typedef typename Jacobian::template LocalView<LFSV_Cache,LFSU_C_Cache> Jacobian_SC_View;
  typedef typename Jacobian::template LocalView<LFSV_C_Cache,LFSU_Cache> Jacobian_CS_View;
  typedef typename Jacobian::template LocalView<LFSV_C_Cache,LFSU_C_Cache> Jacobian_CC_View;

  typedef JacobianAssemblerEngineData<
    Domain_View,Domain_C_View,
    Jacobian_SS_View,Jacobian_SC_View,
    Jacobian_CS_View,Jacobian_CC_View
    > EngineData;

  bool requireSkeleton() const
  {
    return
      requireUVSkeleton() ||
      requireUVEnrichedCoupling() ||
      requireUVBoundary();
  }

  bool requireSkeletonTwoSided() const
  {
    return
      requireUVBoundary() ||
      requireUVEnrichedCoupling() ||
      any_child<typename LocalAssembler::SubProblems,do_skeleton_two_sided<> >::value ||
      any_child<typename LocalAssembler::Couplings,do_alpha_coupling<> >::value;
  }

  bool requireUVVolume() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_alpha_volume<> >::value;
  }

  bool requireUVSkeleton() const
  {
    return
      any_child<typename LocalAssembler::SubProblems,do_alpha_skeleton<> >::value ||
      any_child<typename LocalAssembler::SubProblems,do_alpha_boundary<> >::value ||
      any_child<typename LocalAssembler::Couplings,do_alpha_coupling<> >::value;
  }

  bool requireUVBoundary() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_alpha_boundary<> >::value;
  }

  bool requireUVEnrichedCoupling() const
  {
    return any_child<typename LocalAssembler::Couplings,do_alpha_enriched_coupling<> >::value;
  }

  bool requireUVVolumePostSkeleton() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_alpha_volume_post_skeleton<> >::value;
  }


  template<typename EG>
  void onBindLFSUV(const EG& eg, const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache)
  {
    if (localAssembler().readData())
      {
        data().x_s.resize(lfsu_s_cache.size());
        // initialize local jacobian matrix
        data().a_ss.assign(lfsv_s_cache.size(),lfsu_s_cache.size(),0.0);
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
        // initialize local jacobian matrix
        data().a_nn.assign(lfsv_n_cache.size(),lfsu_n_cache.size(),0.0);
        data().a_sn.assign(lfsv_s_cache.size(),lfsu_n_cache.size(),0.0);
        data().a_ns.assign(lfsv_n_cache.size(),lfsu_s_cache.size(),0.0);
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
        // initialize local jacobian matrix
        data().a_cc.assign(lfsv_c_cache.size(),lfsu_c_cache.size(),0.0);
        data().a_sc.assign(lfsv_s_cache.size(),lfsu_c_cache.size(),0.0);
        data().a_cs.assign(lfsv_c_cache.size(),lfsu_s_cache.size(),0.0);
        data().a_nc.assign(lfsv_n_cache.size(),lfsu_c_cache.size(),0.0);
        data().a_cn.assign(lfsv_c_cache.size(),lfsu_n_cache.size(),0.0);
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
  void onUnbindLFSUV(const EG& eg, const LFSU_Cache & lfsu_s_cache, const LFSV_Cache & lfsv_s_cache)
  {
    // write back local jacobian contributions
    if (localAssembler().writeData())
      {
        if (data().a_ss_view.modified())
          {
            data().jacobian_ss_view.bind(lfsv_s_cache,lfsu_s_cache);
            localAssembler().etadd(data().a_ss,data().jacobian_ss_view);
            data().jacobian_ss_view.commit();
          }
        data().a_ss_view.resetModified();
      }
  }

  template<typename IG>
  void onUnbindLFSUVOutside(const IG& ig,
                            const LFSU_Cache & lfsu_s_cache, const LFSV_Cache & lfsv_s_cache,
                            const LFSU_Cache & lfsu_n_cache, const LFSV_Cache & lfsv_n_cache)
  {
    // write back local jacobian contributions
    if (localAssembler().writeData())
      {
        if (data().a_nn_view.modified())
          {
            data().jacobian_nn_view.bind(lfsv_n_cache,lfsu_n_cache);
            localAssembler().etadd(data().a_nn,data().jacobian_nn_view);
            data().jacobian_nn_view.commit();
          }
        data().a_nn_view.resetModified();
        if (data().a_sn_view.modified())
          {
            data().jacobian_sn_view.bind(lfsv_s_cache,lfsu_n_cache);
            localAssembler().etadd(data().a_sn,data().jacobian_sn_view);
            data().jacobian_sn_view.commit();
          }
        data().a_sn_view.resetModified();
        if (data().a_ns_view.modified())
          {
            data().jacobian_ns_view.bind(lfsv_n_cache,lfsu_s_cache);
            localAssembler().etadd(data().a_ns,data().jacobian_ns_view);
            data().jacobian_ns_view.commit();
          }
        data().a_ns_view.resetModified();
      }
  }

  template<typename IG>
  void onUnbindLFSUVCoupling(const IG& ig,
                             const LFSU_Cache & lfsu_s_cache, const LFSV_Cache & lfsv_s_cache,
                             const LFSU_Cache & lfsu_n_cache, const LFSV_Cache & lfsv_n_cache,
                             const LFSU_C_Cache & lfsu_c_cache, const LFSV_C_Cache & lfsv_c_cache)
  {
    // write back local jacobian contributions
    if (localAssembler().writeData())
      {
        if (data().a_cc_view.modified())
          {
            data().jacobian_cc_view.bind(lfsv_c_cache,lfsu_c_cache);
            localAssembler().etadd(data().a_cc,data().jacobian_cc_view);
            data().jacobian_cc_view.commit();
          }
        data().a_cc_view.resetModified();
        if (data().a_sc_view.modified())
          {
            data().jacobian_sc_view.bind(lfsv_s_cache,lfsu_c_cache);
            localAssembler().etadd(data().a_sc,data().jacobian_sc_view);
            data().jacobian_sc_view.commit();
          }
        data().a_sc_view.resetModified();
        if (data().a_cs_view.modified())
          {
            data().jacobian_cs_view.bind(lfsv_c_cache,lfsu_s_cache);
            localAssembler().etadd(data().a_cs,data().jacobian_cs_view);
            data().jacobian_cs_view.commit();
          }
        data().a_cs_view.resetModified();
        if (data().a_nc_view.modified())
          {
            data().jacobian_nc_view.bind(lfsv_n_cache,lfsu_c_cache);
            localAssembler().etadd(data().a_nc,data().jacobian_nc_view);
            data().jacobian_nc_view.commit();
          }
        data().a_nc_view.resetModified();
        if (data().a_cn_view.modified())
          {
            data().jacobian_cn_view.bind(lfsv_c_cache,lfsu_n_cache);
            localAssembler().etadd(data().a_cn,data().jacobian_cn_view);
            data().jacobian_cn_view.commit();
          }
        data().a_cn_view.resetModified();
      }
  }


  template<typename EG>
  void assembleUVVolume(const EG& eg, const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache)
  {
    typedef visitor<functors::invoke_jacobian_volume,do_alpha_volume<> > Visitor;
    data().a_ss_view.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_eg(eg),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_x(data().x_s),wrap_a_ss(data().a_ss_view)
      )
    );
  }

  template<typename IG>
  void assembleUVSkeleton(const IG& ig,
                          const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache,
                          const LFSU_Cache& lfsu_n_cache, const LFSV_Cache& lfsv_n_cache)
  {
    data().a_ss_view.setWeight(localAssembler().weight());
    data().a_sn_view.setWeight(localAssembler().weight());
    data().a_ns_view.setWeight(localAssembler().weight());
    data().a_nn_view.setWeight(localAssembler().weight());
    typedef visitor<functors::invoke_jacobian_skeleton_or_boundary,do_alpha_skeleton_or_boundary<> > SubProblemVisitor;
    localAssembler().applyToSubProblems(
      SubProblemVisitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_x_s(data().x_s),
        wrap_lfsu_n_cache(lfsu_n_cache),wrap_lfsv_n_cache(lfsv_n_cache),
        wrap_x_n(data().x_n),
        wrap_a_ss(data().a_ss_view),wrap_a_sn(data().a_sn_view),
        wrap_a_ns(data().a_ns_view),wrap_a_nn(data().a_nn_view)
      )
    );

    typedef visitor<functors::invoke_jacobian_coupling,do_alpha_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(
      CouplingVisitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_x_s(data().x_s),
        wrap_lfsu_n_cache(lfsu_n_cache),wrap_lfsv_n_cache(lfsv_n_cache),
        wrap_x_n(data().x_n),
        wrap_a_ss(data().a_ss_view),wrap_a_sn(data().a_sn_view),
        wrap_a_ns(data().a_ns_view),wrap_a_nn(data().a_nn_view)
      )
    );
  }

  template<typename IG>
  void assembleUVBoundary(const IG& ig, const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache)
  {
    typedef visitor<functors::invoke_jacobian_boundary,do_alpha_boundary<> > Visitor;
    data().a_ss_view.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_x(data().x_s),wrap_a(data().a_ss_view)));
  }

  template<typename IG>
  void assembleUVEnrichedCoupling(const IG& ig,
                                  const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache,
                                  const LFSU_Cache& lfsu_n_cache, const LFSV_Cache& lfsv_n_cache,
                                  const LFSU_C_Cache& lfsu_c_cache, const LFSV_C_Cache& lfsv_c_cache)
  {
    data().a_ss_view.setWeight(localAssembler().weight());
    data().a_sc_view.setWeight(localAssembler().weight());
    data().a_nn_view.setWeight(localAssembler().weight());
    data().a_nc_view.setWeight(localAssembler().weight());
    data().a_cc_view.setWeight(localAssembler().weight());
    data().a_cs_view.setWeight(localAssembler().weight());
    data().a_cn_view.setWeight(localAssembler().weight());
    typedef visitor<functors::invoke_jacobian_enriched_coupling,do_alpha_enriched_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(
      CouplingVisitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_x_s(data().x_s),
        wrap_lfsu_n_cache(lfsu_n_cache),wrap_lfsv_n_cache(lfsv_n_cache),
        wrap_x_n(data().x_n),
        wrap_lfsu_c_cache(lfsu_c_cache),wrap_lfsv_c_cache(lfsv_c_cache),
        wrap_x_c(data().x_c),
        wrap_a_ss(data().a_ss_view),wrap_a_sc(data().a_sc_view),
        wrap_a_nn(data().a_nn_view),wrap_a_nc(data().a_nc_view),
        wrap_a_cc(data().a_cc_view),wrap_a_cs(data().a_cs_view),wrap_a_cn(data().a_cn_view)
      )
    );
  }

  template<typename EG>
  void assembleUVVolumePostSkeleton(const EG& eg, const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache)
  {
    typedef visitor<functors::invoke_jacobian_volume_post_skeleton,do_alpha_volume_post_skeleton<> > Visitor;
    data().a_ss_view.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_eg(eg),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_x(data().x_s),wrap_a(data().a_ss_view)
      )
    );
  }


  void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
  {
    if (localAssembler().writeData())
      {
        // save reference to container for constraints handling
        auto& jacobian = data().jacobian_ss_view.container();

        data().domain_s_view.detach();
        data().domain_n_view.detach();
        data().domain_c_view.detach();

        data().jacobian_ss_view.detach();
        data().jacobian_sn_view.detach();
        data().jacobian_sc_view.detach();
        data().jacobian_nn_view.detach();
        data().jacobian_ns_view.detach();
        data().jacobian_nc_view.detach();
        data().jacobian_cc_view.detach();
        data().jacobian_cs_view.detach();
        data().jacobian_cn_view.detach();

        localAssembler().handle_dirichlet_constraints(gfsv,jacobian);

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

  void setJacobian(Jacobian& a)
  {
    if (localAssembler().readData())
      {
        data().jacobian_ss_view.attach(a);
        data().jacobian_sn_view.attach(a);
        data().jacobian_sc_view.attach(a);
        data().jacobian_nn_view.attach(a);
        data().jacobian_ns_view.attach(a);
        data().jacobian_nc_view.attach(a);
        data().jacobian_cc_view.attach(a);
        data().jacobian_cs_view.attach(a);
        data().jacobian_cn_view.attach(a);
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

  JacobianAssemblerEngine(const LocalAssembler& localAssembler)
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

#endif // DUNE_PDELAB_MULTIDOMAIN_JACOBIANASSEMBLERENGINE_HH
