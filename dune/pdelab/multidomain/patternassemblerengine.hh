// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_PATTERNASSEMBLERENGINE_HH
#define DUNE_PDELAB_MULTIDOMAIN_PATTERNASSEMBLERENGINE_HH

#include <algorithm>

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>
#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>

#include <dune/pdelab/multidomain/datawrappers.hh>
#include <dune/pdelab/multidomain/patternassemblerfunctors.hh>


namespace Dune {

namespace PDELab {

namespace MultiDomain {

namespace {

  template<typename GlobalPattern>
  struct PatternAssemblerEngineData
  {

    GlobalPattern* _globalPattern;

    LocalSparsityPattern pattern_ss;
    LocalSparsityPattern pattern_sn;
    LocalSparsityPattern pattern_sc;

    LocalSparsityPattern pattern_nn;
    LocalSparsityPattern pattern_ns;
    LocalSparsityPattern pattern_nc;

    LocalSparsityPattern pattern_cc;
    LocalSparsityPattern pattern_cs;
    LocalSparsityPattern pattern_cn;

    PatternAssemblerEngineData()
      : _globalPattern(nullptr)
    {}

  };

} // anonymous namespace

template<typename LA>
class PatternAssemblerEngine
  : public LocalAssemblerEngineBase
{

  template<typename>
  friend class PatternAssemblerEngine;

public:

  static const bool needs_constraints_caching = true;

  typedef LA LocalAssembler;
  typedef typename LA::Domain Domain;
  typedef typename LA::Jacobian Jacobian;
  typedef typename LA::Pattern GlobalPattern;

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

  typedef PatternAssemblerEngineData<GlobalPattern> EngineData;

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
      any_child<typename LocalAssembler::Couplings,do_pattern_coupling<> >::value ||
      any_child<LocalAssembler,do_skeleton_two_sided<> >::value;
  }

  bool requireUVVolume() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_pattern_volume<> >::value;
  }

  bool requireUVSkeleton() const
  {
    return
      any_child<typename LocalAssembler::SubProblems,do_pattern_skeleton<> >::value ||
      any_child<typename LocalAssembler::SubProblems,do_pattern_boundary<> >::value ||
      any_child<typename LocalAssembler::Couplings,do_pattern_coupling<> >::value;
  }

  bool requireUVBoundary() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_pattern_boundary<> >::value;
  }

  bool requireUVEnrichedCoupling() const
  {
    return any_child<typename LocalAssembler::Couplings,do_pattern_enriched_coupling<> >::value;
  }

  bool requireUVVolumePostSkeleton() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_pattern_volume_post_skeleton<> >::value;
  }


  template<typename EG>
  void onUnbindLFSUV(const EG& eg, const LFSU_Cache & lfsu_s_cache, const LFSV_Cache & lfsv_s_cache)
  {
    if (localAssembler().writeData())
      {
        addToGlobalPattern(data().pattern_ss,lfsv_s_cache,lfsu_s_cache);
        data().pattern_ss.clear();
      }
  }

  template<typename IG>
  void onUnbindLFSUVOutside(const IG& ig,
                            const LFSU_Cache & lfsu_s_cache, const LFSV_Cache & lfsv_s_cache,
                            const LFSU_Cache & lfsu_n_cache, const LFSV_Cache & lfsv_n_cache)
  {
    if (localAssembler().writeData())
      {
        addToGlobalPattern(data().pattern_nn,lfsv_n_cache,lfsu_n_cache);
        data().pattern_nn.clear();
        addToGlobalPattern(data().pattern_sn,lfsv_s_cache,lfsu_n_cache);
        data().pattern_sn.clear();
        addToGlobalPattern(data().pattern_ns,lfsv_n_cache,lfsu_s_cache);
        data().pattern_ns.clear();
      }
  }

  template<typename IG>
  void onUnbindLFSUVCoupling(const IG& ig,
                             const LFSU_Cache & lfsu_s_cache, const LFSV_Cache & lfsv_s_cache,
                             const LFSU_Cache & lfsu_n_cache, const LFSV_Cache & lfsv_n_cache,
                             const LFSU_C_Cache & lfsu_c_cache, const LFSV_C_Cache & lfsv_c_cache)
  {
    if (localAssembler().writeData())
      {
        addToGlobalPattern(data().pattern_cc,lfsv_c_cache,lfsu_c_cache);
        data().pattern_cc.clear();
        addToGlobalPattern(data().pattern_sc,lfsv_s_cache,lfsu_c_cache);
        data().pattern_sc.clear();
        addToGlobalPattern(data().pattern_cs,lfsv_c_cache,lfsu_s_cache);
        data().pattern_cs.clear();
        addToGlobalPattern(data().pattern_nc,lfsv_n_cache,lfsu_c_cache);
        data().pattern_nc.clear();
        addToGlobalPattern(data().pattern_cn,lfsv_c_cache,lfsu_n_cache);
        data().pattern_cn.clear();
      }
  }


  template<typename EG>
  void assembleUVVolume(const EG& eg, const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache)
  {
    typedef visitor<functors::invoke_pattern_volume,do_pattern_volume<> > Visitor;
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_eg(eg),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_pattern_ss(data().pattern_ss)
      )
    );
  }

  template<typename IG>
  void assembleUVSkeleton(const IG& ig,
                          const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache,
                          const LFSU_Cache& lfsu_n_cache, const LFSV_Cache& lfsv_n_cache)
  {
    typedef visitor<functors::invoke_pattern_skeleton_or_boundary,do_pattern_skeleton_or_boundary<> > SubProblemVisitor;

    localAssembler().applyToSubProblems(
      SubProblemVisitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_lfsu_n_cache(lfsu_n_cache),wrap_lfsv_n_cache(lfsv_n_cache),
        wrap_pattern_ss(data().pattern_ss),
        wrap_pattern_sn(data().pattern_sn),
        wrap_pattern_nn(data().pattern_nn),
        wrap_pattern_ns(data().pattern_ns)
      )
    );

    typedef visitor<functors::invoke_pattern_coupling,do_pattern_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(
      CouplingVisitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_lfsu_n_cache(lfsu_n_cache),wrap_lfsv_n_cache(lfsv_n_cache),
        wrap_pattern_ss(data().pattern_ss),
        wrap_pattern_sn(data().pattern_sn),
        wrap_pattern_nn(data().pattern_nn),
        wrap_pattern_ns(data().pattern_ns)
      )
    );
  }

  template<typename IG>
  void assembleUVBoundary(const IG& ig, const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache)
  {
    typedef visitor<functors::invoke_pattern_boundary,do_pattern_boundary<> > Visitor;
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_pattern_ss(data().pattern_ss)
      )
    );
  }

  template<typename IG>
  void assembleUVEnrichedCoupling(const IG& ig,
                                  const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache,
                                  const LFSU_Cache& lfsu_n_cache, const LFSV_Cache& lfsv_n_cache,
                                  const LFSU_C_Cache& lfsu_c_cache, const LFSV_C_Cache& lfsv_c_cache)
  {
    typedef visitor<functors::invoke_pattern_enriched_coupling,do_pattern_enriched_coupling<> > CouplingVisitor;

    localAssembler().applyToCouplings(
      CouplingVisitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_ig(ig),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_lfsu_n_cache(lfsu_n_cache),wrap_lfsv_n_cache(lfsv_n_cache),
        wrap_lfsu_c_cache(lfsu_c_cache),wrap_lfsv_c_cache(lfsv_c_cache),
        wrap_pattern_ss(data().pattern_ss),
        wrap_pattern_sc(data().pattern_sc),
        wrap_pattern_nn(data().pattern_nn),
        wrap_pattern_nc(data().pattern_nc),
        wrap_pattern_cc(data().pattern_cc),
        wrap_pattern_cs(data().pattern_cs),
        wrap_pattern_cn(data().pattern_cn)
      )
    );
  }

  template<typename EG>
  void assembleUVVolumePostSkeleton(const EG& eg, const LFSU_Cache& lfsu_s_cache, const LFSV_Cache& lfsv_s_cache)
  {
    typedef visitor<functors::invoke_pattern_volume_post_skeleton,do_pattern_volume_post_skeleton<> > Visitor;
    localAssembler().applyToSubProblems(
      Visitor::add_data(
        wrap_operator_type(DefaultOperator()),
        wrap_eg(eg),
        wrap_lfsu_s_cache(lfsu_s_cache),wrap_lfsv_s_cache(lfsv_s_cache),
        wrap_pattern_ss(data().pattern_ss)
      )
    );
  }


  void setPattern(GlobalPattern& globalPattern)
  {
    if (localAssembler().writeData())
      data()._globalPattern = &globalPattern;
  }

  const LocalAssembler& localAssembler() const
  {
    return _localAssembler;
  }

  const typename Traits::TrialGridFunctionSpaceConstraints& trialGridFunctionSpaceConstraints() const
  {
    return localAssembler().trialConstraints();
  }

  const typename Traits::TestGridFunctionSpaceConstraints& testGridFunctionSpaceConstraints() const
  {
    return localAssembler().testConstraints();
  }

  PatternAssemblerEngine(const LocalAssembler& localAssembler)
    : _localAssembler(localAssembler)
    , _engineData(make_shared<EngineData>())
  {}

  template<typename Engine>
  void shareData(Engine& engine)
  {
    _engineData = engine._engineData;
  }

private:

  template<typename LFSV_Cache_, typename LFSU_Cache_>
  void addToGlobalPattern(const LocalSparsityPattern& pattern, const LFSV_Cache_& lfsv_cache, const LFSU_Cache_& lfsu_cache)
  {
    for(size_t k = 0; k < pattern.size(); ++k)
      localAssembler().add_entry(
        *data()._globalPattern,
        lfsv_cache,pattern[k].i(),
        lfsu_cache,pattern[k].j()
      );
  }

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

#endif // DUNE_PDELAB_MULTIDOMAIN_PATTERNASSEMBLERENGINE_HH
