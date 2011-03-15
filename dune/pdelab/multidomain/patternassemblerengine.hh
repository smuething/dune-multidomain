// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_PATTERNASSEMBLERENGINE_HH
#define DUNE_PDELAB_MULTIDOMAIN_PATTERNASSEMBLERENGINE_HH

#include <algorithm>

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>
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

  typedef LA LocalAssembler;
  typedef typename LA::Domain Domain;
  typedef typename LA::Jacobian Jacobian;
  typedef typename LA::Pattern GlobalPattern;

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


  template<typename EG, typename LFSU, typename LFSV>
  void onBindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    if (localAssembler().readData())
      data().pattern_ss.clear();
  }

  template<typename EG, typename LFSU, typename LFSV>
  void onUnbindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    if (localAssembler().writeData())
      addToGlobalPattern(data().pattern_ss,lfsv,lfsu);
  }


  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onBindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    if (localAssembler().readData())
      {
        data().pattern_sn.clear();
        data().pattern_ns.clear();
        data().pattern_nn.clear();
      }
  }

  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onUnbindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    if (localAssembler().writeData())
      addToGlobalPattern(data().pattern_nn,lfsv_n,lfsu_n);
  }


  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onBindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
    if (localAssembler().readData())
      {
        data().pattern_cc.clear();
        data().pattern_sc.clear();
        data().pattern_nc.clear();
        data().pattern_cs.clear();
        data().pattern_cn.clear();
      }
  }

  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onUnbindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
    if (localAssembler().writeData())
      addToGlobalPattern(data().pattern_cc,lfsv_c,lfsu_c);
  }


  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolume(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<functors::invoke_pattern_volume,do_pattern_volume<> > Visitor;
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(DefaultOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),
                                                          wrap_pattern_ss(data().pattern_ss)));
  }

  template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
  void assembleUVSkeleton(const IG& ig,
                          const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                          const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    typedef visitor<functors::invoke_pattern_skeleton_or_boundary,do_pattern_skeleton_or_boundary<> > SubProblemVisitor;

    localAssembler().applyToSubProblems(SubProblemVisitor::add_data(wrap_operator_type(DefaultOperator()),wrap_ig(ig),
                                                                    wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),
                                                                    wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),
                                                                    wrap_pattern_ss(data().pattern_ss),
                                                                    wrap_pattern_sn(data().pattern_sn),
                                                                    wrap_pattern_nn(data().pattern_nn),
                                                                    wrap_pattern_ns(data().pattern_ns)));

    typedef visitor<functors::invoke_pattern_coupling,do_pattern_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(DefaultOperator()),wrap_ig(ig),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),
                                                                wrap_pattern_ss(data().pattern_ss),
                                                                wrap_pattern_sn(data().pattern_sn),
                                                                wrap_pattern_nn(data().pattern_nn),
                                                                wrap_pattern_ns(data().pattern_ns)));

    if (localAssembler().writeData())
      {
        addToGlobalPattern(data().pattern_sn,lfsv_s,lfsu_n);
        addToGlobalPattern(data().pattern_ns,lfsv_n,lfsu_s);
      }
  }

  template<typename IG, typename LFSU, typename LFSV>
  void assembleUVBoundary(const IG& ig, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<functors::invoke_pattern_boundary,do_pattern_boundary<> > Visitor;
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(DefaultOperator()),wrap_ig(ig),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_pattern_ss(data().pattern_ss)));
  }

  template<typename IG,
           typename LFSU_S, typename LFSV_S,
           typename LFSU_N, typename LFSV_N,
           typename LFSU_C, typename LFSV_C>
  void assembleUVEnrichedCoupling(const IG& ig,
                                  const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                                  const LFSU_N& lfsu_n, const LFSV_N& lfsv_n,
                                  const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
    typedef visitor<functors::invoke_pattern_enriched_coupling,do_pattern_enriched_coupling<> > CouplingVisitor;

    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(DefaultOperator()),wrap_ig(ig),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),
                                                                wrap_lfsu_c(lfsu_c),wrap_lfsv_c(lfsv_c),
                                                                wrap_pattern_ss(data().pattern_ss),
                                                                wrap_pattern_sc(data().pattern_sc),
                                                                wrap_pattern_nn(data().pattern_nn),
                                                                wrap_pattern_nc(data().pattern_nc),
                                                                wrap_pattern_cc(data().pattern_cc),
                                                                wrap_pattern_cs(data().pattern_cs),
                                                                wrap_pattern_cn(data().pattern_cn)));

    if (localAssembler().writeData())
      {
        addToGlobalPattern(data().pattern_sc,lfsv_s,lfsu_c);
        addToGlobalPattern(data().pattern_cs,lfsv_c,lfsu_s);
        addToGlobalPattern(data().pattern_nc,lfsv_n,lfsu_c);
        addToGlobalPattern(data().pattern_cn,lfsv_c,lfsu_n);
      }
  }

  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolumePostSkeleton(const EG& eg, const LFSV& lfsv, const LFSU& lfsu)
  {
    typedef visitor<functors::invoke_pattern_volume_post_skeleton,do_pattern_volume_post_skeleton<> > Visitor;
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(DefaultOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),
                                                          wrap_pattern_ss(data().pattern_ss)));
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

  template<typename LFSU, typename LFSV>
  void addToGlobalPattern(const LocalSparsityPattern& pattern, const LFSV& lfsv, const LFSU& lfsu)
  {
    for(auto it = pattern.begin(); it != pattern.end(); ++it)
      localAssembler().add_entry(*data()._globalPattern,lfsv.globalIndex(it->i()),lfsu.globalIndex(it->j()));
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
