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


template<typename LA>
class PatternAssemblerEngine
  : public LocalAssemblerEngineBase
{

public:

  typedef LA LocalAssembler;
  typedef typename LA::Domain Domain;
  typedef typename LA::Jacobian Jacobian;
  typedef typename LA::Pattern GlobalPattern;


  bool requireIntersections() const
  {
    return
      requireUVSkeleton() ||
      requireUVEnrichedCoupling() ||
      requireUVBoundary();
  }

  bool requireIntersectionsTwoSided() const
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
    pattern_ss.clear();
  }

  template<typename EG, typename LFSU, typename LFSV>
  void onUnbindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    addToGlobalPattern(pattern_ss,lfsv,lfsu);
  }


  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onBindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    pattern_sn.clear();
    pattern_ns.clear();
    pattern_nn.clear();
  }

  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onUnbindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    addToGlobalPattern(pattern_nn,lfsv_n,lfsu_n);
  }


  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onBindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
    pattern_cc.clear();
    pattern_sc.clear();
    pattern_nc.clear();
    pattern_cs.clear();
    pattern_cn.clear();
  }

  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onUnbindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
    addToGlobalPattern(pattern_cc,lfsv_c,lfsu_c);
  }


  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolume(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<functors::invoke_pattern_volume,do_pattern_volume<> > Visitor;
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(DefaultOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),
                                                          wrap_pattern_ss(pattern_ss)));
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
                                                                    wrap_pattern_ss(pattern_ss),
                                                                    wrap_pattern_sn(pattern_sn),
                                                                    wrap_pattern_nn(pattern_nn),
                                                                    wrap_pattern_ns(pattern_ns)));

    typedef visitor<functors::invoke_pattern_coupling,do_pattern_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(DefaultOperator()),wrap_ig(ig),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),
                                                                wrap_pattern_ss(pattern_ss),
                                                                wrap_pattern_sn(pattern_sn),
                                                                wrap_pattern_nn(pattern_nn),
                                                                wrap_pattern_ns(pattern_ns)));

    addToGlobalPattern(pattern_sn,lfsv_s,lfsu_n);
    addToGlobalPattern(pattern_ns,lfsv_n,lfsu_s);
  }

  template<typename IG, typename LFSU, typename LFSV>
  void assembleUVBoundary(const IG& ig, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<functors::invoke_pattern_boundary,do_pattern_boundary<> > Visitor;
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(DefaultOperator()),wrap_ig(ig),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_pattern_ss(pattern_ss)));
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
                                                                wrap_pattern_ss(pattern_ss),
                                                                wrap_pattern_sc(pattern_sc),
                                                                wrap_pattern_nn(pattern_nn),
                                                                wrap_pattern_nc(pattern_nc),
                                                                wrap_pattern_cc(pattern_cc),
                                                                wrap_pattern_cs(pattern_cs),
                                                                wrap_pattern_cn(pattern_cn)));
    addToGlobalPattern(pattern_sc,lfsv_s,lfsu_c);
    addToGlobalPattern(pattern_cs,lfsv_c,lfsu_s);
    addToGlobalPattern(pattern_nc,lfsv_n,lfsu_c);
    addToGlobalPattern(pattern_cn,lfsv_c,lfsu_n);
  }

  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolumePostSkeleton(const EG& eg, const LFSV& lfsv, const LFSU& lfsu)
  {
    typedef visitor<functors::invoke_pattern_volume_post_skeleton,do_pattern_volume_post_skeleton<> > Visitor;
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(DefaultOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),
                                                          wrap_pattern_ss(pattern_ss)));
  }


  void setPattern(GlobalPattern& globalPattern)
  {
    _globalPattern = &globalPattern;
  }

  const LocalAssembler& localAssembler() const
  {
    return _localAssembler;
  }

  PatternAssemblerEngine(const LocalAssembler& localAssembler)
    : _localAssembler(localAssembler)
  {}

private:

  template<typename LFSU, typename LFSV>
  void addToGlobalPattern(const LocalSparsityPattern& pattern, const LFSV& lfsv, const LFSU& lfsu)
  {
    for(auto it = pattern.begin(); it != pattern.end(); ++it)
      localAssembler().add_entry(*_globalPattern,lfsv.globalIndex(it->i()),lfsu.globalIndex(it->j()));
  }

  GlobalPattern* _globalPattern;

  const LocalAssembler& _localAssembler;

  LocalSparsityPattern pattern_ss;
  LocalSparsityPattern pattern_sn;
  LocalSparsityPattern pattern_sc;

  LocalSparsityPattern pattern_nn;
  LocalSparsityPattern pattern_ns;
  LocalSparsityPattern pattern_nc;

  LocalSparsityPattern pattern_cc;
  LocalSparsityPattern pattern_cs;
  LocalSparsityPattern pattern_cn;

};



} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_PATTERNASSEMBLERENGINE_HH
