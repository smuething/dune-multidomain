// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_PATTERNASSEMBLERENGINE_HH
#define DUNE_PDELAB_MULTIDOMAIN_PATTERNASSEMBLERENGINE_HH

#include <algorithm>

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>

#include <dune/pdelab/multidomain/datawrappers.hh>
#include <dune/pdelab/multidomain/patternassemblerfunctors.hh>


namespace Dune {

namespace PDELab {

namespace MultiDomain {


template<typename LA>
class PatternAssemblerEngine
{

public:

  typedef LA LocalAssembler;
  typedef typename LA::Domain Domain;
  typedef typename LA::Jacobian Jacobian;
  typedef typename LA::Pattern GlobalPattern;


  bool requireIntersections() const
  {
    return requireUVSkeleton() || requireUVEnrichedCoupling() || requireUVBoundary();
  }

  bool requireIntersectionsTwoSided() const
  {
    return requireUVBoundary() ||
      requireUVEnrichedCoupling() || any_child<typename LocalAssembler::Couplings,do_pattern_coupling<> >::value ||
      any_child<typename LocalAssembler::SubProblems,do_skeleton_two_sided<> >::value;
  }

  bool requireUVVolume() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_pattern_volume<> >::value;
  }

  bool requireVVolume() const
  {
    return false;
  }

  bool requireUVSkeleton() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_pattern_skeleton<> >::value ||
      any_child<typename LocalAssembler::SubProblems,do_pattern_boundary<> >::value ||
      any_child<typename LocalAssembler::Couplings,do_pattern_coupling<> >::value;
  }

  bool requireVSkeleton() const
  {
    return false;
  }

  bool requireUVBoundary() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_pattern_boundary<> >::value;
  }

  bool requireVBoundary() const
  {
    return false;
  }

  bool requireUVEnrichedCoupling() const
  {
    return any_child<typename LocalAssembler::Couplings,do_pattern_enriched_coupling<> >::value;
  }

  bool requireVEnrichedCoupling() const
  {
    return false;
  }

  bool requireUVVolumePostSkeleton() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_pattern_volume_post_skeleton<> >::value;
  }

  bool requireVVolumePostSkeleton() const
  {
    return false;
  }

  template<typename EG, typename LFSU, typename LFSV>
  void onBindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    pattern_ss.clear();
  }

  template<typename EG, typename LFSU, typename LFSV>
  void onUnbindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    addToGlobalPattern(pattern_ss,lfsu,lfsv);
  }

  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onBindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    pattern_sn.clear();
    pattern_ns.clear();
  }

  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onUnbindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
  }

  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onBindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
  }

  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onUnbindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
  }

  template<typename EG, typename LFSV>
  void onBindLFSV(const EG& eg, const LFSV& lfsv)
  {
  }

  template<typename EG, typename LFSV>
  void onUnbindLFSV(const EG& eg, const LFSV& lfsv)
  {
  }

  template<typename IG, typename LFSV_N>
  void onBindLFSVOutside(const IG& ig, const LFSV_N& lfsv_n)
  {
  }

  template<typename IG, typename LFSV_N>
  void onUnbindLFSVOutside(const IG& ig, const LFSV_N& lfsv_n)
  {
  }

  template<typename IG, typename LFSV_C>
  void onBindLFSVCoupling(const IG& ig, const LFSV_C& lfsv_c)
  {
  }

  template<typename IG, typename LFSV_C>
  void onUnbindLFSVCoupling(const IG& ig, const LFSV_C& lfsv_c)
  {
  }

  template<typename LFSU>
  void loadCoefficientsLFSUInside(const LFSU& lfsu_s)
  {
  }

  template<typename LFSU_N>
  void loadCoefficientsLFSUOutside(const LFSU_N& lfsu_n)
  {
  }

  template<typename LFSU_C>
  void loadCoefficientsLFSUCoupling(const LFSU_C& lfsu_c)
  {
  }

  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolume(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<invoke_pattern_volume,do_pattern_volume<> > Visitor;
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),
                                                          wrap_pattern_ss(pattern_ss)));
  }

  template<typename EG, typename LFSV>
  void assembleVVolume(const EG& eg, const LFSV& lfsv)
  {
  }


  template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
  void assembleUVSkeleton(const IG& ig,
                          const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                          const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    typedef visitor<invoke_pattern_skeleton_or_boundary,do_pattern_skeleton_or_boundary<> > SubProblemVisitor;

    pattern_sn.clear();
    pattern_ns.clear();

    localAssembler().applyToSubProblems(SubProblemVisitor::add_data(wrap_operator_type(SpatialOperator()),wrap_ig(ig),
                                                                    store_neighbor_accessed(false),
                                                                    wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),
                                                                    wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),
                                                                    wrap_pattern_ss(pattern_ss),
                                                                    wrap_pattern_sn(pattern_sn),
                                                                    wrap_pattern_ns(pattern_ns)));

    typedef visitor<invoke_pattern_coupling,do_pattern_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(CouplingOperator()),wrap_ig(ig),
                                                                store_neighbor_accessed(false),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),
                                                                wrap_pattern_sn(pattern_sn),
                                                                wrap_pattern_ns(pattern_ns)));

    addToGlobalPattern(pattern_sn,lfsu_s,lfsv_n);
    addToGlobalPattern(pattern_ns,lfsu_n,lfsv_s);
  }

  template<typename IG, typename LFSV_S, typename LFSV_N>
  void assembleVSkeleton(const IG& ig,
                         const LFSV_S& lfsv_s,
                         const LFSV_N& lfsv_n)
  {
  }


  template<typename IG, typename LFSU, typename LFSV>
  void assembleUVBoundary(const IG& ig, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<invoke_pattern_boundary,do_pattern_boundary<> > Visitor;
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_ig(ig),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_pattern_ss(pattern_ss)));
  }

  template<typename IG, typename LFSV>
  void assembleVBoundary(const IG& ig, const LFSV& lfsv)
  {
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
    typedef visitor<invoke_pattern_enriched_coupling,do_pattern_enriched_coupling<> > CouplingVisitor;

    pattern_sc.clear();
    pattern_cs.clear();
    pattern_nc.clear();
    pattern_cn.clear();

    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(CouplingOperator()),wrap_ig(ig),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),
                                                                wrap_lfsu_c(lfsu_c),wrap_lfsv_c(lfsv_c),
                                                                wrap_pattern_sc(pattern_sc),
                                                                wrap_pattern_cs(pattern_cs),
                                                                wrap_pattern_nc(pattern_nc),
                                                                wrap_pattern_cn(pattern_cn)));
    addToGlobalPattern(pattern_sc,lfsu_s,lfsv_c);
    addToGlobalPattern(pattern_cs,lfsu_c,lfsv_s);
    addToGlobalPattern(pattern_nc,lfsu_n,lfsv_c);
    addToGlobalPattern(pattern_cn,lfsu_c,lfsv_n);
  }

  template<typename IG,
           typename LFSV_S,
           typename LFSV_N,
           typename LFSV_C>
  void assembleVEnrichedCoupling(const IG& ig,
                                 const LFSV_S& lfsv_s,
                                 const LFSV_N& lfsv_n,
                                 const LFSV_C& lfsv_c)
  {
  }


  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolumePostSkeleton(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<invoke_pattern_volume_post_skeleton,do_pattern_volume_post_skeleton<> > Visitor;
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),
                                                          wrap_pattern_ss(pattern_ss)));
  }

  template<typename EG, typename LFSV>
  void assembleVVolumePostSkeleton(const EG& eg, const LFSV& lfsv)
  {
  }


  void preAssembly()
  {}

  void postAssembly()
  {}

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
  void addToGlobalPattern(const LocalSparsityPattern pattern, const LFSU& lfsu, const LFSV& lfsv)
  {
    for(auto it = pattern.begin(); it != pattern.end(); ++it)
      localAssembler().add_entry(*_globalPattern,lfsv.globalIndex(it->i()),lfsu.globalIndex(it->j()));
  }

  GlobalPattern* _globalPattern;

  const LocalAssembler& _localAssembler;

  LocalSparsityPattern pattern_ss;

  LocalSparsityPattern pattern_sn;
  LocalSparsityPattern pattern_ns;

  LocalSparsityPattern pattern_sc;
  LocalSparsityPattern pattern_cs;

  LocalSparsityPattern pattern_nc;
  LocalSparsityPattern pattern_cn;

};



} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_PATTERNASSEMBLERENGINE_HH
