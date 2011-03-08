// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_RESIDUALASSEMBLERENGINE_HH
#define DUNE_PDELAB_MULTIDOMAIN_RESIDUALASSEMBLERENGINE_HH

#include <algorithm>

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>

#include <dune/pdelab/multidomain/datawrappers.hh>
#include <dune/pdelab/multidomain/residualassemblerfunctors.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


template<typename LA>
class ResidualAssemblerEngine
{

public:

  typedef LA LocalAssembler;
  typedef typename LA::Range Range;
  typedef typename LA::Domain Domain;

  bool requireIntersections() const
  {
    return requireUVSkeleton() || requireVSkeleton() ||
      requireUVEnrichedCoupling() || requireVEnrichedCoupling() ||
      requireUVBoundary() || requireVBoundary();
  }

  bool requireIntersectionsTwoSided() const
  {
    return requireUVBoundary() || requireVBoundary() ||
      any_child<typename LocalAssembler::SubProblems,do_skeleton_two_sided<> >::value;
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
    return any_child<typename LocalAssembler::SubProblems,do_alpha_skeleton<> >::value ||
      any_child<typename LocalAssembler::SubProblems,do_alpha_boundary<> >::value ||
      any_child<typename LocalAssembler::Couplings,do_alpha_coupling<> >::value ||
      any_child<typename LocalAssembler::Couplings,do_alpha_enriched_coupling<> >::value;
  }

  bool requireVSkeleton() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_lambda_skeleton<> >::value ||
      any_child<typename LocalAssembler::SubProblems,do_lambda_boundary<> >::value ||
      any_child<typename LocalAssembler::Couplings,do_lambda_coupling<> >::value ||
      any_child<typename LocalAssembler::Couplings,do_lambda_enriched_coupling<> >::value;
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

  template<typename EG, typename LFSU, typename LFSV>
  void onBindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    x_s.resize(lfsu.size());
  }

  template<typename EG, typename LFSV>
  void onBindLFSV(const EG& eg, const LFSV& lfsv)
  {
    // clear local residual
    _r_s.resize(lfsv.size());
    std::fill(_r_s.begin(),_r_s.end(),0.0);
  }

  template<typename EG, typename LFSU, typename LFSV>
  void onUnbindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {}

  template<typename EG, typename LFSV_S>
  void onUnbindLFSV(const EG& eg, const LFSV_S& lfsv_s)
  {
    // accumulate local residual into global residual
    if (r_s.modified())
      lfsv_s.vadd(_r_s,*r);
    r_s.resetModified();
  }


  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onBindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    x_n.resize(lfsu_n.size());
  }

  template<typename IG, typename LFSV_N>
  void onBindLFSVOutside(const IG& ig, const LFSV_N& lfsv_n)
  {
    // clear local residual
    _r_n.resize(lfsv_n.size());
    std::fill(_r_n.begin(),_r_n.end(),0.0);
  }

  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onUnbindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {}

  template<typename IG, typename LFSV_N>
  void onUnbindLFSVOutside(const IG& ig, const LFSV_N& lfsv_n)
  {
    // accumulate local residual into global residual
    if (r_n.modified())
      lfsv_n.vadd(_r_n,*r);
    r_n.resetModified();
  }


  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onBindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
    x_c.resize(lfsu_c.size());
  }

  template<typename IG, typename LFSV_C>
  void onBindLFSVCoupling(const IG& ig, const LFSV_C& lfsv_c)
  {
    // clear local residual
    _r_c.resize(lfsv_c.size());
    std::fill(_r_c.begin(),_r_c.end(),0.0);
  }

  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onUnbindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {}

  template<typename IG, typename LFSV_C>
  void onUnbindLFSVCoupling(const IG& ig, const LFSV_C& lfsv_c)
  {
    // accumulate local residual into global residual
    if (r_c.modified())
      lfsv_c.vadd(_r_c,*r);
    r_c.resetModified();
  }

  template<typename LFSU>
  void loadCoefficientsLFSUInside(const LFSU& lfsu_s)
  {
    // read local data
    lfsu_s.vread(*x,x_s);
  }

  template<typename LFSU_N>
  void loadCoefficientsLFSUOutside(const LFSU_N& lfsu_n)
  {
    // read local data
    lfsu_n.vread(*x,x_n);
  }

  template<typename LFSU_C>
  void loadCoefficientsLFSUCoupling(const LFSU_C& lfsu_c)
  {
    // read local data
    lfsu_c.vread(*x,x_c);
  }


  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolume(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<functors::alpha_volume,do_alpha_volume<> > Visitor;
    r_s.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(x_s),wrap_r(r_s)));
  }

  template<typename EG, typename LFSV>
  void assembleVVolume(const EG& eg, const LFSV& lfsv)
  {
    typedef visitor<functors::lambda_volume,do_lambda_volume<> > Visitor;
    r_s.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_eg(eg),
                                                          wrap_lfsv(lfsv),wrap_r(r_s)));
  }


  template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
  void assembleUVSkeleton(const IG& ig,
                          const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                          const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    typedef visitor<functors::alpha_skeleton_or_boundary,do_alpha_skeleton_or_boundary<> > SubProblemVisitor;
    r_s.setWeight(localAssembler().weight());
    r_n.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(SubProblemVisitor::add_data(wrap_operator_type(SpatialOperator()),wrap_ig(ig),
                                                                    wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(x_s),wrap_r_s(r_s),
                                                                    wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(x_n),wrap_r_n(r_n)));

    typedef visitor<functors::alpha_coupling,do_alpha_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(CouplingOperator()),wrap_ig(ig),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(x_s),wrap_r_s(r_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(x_n),wrap_r_n(r_n)));
  }

  template<typename IG, typename LFSV_S, typename LFSV_N>
  void assembleVSkeleton(const IG& ig,
                         const LFSV_S& lfsv_s,
                         const LFSV_N& lfsv_n)
  {
    typedef visitor<functors::lambda_skeleton_or_boundary,do_lambda_skeleton_or_boundary<> > Visitor;
    r_s.setWeight(localAssembler().weight());
    r_n.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_ig(ig),
                                                          wrap_lfsv_s(lfsv_s),wrap_r_s(r_s),
                                                          wrap_lfsv_n(lfsv_n),wrap_r_n(r_n)));

    typedef visitor<functors::lambda_coupling,do_lambda_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(CouplingOperator()),wrap_ig(ig),
                                                                wrap_lfsv_s(lfsv_s),wrap_r_s(r_s),
                                                                wrap_lfsv_n(lfsv_n),wrap_r_n(r_n)));
  }


  template<typename IG, typename LFSU, typename LFSV>
  void assembleUVBoundary(const IG& ig, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<functors::alpha_boundary,do_alpha_boundary<> > Visitor;
    r_s.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_ig(ig),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(x_s),wrap_r(r_s)));
  }

  template<typename IG, typename LFSV>
  void assembleVBoundary(const IG& ig, const LFSV& lfsv)
  {
    typedef visitor<functors::lambda_boundary,do_lambda_boundary<> > Visitor;
    r_s.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_ig(ig),
                                                          wrap_lfsv(lfsv),wrap_r(r_s)));
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
    typedef visitor<functors::alpha_enriched_coupling,do_alpha_enriched_coupling<> > CouplingVisitor;
    r_s.setWeight(localAssembler().weight());
    r_n.setWeight(localAssembler().weight());
    r_c.setWeight(localAssembler().weight());
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(CouplingOperator()),wrap_ig(ig),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(x_s),wrap_r_s(r_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(x_n),wrap_r_n(r_n),
                                                                wrap_lfsu_c(lfsu_c),wrap_lfsv_c(lfsv_c),wrap_x_c(x_c),wrap_r_c(r_c)));
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
    typedef visitor<functors::lambda_enriched_coupling,do_lambda_enriched_coupling<> > CouplingVisitor;
    r_s.setWeight(localAssembler().weight());
    r_n.setWeight(localAssembler().weight());
    r_c.setWeight(localAssembler().weight());
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(CouplingOperator()),wrap_ig(ig),
                                                                wrap_lfsv_s(lfsv_s),wrap_r_s(r_s),
                                                                wrap_lfsv_n(lfsv_n),wrap_r_n(r_n),
                                                                wrap_lfsv_c(lfsv_c),wrap_r_c(r_c)));
  }


  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolumePostSkeleton(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<functors::alpha_volume_post_skeleton,do_alpha_volume_post_skeleton<> > Visitor;
    r_s.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(x_s),wrap_r(r_s)));
  }

  template<typename EG, typename LFSV>
  void assembleVVolumePostSkeleton(const EG& eg, const LFSV& lfsv)
  {
    typedef visitor<functors::lambda_volume_post_skeleton,do_lambda_volume_post_skeleton<> > Visitor;
    r_s.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_eg(eg),
                                                          wrap_lfsv(lfsv),wrap_r(r_s)));
  }


  void preAssembly()
  {}

  void postAssembly()
  {
    Dune::PDELab::constrain_residual(localAssembler().testConstraints(),*r);
  }

  void setSolution(const Domain& x_)
  {
    x = &x_;
  }

  void setResidual(Range& r_)
  {
    r = &r_;
  }

  const LocalAssembler& localAssembler() const
  {
    return _localAssembler;
  }

  ResidualAssemblerEngine(const LocalAssembler& localAssembler)
    : _localAssembler(localAssembler)
    , x(nullptr)
    , r(nullptr)
    , r_s(_r_s.weightedAccumulationView(1.0))
    , r_n(_r_n.weightedAccumulationView(1.0))
    , r_c(_r_c.weightedAccumulationView(1.0))
  {}

private:

  const LocalAssembler& _localAssembler;

  const Domain* x;
  Range* r;

  typedef LocalVector<typename Domain::ElementType, TrialSpaceTag> SolutionVector;
  typedef LocalVector<typename Range::ElementType, TestSpaceTag> ResidualVector;

  SolutionVector x_s;
  ResidualVector _r_s;
  SolutionVector x_n;
  ResidualVector _r_n;
  SolutionVector x_c;
  ResidualVector _r_c;

  typename ResidualVector::WeightedAccumulationView r_s;
  typename ResidualVector::WeightedAccumulationView r_n;
  typename ResidualVector::WeightedAccumulationView r_c;

};



} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_RESIDUALASSEMBLERENGINE_HH
