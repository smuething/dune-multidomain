// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_RESIDUALASSEMBLERENGINE_HH
#define DUNE_PDELAB_MULTIDOMAIN_RESIDUALASSEMBLERENGINE_HH

#include<map>
#include<tuple>

#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>

#include <dune/pdelab/common/geometrywrapper.hh>
//#include"../gridfunctionspace/gridfunctionspace.hh"
#include <dune/pdelab/gridfunctionspace/constraints.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>

#include <dune/pdelab/multidomain/multidomaingridoperatorspaceutilities.hh>
#include <dune/pdelab/multidomain/operatorapplier.hh>

#include <dune/pdelab/multidomain/multidomaingridoperatorspaceinvocationhelpers.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


template<typename X, typename R>
class ResidualAssemblerEngine
{

  template<typename EG, typename LFSU, typename LFSV>
  void onBindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    // read local data
    xl.resize(lfsu.size);
    lfsu.vread(x,xl);
  }

  template<typename EG, typename LFSV>
  void onBindLFSV(const EG& eg, const LFSV& lfsv)
  {
    // clear local residual
    rl.resize(lfsv.size);
    std::fill(rl.begin(),rl.end(),0.0);
  }

  template<typename EG, typename LFSU, typename LFSV>
  void onUnbindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {}

  template<typename EG, typename LFSV>
  void onUnbindLFSV(const EG& eg, const LFSV& lfsv)
  {
    // accumulate local residual into global residual
    lfsv.vadd(rl,r);
  }


  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onBindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    // read local data
    xn.resize(lfsu_n.size);
    lfsu_n.vread(x,xn);
  }

  template<typename IG, typename LFSV_N>
  void onBindLFSVOutside(const IG& ig, const LFSV_N& lfsv_n)
  {
    // clear local residual
    rn.resize(lfsv_n.size);
    std::fill(rn.begin(),rn.end(),0.0);
  }

  template<typename IG, typename LFSU_N>
  void onUnbindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, typename LFSV_N& lfsv_n)
  {}

  template<typename IG, typename LFSV_N>
  void onUnbindLFSVOutside(const IG& ig, const LFSV_N& lfsv_n)
  {
    // accumulate local residual into global residual
    lfsv_n.vadd(rn,r); // TODO: check if necessary
  }


  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onBindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
    // read local data
    xc.resize(lfsu_c.size);
    lfsu_c.vread(x,xc);
  }

  template<typename IG, typename LFSV_C>
  void onBindLFSVCoupling(const IG& ig, const LFSV_C& lfsv_c)
  {
    // clear local residual
    rc.resize(lfsv_c.size);
    std::fill(rc.begin(),rc.end(),0.0);
  }

  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onUnbindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {}

  template<typename IG, typename LFSV_C>
  void onUnbindLFSVCoupling(const IG& ig, const LFSV_C& lfsv_c)
  {
    // accumulate local residual into global residual
    lfsv_c.vadd(rc,r); // TODO: check if necessary
  }


  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolume(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<invoke_alpha_volume,do_alpha_volume> Visitor;
    applyToSubProblems(Visitor::add_data(wrap_operator_type(spatial_operator()),wrap_eg(eg),
                                         wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(x_s),wrap_r(r_s)));
  }

  template<typename EG, typename LFSV>
  void assembleVVolume(const EG& eg, const LFSV& lfsv)
  {
    typedef visitor<invoke_lambda_volume,do_lambda_volume> Visitor;
    applyToSubProblems(Visitor::add_data(wrap_operator_type(spatial_operator()),wrap_eg(eg),
                                         wrap_lfsv(lfsv),wrap_r(r_s)));
  }


  template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
  void assembleUVSkeleton(const IG& ig,
                          const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                          const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    typedef visitor<invoke_alpha_skeleton_or_boundary,do_alpha_skeleton_or_boundary> SubProblemVisitor;
    applyToSubProblems(SubProblemVisitor::add_data(wrap_operator_type(spatial_operator()),wrap_ig(ig),
                                                   store_neighbor_accessed(false),
                                                   wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(x_s),wrap_r_s(r_s),
                                                   wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(x_n),wrap_r_n(r_n)));

    typedef visitor<invoke_alpha_coupling,do_alpha_coupling> CouplingVisitor;
    applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(coupling_operator()),wrap_ig(ig),
                                               store_neighbor_accessed(false),
                                               wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(x_s),wrap_r_s(r_s),
                                               wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(x_n),wrap_r_n(r_n)));
  }

  template<typename IG, typename LFSV_S, typename LFSV_N>
  void assembleVSkeleton(const IG& ig,
                         const LFSV_S& lfsv_s,
                         const LFSV_N& lfsv_n)
  {
    typedef visitor<invoke_lambda_skeleton_or_boundary,do_lambda_skeleton_or_boundary> Visitor;
    applyToSubProblems(Visitor::add_data(wrap_operator_type(spatial_operator()),wrap_ig(ig),
                                         wrap_lfsv_s(lfsv_s),wrap_r_s(r_s),
                                         wrap_lfsv_n(lfsv_n),wrap_r_n(r_n)));
  }


  template<typename IG, typename LFSU, typename LFSV>
  void assembleUVBoundary(const IG& ig, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<invoke_alpha_boundary,do_alpha_boundary> Visitor;
    applyToSubProblems(Visitor::add_data(wrap_operator_type(spatial_operator()),wrap_ig(ig),
                                         wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(x_s),wrap_r(r_s)));
  }

  template<typename IG, typename LFSV>
  void assembleVBoundary(const IG& ig, const LFSV& lfsv)
  {
    typedef visitor<invoke_lambda_boundary,do_lambda_boundary> Visitor;
    applyToSubProblems(Visitor::add_data(wrap_operator_type(spatial_operator()),wrap_ig(ig),
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
    typedef visitor<invoke_alpha_enriched_coupling,do_alpha_enriched_coupling> CouplingVisitor;
    applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(coupling_operator()),wrap_ig(ig),
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
    typedef visitor<invoke_lambda_enriched_coupling,do_lambda_enriched_coupling> CouplingVisitor;
    applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(coupling_operator()),wrap_ig(ig),
                                               wrap_lfsv_s(lfsv_s),wrap_r_s(r_s),
                                               wrap_lfsv_n(lfsv_n),wrap_r_n(r_n),
                                               wrap_lfsv_c(lfsv_c),wrap_r_c(r_c)));
  }


  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolumePostSkeleton(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<invoke_alpha_volume_post_skeleton,do_alpha_volume_post_skeleton> Visitor;
    applyToSubProblems(Visitor::add_data(wrap_operator_type(spatial_operator()),wrap_eg(eg),
                                         wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(x_s),wrap_r(r_s)));
  }

  template<typename EG, typename LFSV>
  void assembleVVolumePostSkeleton(const EG& eg, const LFSV& lfsv)
  {
    typedef visitor<invoke_lambda_volume_post_skeleton,do_lambda_volume_post_skeleton> Visitor;
    applyToSubProblems(Visitor::add_data(wrap_operator_type(spatial_operator()),wrap_eg(eg),
                                         wrap_lfsv(lfsv),wrap_r(r_s)));
  }


  void preAssembly()
  {}

  void postAssembly()
  {
    Dune::PDELab::constrain_residual(*pconstraintsv,r);
  }

  const X& x;
  R& r;

  LocalVector<typename X::ElementType, TrialSpaceTag> x_s;
  LocalVector<typename R::ElementType, TestSpaceTag> r_s;
  LocalVector<typename X::ElementType, TrialSpaceTag> x_n;
  LocalVector<typename R::ElementType, TestSpaceTag> r_n;
  LocalVector<typename X::ElementType, TrialSpaceTag> x_c;
  LocalVector<typename R::ElementType, TestSpaceTag> r_c;

};



} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_RESIDUALASSEMBLERENGINE_HH
