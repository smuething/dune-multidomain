// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_JACOBIANASSEMBLERENGINE_HH
#define DUNE_PDELAB_MULTIDOMAIN_JACOBIANASSEMBLERENGINE_HH

#include <algorithm>

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>
#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>

#include <dune/pdelab/multidomain/datawrappers.hh>
#include <dune/pdelab/multidomain/jacobianassemblerfunctors.hh>


namespace Dune {

namespace PDELab {

namespace MultiDomain {


template<typename LA>
class JacobianAssemblerEngine
  : public LocalAssemblerEngineBase
{

public:

  typedef LA LocalAssembler;
  typedef typename LA::Domain Domain;
  typedef typename LA::Jacobian Jacobian;


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
      any_child<typename LocalAssembler::SubProblems,do_skeleton_two_sided<> >::value;
  }

  bool requireUVVolume() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_alpha_volume<> >::value;
  }

  bool requireUVSkeleton() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_alpha_skeleton<> >::value ||
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

  template<typename EG, typename LFSU, typename LFSV>
  void onBindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    x_s.resize(lfsu.size());
    // initialize local jacobian matrix
    _a_ss.assign(lfsv.size(),lfsu.size(),0.0);
  }

  template<typename EG, typename LFSU, typename LFSV>
  void onUnbindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    // write back local jacobian contributions
    if (a_ss.modified())
      localAssembler().etadd(lfsv,lfsu,_a_ss,*a);
    a_ss.resetModified();
  }

  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onBindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    x_n.resize(lfsu_n.size());
    // initialize local jacobian matrix
    _a_nn.assign(lfsv_n.size(),lfsu_n.size(),0.0);
  }

  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onUnbindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    // write back local jacobian contributions
    if (a_nn.modified())
      localAssembler().etadd(lfsv_n,lfsu_n,_a_nn,*a);
    a_nn.resetModified();
  }

  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onBindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
    x_c.resize(lfsu_c.size());
    // initialize local jacobian matrix
    _a_cc.assign(lfsv_c.size(),lfsu_c.size(),0.0);
  }

  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onUnbindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
    // write back local jacobian contributions
    if (a_cc.modified())
      localAssembler().etadd(lfsv_c,lfsu_c,_a_cc,*a);
    a_cc.resetModified();
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
    typedef visitor<invoke_jacobian_volume,do_alpha_volume<> > Visitor;
    a_ss.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(x_s),wrap_a_ss(a_ss)));
  }

  template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
  void assembleUVSkeleton(const IG& ig,
                          const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                          const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    _a_sn.assign(lfsv_s.size(),lfsu_n.size(),0.0);
    _a_ns.assign(lfsv_n.size(),lfsu_s.size(),0.0);
    a_ss.setWeight(localAssembler().weight());
    a_sn.setWeight(localAssembler().weight());
    a_ns.setWeight(localAssembler().weight());
    a_nn.setWeight(localAssembler().weight());
    typedef visitor<invoke_jacobian_skeleton_or_boundary,do_alpha_skeleton_or_boundary<> > SubProblemVisitor;
    localAssembler().applyToSubProblems(SubProblemVisitor::add_data(wrap_operator_type(SpatialOperator()),wrap_ig(ig),
                                                                    wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(x_s),
                                                                    wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(x_n),
                                                                    wrap_a_ss(a_ss),wrap_a_sn(a_sn),
                                                                    wrap_a_ns(a_ns),wrap_a_nn(a_nn)));

    typedef visitor<invoke_jacobian_coupling,do_alpha_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(CouplingOperator()),wrap_ig(ig),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(x_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(x_n),
                                                                wrap_a_ss(a_ss),wrap_a_sn(a_sn),
                                                                wrap_a_ns(a_ns),wrap_a_nn(a_nn)));
    if (a_sn.modified())
      localAssembler().etadd(lfsv_s,lfsu_n,_a_sn,*a);
    a_sn.resetModified();
    if (a_ns.modified())
      localAssembler().etadd(lfsv_n,lfsu_s,_a_ns,*a);
    a_ns.resetModified();
  }

  template<typename IG, typename LFSU, typename LFSV>
  void assembleUVBoundary(const IG& ig, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<invoke_jacobian_boundary,do_alpha_boundary<> > Visitor;
    a_ss.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_ig(ig),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(x_s),wrap_a(a_ss)));
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
    _a_sc.assign(lfsv_s.size(),lfsu_c.size(),0.0);
    _a_cs.assign(lfsv_c.size(),lfsu_s.size(),0.0);
    _a_nc.assign(lfsv_n.size(),lfsu_c.size(),0.0);
    _a_cn.assign(lfsv_c.size(),lfsu_n.size(),0.0);
    a_ss.setWeight(localAssembler().weight());
    a_sc.setWeight(localAssembler().weight());
    a_nn.setWeight(localAssembler().weight());
    a_nc.setWeight(localAssembler().weight());
    a_cc.setWeight(localAssembler().weight());
    a_cs.setWeight(localAssembler().weight());
    a_cn.setWeight(localAssembler().weight());
    typedef visitor<invoke_jacobian_enriched_coupling,do_alpha_enriched_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(CouplingOperator()),wrap_ig(ig),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(x_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(x_n),
                                                                wrap_lfsu_c(lfsu_c),wrap_lfsv_c(lfsv_c),wrap_x_c(x_c),
                                                                wrap_a_ss(a_ss),wrap_a_sc(a_sc),
                                                                wrap_a_nn(a_nn),wrap_a_nc(a_nc),
                                                                wrap_a_cc(a_cc),wrap_a_cs(a_cs),wrap_a_cn(a_cn)));
    if (a_sc.modified())
      localAssembler().etadd(lfsv_s,lfsu_c,_a_sc,*a);
    a_sc.resetModified();
    if (a_cs.modified())
      localAssembler().etadd(lfsv_c,lfsu_s,_a_cs,*a);
    a_cs.resetModified();
    if (a_nc.modified())
      localAssembler().etadd(lfsv_n,lfsu_c,_a_nc,*a);
    a_nc.resetModified();
    if (a_cn.modified())
      localAssembler().etadd(lfsv_c,lfsu_n,_a_cn,*a);
    a_cn.resetModified();
  }

  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolumePostSkeleton(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<invoke_jacobian_volume_post_skeleton,do_alpha_volume_post_skeleton<> > Visitor;
    a_ss.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(x_s),wrap_a(a_ss)));
  }


  void postAssembly()
  {
    localAssembler().handle_dirichlet_constraints(*a);
  }

  void setSolution(const Domain& x_)
  {
    x = &x_;
  }

  void setJacobian(Jacobian& a_)
  {
    a = &a_;
  }

  const LocalAssembler& localAssembler() const
  {
    return _localAssembler;
  }

  JacobianAssemblerEngine(const LocalAssembler& localAssembler)
    : _localAssembler(localAssembler)
    , a_ss(_a_ss,1.0)
    , a_sn(_a_sn,1.0)
    , a_sc(_a_sc,1.0)
    , a_nn(_a_nn,1.0)
    , a_ns(_a_ns,1.0)
    , a_nc(_a_nc,1.0)
    , a_cc(_a_cc,1.0)
    , a_cs(_a_cs,1.0)
    , a_cn(_a_cn,1.0)
  {}

private:

  const LocalAssembler& _localAssembler;

  const Domain* x;

  typedef LocalVector<typename Domain::ElementType, TrialSpaceTag> SolutionVector;
  typedef LocalMatrix<typename Jacobian::ElementType> LocalJacobian;

  SolutionVector x_s;
  SolutionVector x_n;
  SolutionVector x_c;

  Jacobian* a;

  LocalJacobian _a_ss;
  LocalJacobian _a_sn;
  LocalJacobian _a_sc;

  LocalJacobian _a_nn;
  LocalJacobian _a_ns;
  LocalJacobian _a_nc;

  LocalJacobian _a_cc;
  LocalJacobian _a_cs;
  LocalJacobian _a_cn;

  typename LocalJacobian::WeightedAccumulationView a_ss;
  typename LocalJacobian::WeightedAccumulationView a_sn;
  typename LocalJacobian::WeightedAccumulationView a_sc;

  typename LocalJacobian::WeightedAccumulationView a_nn;
  typename LocalJacobian::WeightedAccumulationView a_ns;
  typename LocalJacobian::WeightedAccumulationView a_nc;

  typename LocalJacobian::WeightedAccumulationView a_cc;
  typename LocalJacobian::WeightedAccumulationView a_cs;
  typename LocalJacobian::WeightedAccumulationView a_cn;

};



} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_JACOBIANASSEMBLERENGINE_HH
