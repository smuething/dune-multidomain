// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_JACOBIANASSEMBLERENGINE_HH
#define DUNE_PDELAB_MULTIDOMAIN_JACOBIANASSEMBLERENGINE_HH

#include <algorithm>

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>

#include <dune/pdelab/multidomain/datawrappers.hh>
#include <dune/pdelab/multidomain/jacobianassemblerfunctors.hh>


namespace Dune {

namespace PDELab {

namespace MultiDomain {


template<typename LA>
class JacobianAssemblerEngine
{

public:

  typedef LA LocalAssembler;
  typedef typename LA::Domain Domain;
  typedef typename LA::Jacobian Jacobian;


  bool requireIntersections() const
  {
    return requireUVSkeleton() || requireUVEnrichedCoupling() || requireUVBoundary();
  }

  bool requireIntersectionsTwoSided() const
  {
    return requireUVBoundary() ||
      any_child<typename LocalAssembler::SubProblems,do_skeleton_two_sided<> >::value;
  }

  bool requireUVVolume() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_alpha_volume<> >::value;
  }

  bool requireVVolume() const
  {
    return false;
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
    return false;
  }

  bool requireUVBoundary() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_alpha_boundary<> >::value;
  }

  bool requireVBoundary() const
  {
    return false;
  }

  bool requireUVEnrichedCoupling() const
  {
    return any_child<typename LocalAssembler::Couplings,do_alpha_enriched_coupling<> >::value;
  }

  bool requireVEnrichedCoupling() const
  {
    return false;
  }

  bool requireUVVolumePostSkeleton() const
  {
    return any_child<typename LocalAssembler::SubProblems,do_alpha_volume_post_skeleton<> >::value;
  }

  bool requireVVolumePostSkeleton() const
  {
    return false;
  }

  template<typename EG, typename LFSU, typename LFSV>
  void onBindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    x_s.resize(lfsu.size());
    // initialize local jacobian matrix
    a_ss.assign(lfsv.size(),lfsu.size(),0.0);
  }

  template<typename EG, typename LFSU, typename LFSV>
  void onUnbindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    // write back local jacobian contributions
    localAssembler().etadd(lfsv,lfsu,a_ss,*a);
  }

  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onBindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    x_n.resize(lfsu_n.size());
    // initialize local jacobian matrix
    a_nn.assign(lfsv_n.size(),lfsu_n.size(),0.0);
  }

  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onUnbindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    // write back local jacobian contributions
    localAssembler().etadd(lfsv_n,lfsu_n,a_nn,*a);
  }

  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onBindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
    x_c.resize(lfsu_c.size());
    // initialize local jacobian matrix
    a_cc.assign(lfsv_c.size(),lfsu_c.size(),0.0);
  }

  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onUnbindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
    // write back local jacobian contributions
    localAssembler().etadd(lfsv_c,lfsu_c,a_cc,*a);
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
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(x_s),wrap_a_ss(a_ss)));
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
    a_sn.assign(lfsv_s.size(),lfsu_n.size(),0.0);
    a_ns.assign(lfsv_n.size(),lfsu_s.size(),0.0);
    typedef visitor<invoke_jacobian_skeleton_or_boundary,do_alpha_skeleton_or_boundary<> > SubProblemVisitor;
    localAssembler().applyToSubProblems(SubProblemVisitor::add_data(wrap_operator_type(SpatialOperator()),wrap_ig(ig),
                                                                    store_neighbor_accessed(false),
                                                                    wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(x_s),
                                                                    wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(x_n),
                                                                    wrap_a_ss(a_ss),wrap_a_sn(a_sn),
                                                                    wrap_a_ns(a_ns),wrap_a_nn(a_nn)));

    typedef visitor<invoke_jacobian_coupling,do_alpha_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(CouplingOperator()),wrap_ig(ig),
                                                                store_neighbor_accessed(false),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(x_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(x_n),
                                                                wrap_a_ss(a_ss),wrap_a_sn(a_sn),
                                                                wrap_a_ns(a_ns),wrap_a_nn(a_nn)));
    localAssembler().etadd(lfsv_s,lfsu_n,a_sn,*a);
    localAssembler().etadd(lfsv_n,lfsu_s,a_ns,*a);
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
    typedef visitor<invoke_jacobian_boundary,do_alpha_boundary<> > Visitor;
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_ig(ig),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(x_s),wrap_a(a_ss)));
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
    a_sc.assign(lfsv_s.size(),lfsu_c.size(),0.0);
    a_cs.assign(lfsv_c.size(),lfsu_s.size(),0.0);
    a_nc.assign(lfsv_n.size(),lfsu_c.size(),0.0);
    a_cn.assign(lfsv_c.size(),lfsu_n.size(),0.0);
    typedef visitor<invoke_jacobian_enriched_coupling,do_alpha_enriched_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(CouplingOperator()),wrap_ig(ig),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(x_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(x_n),
                                                                wrap_lfsu_c(lfsu_c),wrap_lfsv_c(lfsv_c),wrap_x_c(x_c),
                                                                wrap_a_ss(a_ss),wrap_a_sc(a_sc),
                                                                wrap_a_nn(a_nn),wrap_a_nc(a_nc),
                                                                wrap_a_cc(a_cc),wrap_a_cs(a_cs),wrap_a_cn(a_cn)));
    localAssembler().etadd(lfsv_s,lfsu_c,a_sc,*a);
    localAssembler().etadd(lfsv_c,lfsu_s,a_cs,*a);
    localAssembler().etadd(lfsv_n,lfsu_c,a_nc,*a);
    localAssembler().etadd(lfsv_c,lfsu_n,a_cn,*a);
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
    typedef visitor<invoke_jacobian_volume_post_skeleton,do_alpha_volume_post_skeleton<> > Visitor;
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(SpatialOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(x_s),wrap_a(a_ss)));
  }

  template<typename EG, typename LFSV>
  void assembleVVolumePostSkeleton(const EG& eg, const LFSV& lfsv)
  {
  }


  void preAssembly()
  {}

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
  {}

private:

  const LocalAssembler& _localAssembler;

  const Domain* x;

  LocalVector<typename Domain::ElementType, TrialSpaceTag> x_s;
  LocalVector<typename Domain::ElementType, TrialSpaceTag> x_n;
  LocalVector<typename Domain::ElementType, TrialSpaceTag> x_c;

  Jacobian* a;

  LocalMatrix<typename Jacobian::ElementType> a_ss;
  LocalMatrix<typename Jacobian::ElementType> a_sn;
  LocalMatrix<typename Jacobian::ElementType> a_sc;

  LocalMatrix<typename Jacobian::ElementType> a_nn;
  LocalMatrix<typename Jacobian::ElementType> a_ns;
  LocalMatrix<typename Jacobian::ElementType> a_nc;

  LocalMatrix<typename Jacobian::ElementType> a_cc;
  LocalMatrix<typename Jacobian::ElementType> a_cs;
  LocalMatrix<typename Jacobian::ElementType> a_cn;

};



} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_JACOBIANASSEMBLERENGINE_HH
