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

namespace {

  template<typename Domain, typename Jacobian>
  struct JacobianAssemblerEngineData
  {

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

    JacobianAssemblerEngineData()
      : x(nullptr)
      , a(nullptr)
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

  };

} // anonymous namespace

template<typename LA>
class JacobianAssemblerEngine
  : public LocalAssemblerEngineBase
{

  template<typename>
  friend class JacobianAssemblerEngine;

public:

  typedef LA LocalAssembler;
  typedef typename LA::Domain Domain;
  typedef typename LA::Jacobian Jacobian;

  typedef JacobianAssemblerEngineData<Domain,Jacobian> EngineData;

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


  template<typename EG, typename LFSU, typename LFSV>
  void onBindLFSUV(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    if (localAssembler().readData())
      {
        data().x_s.resize(lfsu.size());
        // initialize local jacobian matrix
        data()._a_ss.assign(lfsv.size(),lfsu.size(),0.0);
      }
  }

  template<typename IG, typename LFSU_N, typename LFSV_N>
  void onBindLFSUVOutside(const IG& ig, const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    if (localAssembler().readData())
      {
        data().x_n.resize(lfsu_n.size());
        // initialize local jacobian matrix
        data()._a_nn.assign(lfsv_n.size(),lfsu_n.size(),0.0);
        data()._a_sn.assign(data()._a_ss.nrows(),lfsu_n.size(),0.0);
        data()._a_ns.assign(lfsv_n.size(),data()._a_ss.ncols(),0.0);
      }
  }

  template<typename IG, typename LFSU_C, typename LFSV_C>
  void onBindLFSUVCoupling(const IG& ig, const LFSU_C& lfsu_c, const LFSV_C& lfsv_c)
  {
    if (localAssembler().readData())
      {
        data().x_c.resize(lfsu_c.size());
        // initialize local jacobian matrix
        data()._a_cc.assign(lfsv_c.size(),lfsu_c.size(),0.0);
        data()._a_sc.assign(data()._a_ss.nrows(),lfsu_c.size(),0.0);
        data()._a_cs.assign(lfsv_c.size(),data()._a_ss.ncols(),0.0);
        data()._a_nc.assign(data()._a_nn.nrows(),lfsu_c.size(),0.0);
        data()._a_cn.assign(lfsv_c.size(),data()._a_nn.ncols(),0.0);
      }
  }


  template<typename LFSU>
  void loadCoefficientsLFSUInside(const LFSU& lfsu_s)
  {
    // read local data
    if (localAssembler().readData())
      lfsu_s.vread(*data().x,data().x_s);
  }

  template<typename LFSU_N>
  void loadCoefficientsLFSUOutside(const LFSU_N& lfsu_n)
  {
    // read local data
    if (localAssembler().readData())
      lfsu_n.vread(*data().x,data().x_n);
  }

  template<typename LFSU_C>
  void loadCoefficientsLFSUCoupling(const LFSU_C& lfsu_c)
  {
    // read local data
    if (localAssembler().readData())
      lfsu_c.vread(*data().x,data().x_c);
  }


  template<typename LFSU_S, typename LFSV_S>
  void writeResultsInside(const LFSU_S & lfsu_s, const LFSV_S & lfsv_s)
  {
    // write back local jacobian contributions
    if (localAssembler().writeData())
      {
        if (data().a_ss.modified())
          localAssembler().etadd(lfsv_s,lfsu_s,data()._a_ss,*data().a);
        data().a_ss.resetModified();
      }
  }

  template<typename LFSU_S, typename LFSV_S,
           typename LFSU_N, typename LFSV_N>
  void writeResultsOutside(const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                           const LFSU_N & lfsu_n, const LFSV_N & lfsv_n)
  {
    // write back local jacobian contributions
    if (localAssembler().writeData())
      {
        if (data().a_nn.modified())
          localAssembler().etadd(lfsv_n,lfsu_n,data()._a_nn,*data().a);
        data().a_nn.resetModified();
        if (data().a_sn.modified())
          localAssembler().etadd(lfsv_s,lfsu_n,data()._a_sn,*data().a);
        data().a_sn.resetModified();
        if (data().a_ns.modified())
          localAssembler().etadd(lfsv_n,lfsu_s,data()._a_ns,*data().a);
        data().a_ns.resetModified();
      }
  }

  template<typename LFSU_S, typename LFSV_S,
           typename LFSU_N, typename LFSV_N,
           typename LFSU_C, typename LFSV_C>
  void writeResultsCoupling(const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                            const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                            const LFSU_C & lfsu_c, const LFSV_C & lfsv_c)
  {
    // write back local jacobian contributions
    if (localAssembler().writeData())
      {
        if (data().a_cc.modified())
          localAssembler().etadd(lfsv_c,lfsu_c,data()._a_cc,*data().a);
        data().a_cc.resetModified();
        if (data().a_sc.modified())
          localAssembler().etadd(lfsv_s,lfsu_c,data()._a_sc,*data().a);
        data().a_sc.resetModified();
        if (data().a_cs.modified())
          localAssembler().etadd(lfsv_c,lfsu_s,data()._a_cs,*data().a);
        data().a_cs.resetModified();
        if (data().a_nc.modified())
          localAssembler().etadd(lfsv_n,lfsu_c,data()._a_nc,*data().a);
        data().a_nc.resetModified();
        if (data().a_cn.modified())
          localAssembler().etadd(lfsv_c,lfsu_n,data()._a_cn,*data().a);
        data().a_cn.resetModified();
      }
  }


  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolume(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<functors::invoke_jacobian_volume,do_alpha_volume<> > Visitor;
    data().a_ss.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(DefaultOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(data().x_s),wrap_a_ss(data().a_ss)));
  }

  template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
  void assembleUVSkeleton(const IG& ig,
                          const LFSU_S& lfsu_s, const LFSV_S& lfsv_s,
                          const LFSU_N& lfsu_n, const LFSV_N& lfsv_n)
  {
    data().a_ss.setWeight(localAssembler().weight());
    data().a_sn.setWeight(localAssembler().weight());
    data().a_ns.setWeight(localAssembler().weight());
    data().a_nn.setWeight(localAssembler().weight());
    typedef visitor<functors::invoke_jacobian_skeleton_or_boundary,do_alpha_skeleton_or_boundary<> > SubProblemVisitor;
    localAssembler().applyToSubProblems(SubProblemVisitor::add_data(wrap_operator_type(DefaultOperator()),wrap_ig(ig),
                                                                    wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(data().x_s),
                                                                    wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(data().x_n),
                                                                    wrap_a_ss(data().a_ss),wrap_a_sn(data().a_sn),
                                                                    wrap_a_ns(data().a_ns),wrap_a_nn(data().a_nn)));

    typedef visitor<functors::invoke_jacobian_coupling,do_alpha_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(DefaultOperator()),wrap_ig(ig),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(data().x_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(data().x_n),
                                                                wrap_a_ss(data().a_ss),wrap_a_sn(data().a_sn),
                                                                wrap_a_ns(data().a_ns),wrap_a_nn(data().a_nn)));
  }

  template<typename IG, typename LFSU, typename LFSV>
  void assembleUVBoundary(const IG& ig, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<functors::invoke_jacobian_boundary,do_alpha_boundary<> > Visitor;
    data().a_ss.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(DefaultOperator()),wrap_ig(ig),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(data().x_s),wrap_a(data().a_ss)));
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
    data().a_ss.setWeight(localAssembler().weight());
    data().a_sc.setWeight(localAssembler().weight());
    data().a_nn.setWeight(localAssembler().weight());
    data().a_nc.setWeight(localAssembler().weight());
    data().a_cc.setWeight(localAssembler().weight());
    data().a_cs.setWeight(localAssembler().weight());
    data().a_cn.setWeight(localAssembler().weight());
    typedef visitor<functors::invoke_jacobian_enriched_coupling,do_alpha_enriched_coupling<> > CouplingVisitor;
    localAssembler().applyToCouplings(CouplingVisitor::add_data(wrap_operator_type(DefaultOperator()),wrap_ig(ig),
                                                                wrap_lfsu_s(lfsu_s),wrap_lfsv_s(lfsv_s),wrap_x_s(data().x_s),
                                                                wrap_lfsu_n(lfsu_n),wrap_lfsv_n(lfsv_n),wrap_x_n(data().x_n),
                                                                wrap_lfsu_c(lfsu_c),wrap_lfsv_c(lfsv_c),wrap_x_c(data().x_c),
                                                                wrap_a_ss(data().a_ss),wrap_a_sc(data().a_sc),
                                                                wrap_a_nn(data().a_nn),wrap_a_nc(data().a_nc),
                                                                wrap_a_cc(data().a_cc),wrap_a_cs(data().a_cs),wrap_a_cn(data().a_cn)));
  }

  template<typename EG, typename LFSU, typename LFSV>
  void assembleUVVolumePostSkeleton(const EG& eg, const LFSU& lfsu, const LFSV& lfsv)
  {
    typedef visitor<functors::invoke_jacobian_volume_post_skeleton,do_alpha_volume_post_skeleton<> > Visitor;
    data().a_ss.setWeight(localAssembler().weight());
    localAssembler().applyToSubProblems(Visitor::add_data(wrap_operator_type(DefaultOperator()),wrap_eg(eg),
                                                          wrap_lfsu(lfsu),wrap_lfsv(lfsv),wrap_x(data().x_s),wrap_a(data().a_ss)));
  }


  void postAssembly()
  {
    if (localAssembler().writeData())
      localAssembler().handle_dirichlet_constraints(*data().a);
  }

  void setSolution(const Domain& x_)
  {
    if (localAssembler().readData())
      data().x = &x_;
  }

  void setJacobian(Jacobian& a_)
  {
    if (localAssembler().readData())
      data().a = &a_;
  }

  const LocalAssembler& localAssembler() const
  {
    return _localAssembler;
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
