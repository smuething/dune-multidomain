#ifndef DUNE_MULTIDOMAIN_MULTIDOMAINGRIDOPERATORSPACEINVOCATIONHELPERS_HH
#define DUNE_MULTIDOMAIN_MULTIDOMAINGRIDOPERATORSPACEINVOCATIONHELPERS_HH

namespace Dune {
namespace PDELab {
namespace MultiDomain {

template<typename TReal>
struct SetTime
{

  SetTime(TReal t) :
    time(t)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, Child& child)
  {
    child.setTime(time);
  }

private:
  const TReal time;
};


struct PreStep
{

  template<typename Data, typename Child>
  void operator()(Data& data, Child& child)
  {
    const typename Data::GOS& gos = data.gos();
    child.preStep(gos.time,gos.dt,gos.method->s());
  }

};


struct PostStep
{

  template<typename Data, typename Child>
  void operator()(Data& data, Child& child)
  {
    child.postStep();
  }

};


template<typename TReal,typename StageType>
struct PreStage
{

  PreStage(TReal t, StageType s) :
    time(t),
    stage(s)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, Child& child)
  {
    child.preStage(time,stage);
  }

private:
  const TReal time;
  const StageType stage;
};


struct PostStage
{

  template<typename Data, typename Child>
  void operator()(Data& data, Child& child)
  {
    child.postStage();
  }

};


template<typename TReal>
struct SuggestTimestep
{

  SuggestTimestep(TReal dt) :
    _dt(dt),
    _suggested_dt(std::numeric_limits<TReal>::max())
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, Child& child)
  {
    _suggested_dt = std::min(_suggested_dt,child.suggestTimestep(_dt));
  }

  TReal value() const {
    return _suggested_dt;
  }

private:
  const TReal _dt;
  TReal _suggested_dt;

};


template<typename P, typename Operator=SpatialOperator>
struct BuildVolumePattern
{
  BuildVolumePattern(P& gp) : globalpattern(gp) {}

  template<typename Data, typename Child>
  void operator()(Data& data, Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains()))
      return;
    typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
    typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
    LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsu.bind();
    lfsv.bind();
    LocalSparsityPattern localpattern;
    Operator::extract(child).pattern_volume(lfsu,lfsv,localpattern);

    // translate local to global indices and add to global pattern
    for (size_t k=0; k<localpattern.size(); ++k)
      // TODO: Figure out if we really need the localIndex() call
      data.gos().add_entry(globalpattern,
                           data.lfsu().globalIndex(lfsu.localIndex(localpattern[k].i())),
                           data.lfsv().globalIndex(lfsv.localIndex(localpattern[k].j()))
                           );
  }

  P& globalpattern;
};

template<typename P, typename Operator=SpatialOperator>
struct BuildSkeletonPattern
{
  BuildSkeletonPattern(P& gp) : globalpattern(gp) {}

  template<typename Data, typename Child>
  void operator()(Data& data, Child& child)
  {
    // skip subdomain borders that do not coincide with the overall domain border
    if (!(child.appliesTo(data.elementSubDomains()) && child.appliesTo(data.neighborSubDomains())))
      return;
    typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
    typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
    LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsu.bind();
    lfsv.bind();
    LFSU lfsun(data.lfsun(),child,child.trialGridFunctionSpaceConstraints());
    LFSV lfsvn(data.lfsvn(),child,child.testGridFunctionSpaceConstraints());
    lfsun.bind();
    lfsvn.bind();
    LocalSparsityPattern localpattern_sn, localpattern_ns;
    Operator::extract(child).pattern_skeleton(lfsu,lfsv,lfsun,lfsvn,localpattern_sn,localpattern_ns);

    // translate local to global indices and add to global pattern
    for (size_t k=0; k<localpattern_sn.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.lfsu().globalIndex(lfsu.localIndex(localpattern_sn[k].i())),
                           data.lfsvn().globalIndex(lfsvn.localIndex(localpattern_sn[k].j()))
                           );

    for (size_t k=0; k<localpattern_ns.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.lfsun().globalIndex(lfsun.localIndex(localpattern_ns[k].i())),
                           data.lfsv().globalIndex(lfsv.localIndex(localpattern_ns[k].j()))
                           );
  }

  P& globalpattern;
};


template<typename P, typename Operator = CouplingOperator>
struct BuildCouplingPattern
{
  BuildCouplingPattern(P& gp) : globalpattern(gp) {}

  template<typename Data, typename Child>
  void operator()(Data& data, Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains(),data.neighborSubDomains()))
      return;
    typedef typename Operator::template ExtractType<Child>::Type LOP;
    typedef typename Child::Traits::LocalSubProblem LocalSubProblem;
    typedef typename Child::Traits::RemoteSubProblem RemoteSubProblem;
    typedef typename LocalSubProblem::Traits::TrialLocalFunctionSpace LocalLFSU;
    typedef typename LocalSubProblem::Traits::TestLocalFunctionSpace LocalLFSV;
    typedef typename RemoteSubProblem::Traits::TrialLocalFunctionSpace RemoteLFSU;
    typedef typename RemoteSubProblem::Traits::TestLocalFunctionSpace RemoteLFSV;
    const LocalSubProblem& localSubProblem = child.localSubProblem();
    const RemoteSubProblem& remoteSubProblem = child.remoteSubProblem();
    LocalLFSU local_lfsu(data.lfsu(),localSubProblem,localSubProblem.trialGridFunctionSpaceConstraints());
    LocalLFSV local_lfsv(data.lfsv(),localSubProblem,localSubProblem.testGridFunctionSpaceConstraints());
    local_lfsu.bind();
    local_lfsv.bind();
    RemoteLFSU remote_lfsu(data.lfsun(),remoteSubProblem,remoteSubProblem.trialGridFunctionSpaceConstraints());
    RemoteLFSV remote_lfsv(data.lfsvn(),remoteSubProblem,remoteSubProblem.testGridFunctionSpaceConstraints());
    remote_lfsu.bind();
    remote_lfsv.bind();
    LocalSparsityPattern localpattern_sn, localpattern_ns;
    Operator::extract(child).pattern_coupling(local_lfsu,local_lfsv,remote_lfsu,remote_lfsv,localpattern_sn,localpattern_ns);

    // translate local to global indices and add to global pattern
    // FIXME: this only works because of some kind of miracle!
    for (size_t k=0; k<localpattern_sn.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.lfsu().globalIndex(local_lfsu.localIndex(localpattern_sn[k].i())),
                           data.lfsvn().globalIndex(remote_lfsv.localIndex(localpattern_sn[k].j()))
                           );
    for (size_t k=0; k<localpattern_ns.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.lfsun().globalIndex(remote_lfsu.localIndex(localpattern_ns[k].i())),
                           data.lfsv().globalIndex(local_lfsv.localIndex(localpattern_ns[k].j()))
                           );
  }

  P& globalpattern;
};


template<typename P, typename Operator = CouplingOperator>
struct BuildEnrichedCouplingPattern
{
  BuildEnrichedCouplingPattern(P& gp) : globalpattern(gp) {}

  template<typename Data, typename Child>
  void operator()(Data& data, Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains(),data.neighborSubDomains()))
      return;
    typedef typename Operator::template ExtractType<Child>::Type LOP;
    typedef typename Child::Traits::LocalSubProblem LocalSubProblem;
    typedef typename Child::Traits::RemoteSubProblem RemoteSubProblem;
    typedef typename LocalSubProblem::Traits::TrialLocalFunctionSpace LocalLFSU;
    typedef typename LocalSubProblem::Traits::TestLocalFunctionSpace LocalLFSV;
    typedef typename RemoteSubProblem::Traits::TrialLocalFunctionSpace RemoteLFSU;
    typedef typename RemoteSubProblem::Traits::TestLocalFunctionSpace RemoteLFSV;
    const LocalSubProblem& localSubProblem = child.localSubProblem();
    const RemoteSubProblem& remoteSubProblem = child.remoteSubProblem();
    LocalLFSU local_lfsu(data.lfsu(),localSubProblem,localSubProblem.trialGridFunctionSpaceConstraints());
    LocalLFSV local_lfsv(data.lfsv(),localSubProblem,localSubProblem.testGridFunctionSpaceConstraints());
    local_lfsu.bind();
    local_lfsv.bind();
    RemoteLFSU remote_lfsu(data.lfsun(),remoteSubProblem,remoteSubProblem.trialGridFunctionSpaceConstraints());
    RemoteLFSV remote_lfsv(data.lfsvn(),remoteSubProblem,remoteSubProblem.testGridFunctionSpaceConstraints());
    remote_lfsu.bind();
    remote_lfsv.bind();
    typedef typename Child::template CouplingLocalFunctionSpace<typename Data::LFSU>::Type CouplingLFSU;
    typedef typename Child::template CouplingLocalFunctionSpace<typename Data::LFSV>::Type CouplingLFSV;
    const CouplingLFSU& coupling_lfsu = child.couplingLocalFunctionSpace(data.couplinglfsu());
    const CouplingLFSV& coupling_lfsv = child.couplingLocalFunctionSpace(data.couplinglfsv());
    LocalSparsityPattern localpattern_sn, localpattern_ns, localpattern_sc, localpattern_cs, localpattern_nc, localpattern_cn, localpattern_coupling;

    // handle pattern internal to the coupling
    Operator::extract(child).pattern_enriched_coupling(coupling_lfsu,
                                                       coupling_lfsv,
                                                       localpattern_coupling);

    // handle coupling between first subproblem and enrichment space
    Operator::extract(child).pattern_enriched_coupling_first(local_lfsu,
                                                             local_lfsv,
                                                             coupling_lfsu,
                                                             coupling_lfsv,
                                                             localpattern_sc,
                                                             localpattern_cs);
    // handle coupling between second subproblem and enrichment space
    Operator::extract(child).pattern_enriched_coupling_second(remote_lfsu,
                                                              remote_lfsv,
                                                              coupling_lfsu,
                                                              coupling_lfsv,
                                                              localpattern_nc,
                                                              localpattern_cn);
    // translate local to global indices and add to global pattern
    // FIXME: this only works because of some kind of miracle!
    for (size_t k=0; k<localpattern_sn.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.lfsu().globalIndex(local_lfsu.localIndex(localpattern_sn[k].i())),
                           data.lfsvn().globalIndex(remote_lfsv.localIndex(localpattern_sn[k].j()))
                           );
    for (size_t k=0; k<localpattern_ns.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.lfsun().globalIndex(remote_lfsu.localIndex(localpattern_ns[k].i())),
                           data.lfsv().globalIndex(local_lfsv.localIndex(localpattern_ns[k].j()))
                           );
    for (size_t k=0; k<localpattern_sc.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.lfsu().globalIndex(local_lfsu.localIndex(localpattern_sc[k].i())),
                           data.couplinglfsv().globalIndex(coupling_lfsv.localIndex(localpattern_sc[k].j()))
                           );
    for (size_t k=0; k<localpattern_cs.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.couplinglfsu().globalIndex(coupling_lfsu.localIndex(localpattern_cs[k].i())),
                           data.lfsv().globalIndex(local_lfsv.localIndex(localpattern_cs[k].j()))
                           );
    for (size_t k=0; k<localpattern_nc.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.lfsun().globalIndex(remote_lfsu.localIndex(localpattern_nc[k].i())),
                           data.couplinglfsv().globalIndex(coupling_lfsv.localIndex(localpattern_nc[k].j()))
                           );
    for (size_t k=0; k<localpattern_cn.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.couplinglfsu().globalIndex(coupling_lfsu.localIndex(localpattern_cn[k].i())),
                           data.lfsvn().globalIndex(remote_lfsv.localIndex(localpattern_cn[k].j()))
                           );
    for (size_t k=0; k<localpattern_coupling.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.couplinglfsu().globalIndex(coupling_lfsu.localIndex(localpattern_coupling[k].i())),
                           data.couplinglfsv().globalIndex(coupling_lfsv.localIndex(localpattern_coupling[k].j()))
                           );
  }

  P& globalpattern;
};

namespace LFS {

template<typename Data, typename SubProblem>
struct LFSU
{
  typedef typename SubProblem::Traits::TrialLocalFunctionSpace type;
};

template<typename Data, typename SubProblem>
inline typename LFSU<Data,SubProblem>::type
lfsu(const Data& data, const SubProblem& subProblem)
{
  typedef typename LFSU<Data,SubProblem>::type LFS;
  LFSU lfs(data.lfsu(),subProblem,subProblem.trialGridFunctionSpaceConstraints());
  lfs.bind();
  return lfs;
}

template<typename Data, typename SubProblem>
struct LFSV
{
  typedef typename SubProblem::Traits::TestLocalFunctionSpace type;
};

template<typename Data, typename SubProblem>
inline typename LFSV<Data,SubProblem>::type
lfsv(const Data& data, const SubProblem& subProblem)
{
  typedef typename LFSV<Data,SubProblem>::type LFS;
  LFS lfs(data.lfsv(),subProblem,subProblem.testGridFunctionSpaceConstraints());
  lfs.bind();
  return lfs;
}

template<typename Data, typename SubProblem>
struct LFSU_S
{
  typedef typename SubProblem::Traits::TrialLocalFunctionSpace type;
};

template<typename Data, typename SubProblem>
inline typename LFSU_S<Data,SubProblem>::type
lfsu_s(const Data& data, const SubProblem& subProblem)
{
  typedef typename LFSU_S<Data,SubProblem>::type LFS;
  LFS lfs(data.lfsu_s(),subProblem,subProblem.trialGridFunctionSpaceConstraints());
  lfs.bind();
  return lfs;
}

template<typename Data, typename SubProblem>
struct LFSV_S
{
  typedef typename SubProblem::Traits::TestLocalFunctionSpace type;
};

template<typename Data, typename SubProblem>
inline typename LFSV_S<Data,SubProblem>::type
lfsv_s(const Data& data, const SubProblem& subProblem)
{
  typedef typename LFSV_S<Data,SubProblem>::type LFS;
  LFS lfs(data.lfsv_s(),subProblem,subProblem.testGridFunctionSpaceConstraints());
  lfs.bind();
  return lfs;
}

template<typename Data, typename SubProblem>
struct LFSU_N
{
  typedef typename SubProblem::Traits::TrialLocalFunctionSpace type;
};

template<typename Data, typename SubProblem>
inline typename LFSU_N<Data,SubProblem>::type
lfsu_n(const Data& data, const SubProblem& subProblem)
{
  typedef typename LFSU_N<Data,SubProblem>::type LFS;
  LFS lfs(data.lfsu_n(),subProblem,subProblem.trialGridFunctionSpaceConstraints());
  lfs.bind();
  return lfs;
}

template<typename Data, typename SubProblem>
struct LFSV_N
{
  typedef typename SubProblem::Traits::TestLocalFunctionSpace type;
};

template<typename Data, typename SubProblem>
inline typename LFSV_N<Data,SubProblem>::type
lfsv_n(const Data& data, const SubProblem& subProblem)
{
  typedef typename LFSV_N<Data,SubProblem>::type LFS;
  LFS lfs(data.lfsv_n(),subProblem,subProblem.testGridFunctionSpaceConstraints());
  lfs.bind();
  return lfs;
}

template<typename Data, typename Coupling>
struct LFSU_C
{
  typedef typename Coupling::template CouplingLocalFunctionSpace<typename Data::LFSU_C>::Type type;
};

template<typename Data, typename Coupling>
inline typename LFSU_C<Data,Coupling>::type
lfsu_c(const Data& data, const Coupling& coupling)
{
  return coupling.couplingLocalFunctionSpace(data.lfsu_c());
}

template<typename Data, typename Coupling>
struct LFSV_C
{
  typedef typename Coupling::template CouplingLocalFunctionSpace<typename Data::LFSV_C>::Type type;
};

template<typename Data, typename Coupling>
inline typename LFSV_C<Data,Coupling>::type
lfsv_c(const Data& data, const Coupling& coupling)
{
  return coupling.couplingLocalFunctionSpace(data.lfsv_c());
}

} // namespace LFS


template<typename data_container>
struct invoke_alpha_volume
  : public data_accessor<data_container>
{

  typedef typename data_accessor<data_container> Data;
  using data_accessor<data_container>::data;

  template<typename SubProblem>
  void operator()(const SubProblem& subProblem)
  {
    if (!subProblem.appliesTo(data().eg()))
      return;
    Data::Operator::extract(subProblem).alpha_volume(data().eg(),
                                                     LFS::lfsu(data(),subProblem),
                                                     data().xl(),
                                                     LFS::lfsv(data(),subProblem),
                                                     data().rl());
  }

};


template<typename data_container>
struct invoke_lambda_volume
  : public data_accessor<data_container>
{

  typedef typename data_accessor<data_container> Data;
  using data_accessor<data_container>::data;

  template<typename SubProblem>
  void operator()(const SubProblem& subProblem)
  {
    if (!subProblem.appliesTo(data().eg()))
      return;
    Data::Operator::extract(subProblem).lambda_volume(data().eg(),
                                                      LFS::LFSV(data(),subProblem),
                                                      data().rl());
  }

};


template<typename data_container>
struct invoke_alpha_skeleton_or_boundary
  : public data_accessor<data_container>
{

  typedef typename data_accessor<data_container> Data;
  using data_accessor<data_container>::data;

  template<typename SubProblem>
  void operator()(const SubProblem& subProblem)
  {
    if (!subProblem.appliesTo(data().ig().innerElement()))
      return;
    typedef typename Data::Operator::template ExtractType<SubProblem>::Type LOP;
    typename LFS::LFSU_S<Data,SubProblem>::type lfsu_s = LFS::lfsu_s(data(),subProblem);
    typename LFS::LFSV_S<Data,SubProblem>::type lfsv_s = LFS::lfsv_s(data(),subProblem);

    if (subProblem.appliesTo(data().ig().outerElement()))
      {
        if (_applyOneSided || LOP::doSkeletonTwoSided)
          {
            LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
              alpha_skeleton(Data::Operator::extract(subProblem),data().ig(),
                             lfsu_s,data().xl(),lfsv_s,
                             LFS::lfsu_n(data(),subProblem),data().xn(),LFS::lfsv_n(data(),subProblem),
                             data().rl(),data().rn());
            if(LOP::doAlphaSkeleton)
              data().alpha_skeleton_invoked() = true;
          }
      }
    else
      {
        LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          alpha_boundary(Data::Operator::extract(subProblem),
                         data().ig(),
                         lfsu_s,data().xl(),lfsv_s,data().rl());
      }
  }
};


template<typename data_container>
struct invoke_alpha_boundary
  : public data_accessor<data_container>
{

  typedef typename data_accessor<data_container> Data;
  using data_accessor<data_container>::data;

  template<typename SubProblem>
  void operator()(const SubProblem& subProblem)
  {
    if (!subProblem.appliesTo(data().eg()))
      return;
    Data::Operator::extract(subProblem).alpha_boundary(data().ig(),
                                                       LFS::lfsu(data(),subProblem),
                                                       data().xl(),
                                                       LFS::lfsv(data(),subProblem),
                                                       data().rl());
  }

};


template<typename data_container>
struct invoke_alpha_coupling
  : public data_accessor<data_container>
{

  typedef typename data_accessor<data_container> Data;
  using data_accessor<data_container>::data;

  template<typename Coupling>
  void operator()(const Coupling& coupling)
  {
    if (!coupling.appliesTo(data().ig()))
      return;
    typedef typename Data::Operator::template ExtractType<Coupling>::Type LOP;
    typedef typename Coupling::Traits::LocalSubProblem LocalSubProblem;
    typedef typename Coupling::Traits::RemoteSubProblem RemoteSubProblem;
    const LocalSubProblem& localSubProblem = coupling.localSubProblem();
    const RemoteSubProblem& remoteSubProblem = coupling.remoteSubProblem();
    Data::Operator::extract(coupling).alpha_coupling(data().ig(),
                                                     LFS::lfsu_s(data(),localSubProblem),
                                                     data().xl(),
                                                     LFS::lfsv_s(data(),localSubProblem),
                                                     LFS::lfsu_n(data(),RemoteSubProblem),
                                                     data().xl(),
                                                     LFS::lfsv_n(data(),remoteSubProblem),
                                                     data().rl(),data().rn());
    data().alpha_skeleton_invoked() = true;
  }

};


template<typename data_container>
struct invoke_alpha_enriched_coupling
  : public data_accessor<data_container>
{

  typedef typename data_accessor<data_container> Data;
  using data_accessor<data_container>::data;

  template<typename Coupling>
  void operator()(const Coupling& coupling)
  {
    if (!coupling.appliesTo(data().ig()))
      return;
    typedef typename Data::Operator::template ExtractType<Coupling>::Type LOP;
    typedef typename Coupling::Traits::LocalSubProblem LocalSubProblem;
    typedef typename Coupling::Traits::RemoteSubProblem RemoteSubProblem;
    typename
    Data::Operator::extract(coupling).alpha_enriched_coupling_first(data().ig(),
                                                                    LFS::lfsu_s(data(),localSubProblem),
                                                                    data().xl(),
                                                                    LFS::lfsv_s(data(),localSubProblem),
                                                                    lfsu_c,data().xc(),lfsv_c,
                                                                    data().rl(),data().rc());
    Data::Operator::extract(coupling).alpha_enriched_coupling_second(data().ig(),
                                                                     LFS::lfsu_n(data(),remoteSubProblem),
                                                                     data().xl(),
                                                                     LFS::lfsv_n(data(),remoteSubProblem),
                                                                     lfsu_c,data().xc(),lfsv_c,
                                                                     data().rn(),data().rc());
    data().alphaSkeletonInvoked() = true;
    data().alphaEnrichedCouplingInvoked() = true;
  }

};


template<typename data_container>
struct invoke_lambda_boundary
  : public data_accessor<data_container>
{

  typedef typename data_accessor<data_container> Data;
  using data_accessor<data_container>::data;

  template<typename SubProblem>
  void operator()(const SubProblem& subProblem)
  {
    if (!subProblem.appliesTo(data().ig().inside()))
      return;
    typedef typename SubProblem::Traits::TestLocalFunctionSpace LFSV;
    LFSV lfsv(data().lfsv(),subProblem,subProblem.testGridFunctionSpaceConstraints());
    lfsv.bind();
    Data::Operator::extract(subProblem).lambda_boundary(data().ig(),lfsv,data().rl());
  }

};


template<typename data_container>
struct invoke_alpha_volume_post_skeleton
  : public data_accessor<data_container>
{

  typedef typename data_accessor<data_container> Data;
  using data_accessor<data_container>::data;

  template<typename SubProblem>
  void operator()(const SubProblem& subProblem)
  {
    if (!subProblem.appliesTo(data().eg()))
      return;
    typedef typename SubProblem::Traits::TrialLocalFunctionSpace LFSU;
    typedef typename SubProblem::Traits::TestLocalFunctionSpace LFSV;
    LFSU lfsu(data().lfsu(),subProblem,subProblem.trialGridFunctionSpaceConstraints());
    LFSV lfsv(data().lfsv(),subProblem,subProblem.testGridFunctionSpaceConstraints());
    lfsu.bind();
    lfsv.bind();
    Operator::extract(subProblem).alpha_volume_post_skeleton(data().eg(),lfsu,data().xl(),lfsv,data().rl());
  }

};


template<typename data_container>
struct invoke_lambda_volume_post_skeleton
  : public data_accessor<data_container>
{

  typedef typename data_accessor<data_container> Data;
  using data_accessor<data_container>::data;

  template<typename SubProblem>
  void operator()(const SubProblem& subProblem)
  {
    if (!subProblem.appliesTo(data().eg()))
      return;
    typedef typename SubProblem::Traits::TestLocalFunctionSpace LFSV;
    LFSV lfsv(data().lfsv(),subProblem,subProblem.testGridFunctionSpaceConstraints());
    lfsv.bind();
    Data::Operator::extract(subProblem).lambda_volume_post_skeleton(data().eg(),lfsv,data().rl());
  }

};


template<typename XL, typename YL, typename Operator=SpatialOperator>
struct InvokeJacobianApplyVolume
{

  InvokeJacobianApplyVolume(const XL& xl, YL& yl) :
    x(xl),
    y(yl)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains()))
      return;
    typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
    typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
    LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsu.bind();
    lfsv.bind();
    Operator::extract(child).jacobian_apply_volume(ElementGeometry<typename Data::Element>(data.element()),lfsu,x,lfsv,y);
  }

  const XL& x;
  YL& y;
};


template<typename XL, typename YL, typename Operator=SpatialOperator>
struct InvokeJacobianApplySkeletonOrBoundary
{

  InvokeJacobianApplySkeletonOrBoundary(const XL& xl, const XL& xn, YL& yl, YL& yn, bool applyOneSided) :
    _xl(xl),
    _xn(xn),
    _yl(yl),
    _yn(yn),
    _applyOneSided(applyOneSided)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains()))
      return;
    typedef typename Operator::template ExtractType<Child>::Type LOP;
    typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
    typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
    LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsu.bind();
    lfsv.bind();
    if (child.appliesTo(data.neighborSubDomains()))
      {
        if (_applyOneSided || LOP::doSkeletonTwoSided)
          {
            LFSU lfsun(data.lfsun(),child,child.trialGridFunctionSpaceConstraints());
            LFSV lfsvn(data.lfsvn(),child,child.testGridFunctionSpaceConstraints());
            lfsun.bind();
            lfsvn.bind();
            LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
              jacobian_apply_skeleton(Operator::extract(child),
                                      IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                        data.intersection_index),
                                      lfsu,_xl,lfsv,lfsun,_xn,lfsvn,_yl,_yn);
            if(LOP::doAlphaSkeleton)
              data.setAlphaSkeletonInvoked();
          }
      }
    else
      {
        LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          jacobian_apply_boundary(Operator::extract(child),
                                  IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                    data.intersection_index),
                                  lfsu,_xl,lfsv,_yl);
      }
  }

  const XL& _xl;
  const XL& _xn;
  YL& _yl;
  YL& _yn;
  const bool _applyOneSided;
};


template<typename XL, typename YL, typename Operator=SpatialOperator>
struct InvokeJacobianApplyBoundary
{

  InvokeJacobianApplyBoundary(const XL& xl, YL& yl) :
    x(xl),
    y(yl)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains()))
      return;
    typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
    typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
    LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsu.bind();
    lfsv.bind();
    Operator::extract(child).jacobian_apply_boundary(IntersectionGeometry<typename Data::Intersection>(data.intersection(),data.intersectionIndex()),lfsu,x,lfsv,y);
  }

  const XL& x;
  YL& y;
};


template<typename XL, typename YL, typename Operator=CouplingOperator>
struct InvokeJacobianApplyCoupling
{

  InvokeJacobianApplyCoupling(const XL& xl, const XL& xn, YL& yl, YL& yn) :
    _local_x(xl),
    _remote_x(xn),
    _local_y(yl),
    _remote_y(yn)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains(),data.neighborSubDomains()))
      return;
    typedef typename Operator::template ExtractType<Child>::Type LOP;
    typedef typename Child::Traits::LocalSubProblem LocalSubProblem;
    typedef typename Child::Traits::RemoteSubProblem RemoteSubProblem;
    typedef typename LocalSubProblem::Traits::TrialLocalFunctionSpace LocalLFSU;
    typedef typename LocalSubProblem::Traits::TestLocalFunctionSpace LocalLFSV;
    typedef typename RemoteSubProblem::Traits::TrialLocalFunctionSpace RemoteLFSU;
    typedef typename RemoteSubProblem::Traits::TestLocalFunctionSpace RemoteLFSV;
    const LocalSubProblem& localSubProblem = child.localSubProblem();
    const RemoteSubProblem& remoteSubProblem = child.remoteSubProblem();
    LocalLFSU local_lfsu(data.lfsu(),localSubProblem,localSubProblem.trialGridFunctionSpaceConstraints());
    LocalLFSV local_lfsv(data.lfsv(),localSubProblem,localSubProblem.testGridFunctionSpaceConstraints());
    local_lfsu.bind();
    local_lfsv.bind();
    RemoteLFSU remote_lfsu(data.lfsun(),remoteSubProblem,remoteSubProblem.trialGridFunctionSpaceConstraints());
    RemoteLFSV remote_lfsv(data.lfsvn(),remoteSubProblem,remoteSubProblem.testGridFunctionSpaceConstraints());
    remote_lfsu.bind();
    remote_lfsv.bind();
    Operator::extract(child).jacobian_apply_coupling(IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                                       data.intersectionIndex()),
                                                     local_lfsu,
                                                     _local_x,
                                                     local_lfsv,
                                                     remote_lfsu,
                                                     _remote_x,
                                                     remote_lfsv,
                                                     _local_y,
                                                     _remote_y);
    data.setAlphaSkeletonInvoked();
  }

  const XL& _local_x;
  const XL& _remote_x;
  YL& _local_y;
  YL& _remote_y;
};


template<typename XL, typename YL, typename Operator=CouplingOperator>
struct InvokeJacobianApplyEnrichedCoupling
{

  InvokeJacobianApplyEnrichedCoupling(const XL& xl, const XL& xn, const XL& xc, YL& yl, YL& yn, YL& yc) :
    _local_x(xl),
    _remote_x(xn),
    _coupling_x(xc),
    _local_y(yl),
    _remote_y(yn),
    _coupling_y(yc)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains(),data.neighborSubDomains()))
      return;
    typedef typename Operator::template ExtractType<Child>::Type LOP;
    typedef typename Child::Traits::LocalSubProblem LocalSubProblem;
    typedef typename Child::Traits::RemoteSubProblem RemoteSubProblem;
    typedef typename LocalSubProblem::Traits::TrialLocalFunctionSpace LocalLFSU;
    typedef typename LocalSubProblem::Traits::TestLocalFunctionSpace LocalLFSV;
    typedef typename RemoteSubProblem::Traits::TrialLocalFunctionSpace RemoteLFSU;
    typedef typename RemoteSubProblem::Traits::TestLocalFunctionSpace RemoteLFSV;
    const LocalSubProblem& localSubProblem = child.localSubProblem();
    const RemoteSubProblem& remoteSubProblem = child.remoteSubProblem();
    LocalLFSU local_lfsu(data.lfsu(),localSubProblem,localSubProblem.trialGridFunctionSpaceConstraints());
    LocalLFSV local_lfsv(data.lfsv(),localSubProblem,localSubProblem.testGridFunctionSpaceConstraints());
    local_lfsu.bind();
    local_lfsv.bind();
    RemoteLFSU remote_lfsu(data.lfsun(),remoteSubProblem,remoteSubProblem.trialGridFunctionSpaceConstraints());
    RemoteLFSV remote_lfsv(data.lfsvn(),remoteSubProblem,remoteSubProblem.testGridFunctionSpaceConstraints());
    remote_lfsu.bind();
    remote_lfsv.bind();
    Operator::extract(child).jacobian_apply_enriched_coupling_first(IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                                                      data.intersectionIndex()),
                                                                    local_lfsu,
                                                                    _local_x,
                                                                    local_lfsv,
                                                                    child.couplingLocalFunctionSpace(data.couplinglfsu()),
                                                                    _coupling_x,
                                                                    child.couplingLocalFunctionSpace(data.couplinglfsv()),
                                                                    _local_y,
                                                                    _coupling_y);
    Operator::extract(child).jacobian_apply_enriched_coupling_second(IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                                                       data.intersectionIndex()),
                                                                     remote_lfsu,
                                                                     _remote_x,
                                                                     remote_lfsv,
                                                                     child.couplingLocalFunctionSpace(data.couplinglfsu()),
                                                                     _coupling_x,
                                                                     child.couplingLocalFunctionSpace(data.couplinglfsv()),
                                                                     _remote_y,
                                                                     _coupling_y);

    data.setAlphaSkeletonInvoked();
    data.setAlphaEnrichedCouplingInvoked();
  }

  const XL& _local_x;
  const XL& _remote_x;
  const XL& _coupling_x;
  YL& _local_y;
  YL& _remote_y;
  YL& _coupling_y;
};


template<typename XL, typename YL, typename Operator=SpatialOperator>
struct InvokeJacobianApplyVolumePostSkeleton
{

  InvokeJacobianApplyVolumePostSkeleton(const XL& xl, YL& yl) :
    x(xl),
    y(yl)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains()))
      return;
    typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
    typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
    LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsu.bind();
    lfsv.bind();
    Operator::extract(child).jacobian_apply_volume_post_skeleton(ElementGeometry<typename Data::Element>(data.element()),lfsu,x,lfsv,y);
  }

  const XL& x;
  YL& y;
};


template<typename XL, typename AL, typename Operator=SpatialOperator>
struct InvokeJacobianVolume
{

  InvokeJacobianVolume(const XL& xl, AL& al) :
    x(xl),
    a(al)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains()))
      return;
    typedef typename Child::Traits::LocalTrialFunctionSpace LFSU;
    typedef typename Child::Traits::LocalTestFunctionSpace LFSV;
    LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsu.bind();
    lfsv.bind();
    Operator::extract(child).jacobian_volume(ElementGeometry<typename Data::Element>(data.element()),lfsu,x,lfsv,a);
  }

  const XL& x;
  AL& a;
};


template<typename XL, typename AL, typename Operator=SpatialOperator>
struct InvokeJacobianSkeletonOrBoundary
{

  InvokeJacobianSkeletonOrBoundary(const XL& xl, const XL& xn, AL& al, AL& al_sn, AL& al_ns, AL& al_nn, bool applyOneSided) :
    _xl(xl),
    _xn(xn),
    _al(al),
    _al_sn(al_sn),
    _al_ns(al_ns),
    _al_nn(al_nn),
    _applyOneSided(applyOneSided)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains()))
      return;
    typedef typename Operator::template ExtractType<Child>::Type LOP;
    typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
    typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
    LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsu.bind();
    lfsv.bind();
    if (child.appliesTo(data.neighborSubDomains()))
      {
        if (_applyOneSided || LOP::doSkeletonTwoSided)
          {
            LFSU lfsun(data.lfsun(),child,child.trialGridFunctionSpaceConstraints());
            LFSV lfsvn(data.lfsvn(),child,child.testGridFunctionSpaceConstraints());
            lfsun.bind();
            lfsvn.bind();
            LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
              jacobian_skeleton(Operator::extract(child),
                                IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                  data.intersectionIndex()),
                                lfsu,_xl,lfsv,lfsun,_xn,lfsvn,_al,_al_sn,_al_ns,_al_nn);
            if(LOP::doAlphaSkeleton)
              data.setAlphaSkeletonInvoked();
          }
      }
    else
      {
        LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          jacobian_boundary(Operator::extract(child),
                            IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                              data.intersectionIndex()),
                            lfsu,_xl,lfsv,_al);
      }
  }

  const XL& _xl;
  const XL& _xn;
  AL& _al;
  AL& _al_sn;
  AL& _al_ns;
  AL& _al_nn;
  const bool _applyOneSided;
};


template<typename XL, typename AL, typename Operator=SpatialOperator>
struct InvokeJacobianBoundary
{

  InvokeJacobianBoundary(const XL& xl, AL& al) :
    x(xl),
    a(al)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains()))
      return;
    typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
    typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
    LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsu.bind();
    lfsv.bind();
    Operator::extract(child).jacobian_boundary(IntersectionGeometry<typename Data::Intersection>(data.intersection(),data.intersectionIndex()),lfsu,x,lfsv,a);
  }

  const XL& x;
  AL& a;
};


template<typename XL, typename AL, typename Operator=CouplingOperator>
struct InvokeJacobianCoupling
{

  InvokeJacobianCoupling(const XL& xl, const XL& xn, AL& al, AL& al_sn, AL& al_ns, AL& al_nn) :
    _local_x(xl),
    _remote_x(xn),
    _local_a(al),
    _local_to_remote_a(al_sn),
    _remote_to_local_a(al_ns),
    _remote_a(al_nn)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains(),data.neighborSubDomains()))
      return;
    typedef typename Operator::template ExtractType<Child>::Type LOP;
    typedef typename Child::Traits::LocalSubProblem LocalSubProblem;
    typedef typename Child::Traits::RemoteSubProblem RemoteSubProblem;
    typedef typename LocalSubProblem::Traits::TrialLocalFunctionSpace LocalLFSU;
    typedef typename LocalSubProblem::Traits::TestLocalFunctionSpace LocalLFSV;
    typedef typename RemoteSubProblem::Traits::TrialLocalFunctionSpace RemoteLFSU;
    typedef typename RemoteSubProblem::Traits::TestLocalFunctionSpace RemoteLFSV;
    const LocalSubProblem& localSubProblem = child.localSubProblem();
    const RemoteSubProblem& remoteSubProblem = child.remoteSubProblem();
    LocalLFSU local_lfsu(data.lfsu(),localSubProblem,localSubProblem.trialGridFunctionSpaceConstraints());
    LocalLFSV local_lfsv(data.lfsv(),localSubProblem,localSubProblem.testGridFunctionSpaceConstraints());
    local_lfsu.bind();
    local_lfsv.bind();
    RemoteLFSU remote_lfsu(data.lfsun(),remoteSubProblem,remoteSubProblem.trialGridFunctionSpaceConstraints());
    RemoteLFSV remote_lfsv(data.lfsvn(),remoteSubProblem,remoteSubProblem.testGridFunctionSpaceConstraints());
    remote_lfsu.bind();
    remote_lfsv.bind();
    Operator::extract(child).jacobian_coupling(IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                                 data.intersectionIndex()),
                                               local_lfsu,
                                               _local_x,
                                               local_lfsv,
                                               remote_lfsu,
                                               _remote_x,
                                               remote_lfsv,
                                               _local_a,
                                               _local_to_remote_a,
                                               _remote_to_local_a,
                                               _remote_a);
    data.setAlphaSkeletonInvoked();
  }

  const XL& _local_x;
  const XL& _remote_x;
  AL& _local_a;
  AL& _local_to_remote_a;
  AL& _remote_to_local_a;
  AL& _remote_a;
};


template<typename XL, typename AL, typename Operator=CouplingOperator>
struct InvokeJacobianEnrichedCoupling
{

  InvokeJacobianEnrichedCoupling(const XL& xl, const XL& xn, const XL& xc, AL& al, AL& al_nn, AL& al_sc, AL&  al_cs, AL& al_nc, AL& al_cn, AL& al_cc) :
    _local_x(xl),
    _remote_x(xn),
    _coupling_x(xc),
    _local_a(al),
    _remote_a(al_nn),
    _local_to_coupling_a(al_sc),
    _coupling_to_local_a(al_cs),
    _remote_to_coupling_a(al_nc),
    _coupling_to_remote_a(al_cn),
    _coupling_a(al_cc)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains(),data.neighborSubDomains()))
      return;
    typedef typename Operator::template ExtractType<Child>::Type LOP;
    typedef typename Child::Traits::LocalSubProblem LocalSubProblem;
    typedef typename Child::Traits::RemoteSubProblem RemoteSubProblem;
    typedef typename LocalSubProblem::Traits::TrialLocalFunctionSpace LocalLFSU;
    typedef typename LocalSubProblem::Traits::TestLocalFunctionSpace LocalLFSV;
    typedef typename RemoteSubProblem::Traits::TrialLocalFunctionSpace RemoteLFSU;
    typedef typename RemoteSubProblem::Traits::TestLocalFunctionSpace RemoteLFSV;
    const LocalSubProblem& localSubProblem = child.localSubProblem();
    const RemoteSubProblem& remoteSubProblem = child.remoteSubProblem();
    LocalLFSU local_lfsu(data.lfsu(),localSubProblem,localSubProblem.trialGridFunctionSpaceConstraints());
    LocalLFSV local_lfsv(data.lfsv(),localSubProblem,localSubProblem.testGridFunctionSpaceConstraints());
    local_lfsu.bind();
    local_lfsv.bind();
    RemoteLFSU remote_lfsu(data.lfsun(),remoteSubProblem,remoteSubProblem.trialGridFunctionSpaceConstraints());
    RemoteLFSV remote_lfsv(data.lfsvn(),remoteSubProblem,remoteSubProblem.testGridFunctionSpaceConstraints());
    remote_lfsu.bind();
    remote_lfsv.bind();
    Operator::extract(child).jacobian_enriched_coupling_first(IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                                                data.intersectionIndex()),
                                                              local_lfsu,
                                                              _local_x,
                                                              local_lfsv,
                                                              child.couplingLocalFunctionSpace(data.couplinglfsu()),
                                                              _coupling_x,
                                                              child.couplingLocalFunctionSpace(data.couplinglfsv()),
                                                              _local_a,
                                                              _local_to_coupling_a,
                                                              _coupling_to_local_a,
                                                              _coupling_a);
    Operator::extract(child).jacobian_enriched_coupling_second(IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                                                 data.intersectionIndex()),
                                                               remote_lfsu,
                                                               _remote_x,
                                                               remote_lfsv,
                                                               child.couplingLocalFunctionSpace(data.couplinglfsu()),
                                                               _coupling_x,
                                                               child.couplingLocalFunctionSpace(data.couplinglfsv()),
                                                                _remote_a,
                                                               _remote_to_coupling_a,
                                                               _coupling_to_remote_a,
                                                               _coupling_a);
    data.setAlphaSkeletonInvoked();
    data.setAlphaEnrichedCouplingInvoked();
  }

  const XL& _local_x;
  const XL& _remote_x;
  const XL& _coupling_x;
  AL& _local_a;
  AL& _remote_a;
  AL& _local_to_coupling_a;
  AL& _coupling_to_local_a;
  AL& _remote_to_coupling_a;
  AL& _coupling_to_remote_a;
  AL& _coupling_a;
};


template<typename XL, typename AL, typename Operator=SpatialOperator>
struct InvokeJacobianVolumePostSkeleton
{

  InvokeJacobianVolumePostSkeleton(const XL& xl, AL& al) :
    x(xl),
    a(al)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains()))
      return;
    typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
    typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
    LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsu.bind();
    lfsv.bind();
    Operator::extract(child).jacobian_volume_post_skeleton(ElementGeometry<typename Data::Element>(data.element()),lfsu,x,lfsv,a);
  }

  const XL& x;
  AL& a;
};

}
}
}

#endif // DUNE_MULTIDOMAIN_MULTIDOMAINGRIDOPERATORSPACEINVOCATIONHELPERS_HH
