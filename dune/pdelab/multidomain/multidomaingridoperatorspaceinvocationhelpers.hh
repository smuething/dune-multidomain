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
                           lfsv.globalIndex(lfsv.localIndex(localpattern[k].i())),
                           lfsu.globalIndex(lfsu.localIndex(localpattern[k].j()))
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
                           data.lfsv().globalIndex(data.lfsv().localIndex(localpattern_sn[k].i())),
                           data.lfsun().globalIndex(data.lfsun().localIndex(localpattern_sn[k].j()))
                           );

    for (size_t k=0; k<localpattern_ns.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.lfsvn().globalIndex(data.lfsvn().localIndex(localpattern_ns[k].i())),
                           data.lfsu().globalIndex(data.lfsu().localIndex(localpattern_ns[k].j()))
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
                           data.lfsv().globalIndex(data.lfsv().localIndex(localpattern_sn[k].i())),
                           data.lfsun().globalIndex(data.lfsun().localIndex(localpattern_sn[k].j()))
                           );
    for (size_t k=0; k<localpattern_ns.size(); ++k)
      data.gos().add_entry(globalpattern,
                           data.lfsvn().globalIndex(data.lfsvn().localIndex(localpattern_ns[k].i())),
                           data.lfsu().globalIndex(data.lfsu().localIndex(localpattern_ns[k].j()))
                           );
  }

  P& globalpattern;
};


template<typename XL, typename RL, typename Operator=SpatialOperator>
struct InvokeAlphaVolume
{

  InvokeAlphaVolume(const XL& xl, RL& rl) :
    x(xl),
    r(rl)
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
    Operator::extract(child).alpha_volume(ElementGeometry<typename Data::Element>(data.element()),lfsu,x,lfsv,r);
  }

  const XL& x;
  RL& r;
};


template<typename RL, typename Operator=SpatialOperator>
struct InvokeLambdaVolume
{

  InvokeLambdaVolume(RL& rl) :
    r(rl)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains()))
      return;
    typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsv.bind();
    Operator::extract(child).lambda_volume(ElementGeometry<typename Data::Element>(data.element()),lfsv,r);
  }

  RL& r;
};


template<typename XL, typename RL, typename Operator=SpatialOperator>
struct InvokeAlphaSkeletonOrBoundary
{

  InvokeAlphaSkeletonOrBoundary(const XL& xl, const XL& xn, RL& rl, RL& rn, bool applyOneSided) :
    _xl(xl),
    _xn(xn),
    _rl(rl),
    _rn(rn),
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
            LocalAssemblerCallSwitch<Child,LOP::doAlphaSkeleton>::
              alpha_skeleton(Operator::extract(child),
                             IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                               data.intersection_index),
                             lfsu,_xl,lfsv,lfsun,_xn,lfsvn,_rl,_rn);
            if(LOP::doAlphaSkeleton)
              data.setAlphaSkeletonInvoked();
          }
      }
    else
      {
        LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          alpha_boundary(Operator::extract(child),
                         IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                           data.intersection_index),
                         lfsu,_xl,lfsv,_rl);
        LocalAssemblerCallSwitch<LOP,LOP::doLambdaBoundary>::
          lambda_boundary(Operator::extract(child),
                          IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                            data.intersection_index),
                          lfsv,_rl);
      }
  }

  const XL& _xl;
  const XL& _xn;
  RL& _rl;
  RL& _rn;
  const bool _applyOneSided;
};


template<typename XL, typename RL, typename Operator=SpatialOperator>
struct InvokeAlphaBoundary
{

  InvokeAlphaBoundary(const XL& xl, RL& rl) :
    x(xl),
    r(rl)
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
    Operator::extract(child).alpha_boundary(IntersectionGeometry<typename Data::Intersection>(data.intersection(),data.intersectionIndex()),lfsu,x,lfsv,r);
  }

  const XL& x;
  RL& r;
};


template<typename XL, typename RL, typename Operator=CouplingOperator>
struct InvokeAlphaCoupling
{

  InvokeAlphaCoupling(const XL& xl, const XL& xn, RL& rl, RL& rn) :
    _local_x(xl),
    _remote_x(xn),
    _local_r(rl),
    _remote_r(rn)
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
    Operator::extract(child).alpha_coupling(IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                              data.intersectionIndex()),
                                            local_lfsu,
                                            _local_x,
                                            local_lfsv,
                                            remote_lfsu,
                                            _remote_x,
                                            remote_lfsv,
                                            _local_r,
                                            _remote_r);
    data.setAlphaSkeletonInvoked();
  }

  const XL& _local_x;
  const XL& _remote_x;
  RL& _local_r;
  RL& _remote_r;
};


template<typename RL, typename Operator=SpatialOperator>
struct InvokeLambdaBoundary
{

  InvokeLambdaBoundary(RL& rl) :
    r(rl)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains()))
      return;
    typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsv.bind();
    Operator::extract(child).lambda_boundary(IntersectionGeometry<typename Data::Intersection>(data.intersection(),data.intersectionIndex()),lfsv,r);
  }

  RL& r;
};


template<typename XL, typename RL, typename Operator=SpatialOperator>
struct InvokeAlphaVolumePostSkeleton
{

  InvokeAlphaVolumePostSkeleton(const XL& xl, RL& rl) :
    x(xl),
    r(rl)
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
    Operator::extract(child).alpha_volume_post_skeleton(ElementGeometry<typename Data::Element>(data.element()),lfsu,x,lfsv,r);
  }

  const XL& x;
  RL& r;
};


template<typename RL, typename Operator=SpatialOperator>
struct InvokeLambdaVolumePostSkeleton
{

  InvokeLambdaVolumePostSkeleton(RL& rl) :
    r(rl)
  {}

  template<typename Data, typename Child>
  void operator()(Data& data, const Child& child)
  {
    if (!child.appliesTo(data.elementSubDomains()))
      return;
    typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
    LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
    lfsv.bind();
    Operator::extract(child).lambda_volume_post_skeleton(ElementGeometry<typename Data::Element>(data.element()),lfsv,r);
  }

  RL& r;
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