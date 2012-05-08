#ifndef DUNE_PDELAB_MULTIDOMAIN_RESIDUALASSEMBLERFUNCTORS_HH
#define DUNE_PDELAB_MULTIDOMAIN_RESIDUALASSEMBLERFUNCTORS_HH

#include <dune/pdelab/multidomain/datawrappers.hh>
#include <dune/pdelab/multidomain/visitor.hh>
#include <dune/pdelab/multidomain/operatorcallguards.hh>
#include <dune/pdelab/multidomain/localfunctionspaceutility.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include <dune/pdelab/constraints/constraints.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {

namespace functors {

  template<typename Data>
  struct alpha_volume
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename SubProblem>
    void operator()(const SubProblem& subProblem)
    {
      if (!subProblem.appliesTo(data().eg()))
        return;
      Data::Operator::extract(subProblem).alpha_volume(data().eg(),
                                                       LFS::lfsu(data(),subProblem),
                                                       data().x(),
                                                       LFS::lfsv(data(),subProblem),
                                                       data().r());
    }

  };


  template<typename Data>
  struct lambda_volume
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename SubProblem>
    void operator()(const SubProblem& subProblem)
    {
      if (!subProblem.appliesTo(data().eg()))
        return;
      Data::Operator::extract(subProblem).lambda_volume(data().eg(),
                                                        LFS::lfsv(data(),subProblem),
                                                        data().r());
    }

  };


  template<typename Data>
  struct alpha_skeleton_or_boundary
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename SubProblem>
    void operator()(const SubProblem& subProblem)
    {
      typedef typename Data::Operator::template ExtractType<SubProblem>::Type LOP;
      if (!subProblem.appliesTo(data().ig().insideElement()))
        return;
      if (subProblem.appliesTo(data().ig().outsideElement()))
        {
          if (!(LOP::doSkeletonTwoSided || data().ig().oneSidedDirection()))
            return;
          guarded_call::alpha_skeleton(Data::Operator::extract(subProblem),data().ig(),
                                       LFS::lfsu_s(data(),subProblem),data().x_s(),LFS::lfsv_s(data(),subProblem),
                                       LFS::lfsu_n(data(),subProblem),data().x_n(),LFS::lfsv_n(data(),subProblem),
                                       data().r_s(),data().r_n());
        }
      else
        {
          guarded_call::alpha_boundary(Data::Operator::extract(subProblem),
                                       data().ig(),
                                       LFS::lfsu_s(data(),subProblem),data().x_s(),LFS::lfsv_s(data(),subProblem),
                                       data().r_s());
        }
    }
  };


  template<typename Data>
  struct lambda_skeleton_or_boundary
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename SubProblem>
    void operator()(const SubProblem& subProblem)
    {
      typedef typename Data::Operator::template ExtractType<SubProblem>::Type LOP;
      if (!subProblem.appliesTo(data().ig().insideElement()))
        return;
      if (subProblem.appliesTo(data().ig().outsideElement()))
        {
          if (!(LOP::doSkeletonTwoSided || data().ig().oneSidedDirection()))
            return;
          guarded_call::lambda_skeleton(Data::Operator::extract(subProblem),data().ig(),
                                        LFS::lfsv_s(data(),subProblem),
                                        LFS::lfsv_n(data(),subProblem),
                                        data().r_s(),data().r_n());
        }
      else
        {
          guarded_call::lambda_boundary(Data::Operator::extract(subProblem),data().ig(),
                                        LFS::lfsv_s(data(),subProblem),
                                        data().r_s());
        }
    }
  };


  template<typename Data>
  struct alpha_boundary
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename SubProblem>
    void operator()(const SubProblem& subProblem)
    {
      if (!subProblem.appliesTo(data().ig().insideElement()))
        return;
      Data::Operator::extract(subProblem).alpha_boundary(data().ig(),
                                                         LFS::lfsu(data(),subProblem),
                                                         data().x(),
                                                         LFS::lfsv(data(),subProblem),
                                                         data().r());
    }

  };


  template<typename Data>
  struct lambda_boundary
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename SubProblem>
    void operator()(const SubProblem& subProblem)
    {
      if (!subProblem.appliesTo(data().ig().insideElement()))
        return;
      Data::Operator::extract(subProblem).lambda_boundary(data().ig(),
                                                          LFS::lfsv(data(),subProblem),
                                                          data().r());
    }

  };


  template<typename Data>
  struct alpha_coupling
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

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
                                                       data().x_s(),
                                                       LFS::lfsv_s(data(),localSubProblem),
                                                       LFS::lfsu_n(data(),remoteSubProblem),
                                                       data().x_n(),
                                                       LFS::lfsv_n(data(),remoteSubProblem),
                                                       data().r_s(),data().r_n());
    }

  };


  template<typename Data>
  struct lambda_coupling
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

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
                                                       LFS::lfsv_s(data(),localSubProblem),
                                                       LFS::lfsv_n(data(),remoteSubProblem),
                                                       data().r_s(),data().r_n());
    }

  };


  template<typename Data>
  struct alpha_enriched_coupling
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

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
      typename LFS::LFSU_C<Data,Coupling>::type lfsu_c = LFS::lfsu_c(data(),coupling);
      typename LFS::LFSV_C<Data,Coupling>::type lfsv_c = LFS::lfsv_c(data(),coupling);
      guarded_call::alpha_enriched_coupling_first(Data::Operator::extract(coupling),
                                                  data().ig(),
                                                  LFS::lfsu_s(data(),localSubProblem),
                                                  data().x_s(),
                                                  LFS::lfsv_s(data(),localSubProblem),
                                                  lfsu_c,data().x_c(),lfsv_c,
                                                  data().r_s(),data().r_c());
      guarded_call::alpha_enriched_coupling_second(Data::Operator::extract(coupling),
                                                   data().ig(),
                                                   LFS::lfsu_n(data(),remoteSubProblem),
                                                   data().x_n(),
                                                   LFS::lfsv_n(data(),remoteSubProblem),
                                                   lfsu_c,data().x_c(),lfsv_c,
                                                   data().r_n(),data().r_c());
      guarded_call::alpha_enriched_coupling(Data::Operator::extract(coupling),
                                            data().ig(),
                                            lfsu_c,data().x_c(),lfsv_c,
                                            data().r_c());
    }

  };


  template<typename Data>
  struct lambda_enriched_coupling
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

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
      typename LFS::LFSV_C<Data,Coupling> lfsv_c = LFS::lfsv_c(data(),coupling);
      guarded_call::lambda_enriched_coupling_first(Data::Operator::extract(coupling),
                                                   data().ig(),
                                                   LFS::lfsv_s(data(),localSubProblem),
                                                   lfsv_c,
                                                   data().r_s(),data().r_c());
      guarded_call::lambda_enriched_coupling_second(Data::Operator::extract(coupling),
                                                    data().ig(),
                                                    LFS::lfsv_n(data(),remoteSubProblem),
                                                    lfsv_c,
                                                    data().r_n(),data().r_c());
      guarded_call::lambda_enriched_coupling(Data::Operator::extract(coupling),
                                             data().ig(),
                                             lfsv_c,
                                             data().r_c());
    }

  };


  template<typename Data>
  struct alpha_volume_post_skeleton
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename SubProblem>
    void operator()(const SubProblem& subProblem)
    {
      if (!subProblem.appliesTo(data().eg()))
        return;
      Data::Operator::extract(subProblem).alpha_volume_post_skeleton(data().eg(),
                                                                     LFS::lfsu(data(),subProblem),
                                                                     data().x(),
                                                                     LFS::lfsv(data(),subProblem),
                                                                     data().r());
    }

  };


  template<typename Data>
  struct lambda_volume_post_skeleton
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename SubProblem>
    void operator()(const SubProblem& subProblem)
    {
      if (!subProblem.appliesTo(data().eg()))
        return;
      Data::Operator::extract(subProblem).lambda_volume_post_skeleton(data().eg(),
                                                                      LFS::lfsv(data,subProblem),
                                                                      data().r());
    }

  };

}

}
}
}

#endif // DUNE_PDELAB_MULTIDOMAIN_RESIDUALASSEMBLERFUNCTORS_HH
