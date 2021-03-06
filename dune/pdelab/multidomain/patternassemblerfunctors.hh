#ifndef DUNE_PDELAB_MULTIDOMAIN_PATTERNASSEMBLERFUNCTORS_HH
#define DUNE_PDELAB_MULTIDOMAIN_PATTERNASSEMBLERFUNCTORS_HH

#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/multidomain/datawrappers.hh>
#include <dune/pdelab/multidomain/localfunctionspaceutility.hh>
#include <dune/pdelab/multidomain/operatorcallguards.hh>
#include <dune/pdelab/multidomain/visitor.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {

namespace functors {

  template<typename Data>
  struct invoke_pattern_volume
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename SubProblem>
    void operator()(const SubProblem& subProblem)
    {
      if (!subProblem.appliesTo(data().eg()))
        return;
      Data::Operator::extract(subProblem).pattern_volume(LFS::lfsu(data(),subProblem),
                                                         LFS::lfsv(data(),subProblem),
                                                         data().pattern_ss());
    }

  };


  template<typename Data>
  struct invoke_pattern_skeleton_or_boundary
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
          guarded_call::pattern_skeleton(Data::Operator::extract(subProblem),
                                         LFS::lfsu_s(data(),subProblem),LFS::lfsv_s(data(),subProblem),
                                         LFS::lfsu_n(data(),subProblem),LFS::lfsv_n(data(),subProblem),
                                         data().pattern_sn(),data().pattern_ns());
        }
      else
        {
          guarded_call::pattern_boundary(Data::Operator::extract(subProblem),
                                         LFS::lfsu_s(data(),subProblem),LFS::lfsv_s(data(),subProblem),
                                         data().pattern_ss());
        }
    }

  };


  template<typename Data>
  struct invoke_pattern_boundary
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename SubProblem>
    void operator()(const SubProblem& subProblem)
    {
      if (!subProblem.appliesTo(data().ig().insideElement()))
        return;
      Data::Operator::extract(subProblem).pattern_boundary(LFS::lfsu(data(),subProblem),LFS::lfsv(data(),subProblem),
                                                           data().pattern_ss());
    }

  };


  template<typename Data>
  struct invoke_pattern_coupling
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Coupling>
    void operator()(const Coupling& coupling)
    {
      if (!coupling.appliesTo(data().ig()))
        return;
      typedef typename Coupling::Traits::LocalSubProblem LocalSubProblem;
      typedef typename Coupling::Traits::RemoteSubProblem RemoteSubProblem;
      const LocalSubProblem& localSubProblem = coupling.localSubProblem();
      const RemoteSubProblem& remoteSubProblem = coupling.remoteSubProblem();
      Data::Operator::extract(coupling).pattern_coupling(LFS::lfsu_s(data(),localSubProblem),
                                                         LFS::lfsv_s(data(),localSubProblem),
                                                         LFS::lfsu_n(data(),remoteSubProblem),
                                                         LFS::lfsv_n(data(),remoteSubProblem),
                                                         data().pattern_ss(),data().pattern_sn(),
                                                         data().pattern_ns(),data().pattern_nn());
    }

  };


  template<typename Data>
  struct invoke_pattern_enriched_coupling
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Coupling>
    void operator()(const Coupling& coupling)
    {
      if (!coupling.appliesTo(data().ig()))
        return;
      typedef typename Coupling::Traits::LocalSubProblem LocalSubProblem;
      typedef typename Coupling::Traits::RemoteSubProblem RemoteSubProblem;
      const LocalSubProblem& localSubProblem = coupling.localSubProblem();
      const RemoteSubProblem& remoteSubProblem = coupling.remoteSubProblem();
      typename LFS::LFSU_C<Data,Coupling>::type lfsu_c = LFS::lfsu_c(data(),coupling);
      typename LFS::LFSV_C<Data,Coupling>::type lfsv_c = LFS::lfsv_c(data(),coupling);
      guarded_call::pattern_enriched_coupling_first(Data::Operator::extract(coupling),
                                                    LFS::lfsu_s(data(),localSubProblem),
                                                    LFS::lfsv_s(data(),localSubProblem),
                                                    lfsu_c,lfsv_c,
                                                    data().pattern_ss(),data().pattern_sc(),
                                                    data().pattern_cs(),data().pattern_cc());
      guarded_call::pattern_enriched_coupling_second(Data::Operator::extract(coupling),
                                                     LFS::lfsu_n(data(),remoteSubProblem),
                                                     LFS::lfsv_n(data(),remoteSubProblem),
                                                     lfsu_c,lfsv_c,
                                                     data().pattern_nn(),data().pattern_nc(),
                                                     data().pattern_cn(),data().pattern_cc());
      guarded_call::pattern_enriched_coupling(Data::Operator::extract(coupling),
                                              lfsu_c,lfsv_c,
                                              data().pattern_cc());
    }

  };


  template<typename Data>
  struct invoke_pattern_volume_post_skeleton
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename SubProblem>
    void operator()(const SubProblem& subProblem)
    {
      if (!subProblem.appliesTo(data().eg()))
        return;
      Data::Operator::extract(subProblem).pattern_volume_post_skeleton(LFS::lfsu(data(),subProblem),
                                                                       LFS::lfsv(data(),subProblem),
                                                                       data().pattern_ss());
    }

  };

}

}
}
}

#endif // DUNE_PDELAB_MULTIDOMAIN_PATTERNASSEMBLERFUNCTORS_HH
