#ifndef DUNE_PDELAB_MULTIDOMAIN_JACOBIANASSEMBLERFUNCTORS_HH
#define DUNE_PDELAB_MULTIDOMAIN_JACOBIANASSEMBLERFUNCTORS_HH

namespace Dune {
namespace PDELab {
namespace MultiDomain {


template<typename Data>
struct invoke_jacobian_volume
  : public data_accessor<Data>
{

  using data_accessor<Data>::data;

  template<typename SubProblem>
  void operator()(const SubProblem& subProblem)
  {
    if (!subProblem.appliesTo(data().eg()))
      return;
    Data::Operator::extract(subProblem).jacobian_volume(data().eg(),
                                                        LFS::lfsu(data(),subProblem),
                                                        data().x(),
                                                        LFS::lfsv(data(),subProblem),
                                                        data().a());
  }

};


template<typename Data>
struct invoke_jacobian_skeleton_or_boundary
  : public data_accessor<Data>
{

  using data_accessor<Data>::data;

  template<typename SubProblem>
  void operator()(const SubProblem& subProblem)
  {
    if (!subProblem.appliesTo(data().ig().insideElement()))
      return;
    if (subProblem.appliesTo(data().ig().outsideElement()))
      {
        if (_applyOneSided || LOP::doSkeletonTwoSided)
          {
            LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
              jacobian_skeleton(Data::Operator::extract(child),data().ig()
                                LFS::lfsu_s(data(),subProblem),data().x_s(),LFS::lfsv_s(data(),subProblem),
                                LFS::lfsu_n(data(),subProblem),data().x_n(),LFS::lfsv_n(data(),subProblem),
                                data().a_ss(),data().a_sn(),data().a_ns(),data().a_nn());
            if(LOP::doAlphaSkeleton)
              data().neighbor_accessed() = true;
          }
      }
    else
      {
        LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
          jacobian_boundary(Data::Operator::extract(child),data().ig(),
                            LFS::lfsu_s(data(),subProblem),data().x_s(),LFS::lfsv_s(data(),subProblem),
                            data().a_ss());
      }
  }

};


template<typename Data>
struct invoke_jacobian_boundary
  : public data_accessor<Data>
{

  using data_accessor<Data>::data;

  template<typename SubProblem>
  void operator()(const SubProblem& subProblem)
  {
    if (!subProblem.appliesTo(data().ig().insideElement()))
      return;
    Data::Operator::extract(subProblem).jacobian_boundary(data().ig(),
                                                          LFS::lfsu(data(),subProblem),data().x(),LFS::lfsv(data(),subProblem),
                                                          data().a());
  }

};


template<typename Data>
struct invoke_jacobian_coupling
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
    Data::Operator::extract(subProblem).jacobian_coupling(data().ig(),
                                                          LFS::lfsu_s(data(),localSubProblem),
                                                          data().x_s(),
                                                          LFS::lfsv_s(data(),localSubProblem),
                                                          LFS::lfsu_n(data(),remoteSubProblem),
                                                          data().x_n(),
                                                          LFS::lfsv_n(data(),remoteSubProblem),
                                                          data().a_ss(),data().a_sn(),data().a_ns(),data.a_n());
    data().neighbor_accessed() = true;
  }

};


template<typename Data>
struct invoke_jacobian_enriched_coupling
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
    typename LFS::LFSU_C<Data,SubProblem> lfsu_c = LFS::lfsu_c(data(),coupling);
    typename LFS::LFSV_C<Data,SubProblem> lfsv_c = LFS::lfsv_c(data(),coupling);
    Operator::extract(coupling).jacobian_enriched_coupling_first(data().ig(),
                                                                 LFS::lfsu_s(data(),localSubProblem),
                                                                 data().x_s(),
                                                                 LFS::lfsv_s(data(),localSubProblem),
                                                                 lfsu_c,data().x_c(),lfsv_c,
                                                                 data().a_ss(),data().a_sc(),data().a_cs(),data().a_cc());
    Operator::extract(coupling).jacobian_enriched_coupling_second(data().ig(),
                                                                  LFS::lfsu_n(data(),remoteSubProblem),
                                                                  data().x_n(),
                                                                  LFS::lfsv_n(data(),remoteSubProblem),
                                                                  lfsu_c,data().x_c(),lfsv_c,
                                                                  data().a_nn(),data().a_nc(),data().a_cn(),data().a_cc());
    data().neighbor_accessed() = true;
    data().coupling_accessed() = true;
  }

};


template<typename Data>
struct invoke_jacobian_volume_post_skeleton
  : public data_accessor<Data>
{

  using data_accessor<Data>::data;

  template<typename SubProblem>
  void operator()(const SubProblem& subProblem)
  {
    if (!subProblem.appliesTo(data().ig().insideElement()))
      return;
    Operator::extract(subProblem).jacobian_volume_post_skeleton(data().eg(),
                                                                LFS::lfsu(data(),subProblem),
                                                                data().x(),
                                                                LFS::lfsv(data(),subProblem),
                                                                data().a());
  }

};

}
}
}

#endif // DUNE_PDELAB_MULTIDOMAIN_JACOBIANASSEMBLERFUNCTORS_HH
