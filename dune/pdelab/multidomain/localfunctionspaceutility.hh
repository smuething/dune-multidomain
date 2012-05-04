#ifndef DUNE_PDELAB_MULTIDOMAIN_LOCALFUNCTIONSPACEUTILITY_HH
#define DUNE_PDELAB_MULTIDOMAIN_LOCALFUNCTIONSPACEUTILITY_HH

namespace Dune {
namespace PDELab {
namespace MultiDomain {

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
    return LFS(data.lfsu_s_cache().localFunctionSpace(),subProblem);
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
    return LFS(data.lfsv_s_cache().localFunctionSpace(),subProblem);
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
    return LFS(data.lfsu_s_cache().localFunctionSpace(),subProblem);
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
    return LFS(data.lfsv_s_cache().localFunctionSpace(),subProblem);
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
    return LFS(data.lfsu_n_cache().localFunctionSpace(),subProblem);
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
    return LFS(data.lfsv_n_cache().localFunctionSpace(),subProblem);
  }

  template<typename Data, typename Coupling>
  struct LFSU_C
  {
    typedef const typename Coupling::template CouplingLocalFunctionSpace<typename Data::LFSU_C_Cache::LocalFunctionSpace>::Type& type;
  };

  template<typename Data, typename Coupling>
  inline typename LFSU_C<Data,Coupling>::type
  lfsu_c(const Data& data, const Coupling& coupling)
  {
    return coupling.couplingLocalFunctionSpace(data.lfsu_c_cache().localFunctionSpace());
  }

  template<typename Data, typename Coupling>
  struct LFSV_C
  {
    typedef const typename Coupling::template CouplingLocalFunctionSpace<typename Data::LFSV_C_Cache::LocalFunctionSpace>::Type& type;
  };

  template<typename Data, typename Coupling>
  inline typename LFSV_C<Data,Coupling>::type
  lfsv_c(const Data& data, const Coupling& coupling)
  {
    return coupling.couplingLocalFunctionSpace(data.lfsv_c_cache().localFunctionSpace());
  }

} // namespace LFS

}
}
}

#endif // DUNE_PDELAB_MULTIDOMAIN_LOCALFUNCTIONSPACEUTILITY_HH
