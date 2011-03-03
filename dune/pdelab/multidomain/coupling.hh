#ifndef DUNE_MULTIDOMAIN_COUPLING_HH
#define DUNE_MULTIDOMAIN_COUPLING_HH

#include <dune/pdelab/common/typetree.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {

template<
  typename LocalSubProblem_,
  typename RemoteSubProblem_,
  typename CouplingOperator_
  >
struct CouplingTraits
{
  typedef LocalSubProblem_ LocalSubProblem;
  typedef RemoteSubProblem_ RemoteSubProblem;
  typedef CouplingOperator_ CouplingOperator;
};


struct CouplingTag;

template<
  typename LocalSubProblem,
  typename RemoteSubProblem,
  typename CouplingOperator
  >
class Coupling
  : public Dune::PDELab::TypeTree::LeafNode
{

public:

  typedef CouplingTag MultiDomainComponentTag;

  typedef CouplingTraits<LocalSubProblem,RemoteSubProblem,CouplingOperator> Traits;

  Coupling(const LocalSubProblem& localSubProblem,
           const RemoteSubProblem& remoteSubProblem,
           CouplingOperator& couplingOperator)
    : _localSubProblem(localSubProblem)
    , _remoteSubProblem(remoteSubProblem)
    , _operator(couplingOperator)
  {}

  const CouplingOperator& couplingOperator() const {
    return _operator;
  }

  const LocalSubProblem& localSubProblem() const {
    return _localSubProblem;
  }

  const RemoteSubProblem& remoteSubProblem() const {
    return _remoteSubProblem;
  }

  template<typename SDS1, typename SDS2>
  bool appliesTo(const SDS1& localSubDomains, const SDS2& remoteSubDomains) const {
    return localSubProblem().appliesTo(localSubDomains) && remoteSubProblem().appliesTo(remoteSubDomains);
  }

private:
  const LocalSubProblem& _localSubProblem;
  const RemoteSubProblem& _remoteSubProblem;

protected:
  CouplingOperator& _operator;
};


template<
  typename LocalSubProblem,
  typename RemoteSubProblem,
  typename CouplingOperator
  >
struct is_coupling<Coupling<LocalSubProblem,RemoteSubProblem,CouplingOperator> >
{
  static const bool value = true;
};


template<
  typename LocalSubProblem_,
  typename RemoteSubProblem_,
  typename CouplingOperator_,
  std::size_t couplingLFSIndex_
  >
struct EnrichedCouplingTraits
  : public CouplingTraits<LocalSubProblem_,
                          RemoteSubProblem_,
                          CouplingOperator_
                          >
{
  static const std::size_t couplingLFSIndex = couplingLFSIndex_;
};


struct EnrichedCouplingTag;

template<
  typename LocalSubProblem,
  typename RemoteSubProblem,
  typename CouplingOperator,
  std::size_t couplingLFSIndex
  >
class EnrichedCoupling
  : public Dune::PDELab::TypeTree::LeafNode
{

public:

  typedef EnrichedCouplingTag MultiDomainComponentTag;

  typedef EnrichedCouplingTraits<LocalSubProblem,RemoteSubProblem,CouplingOperator,couplingLFSIndex> Traits;

  EnrichedCoupling(const LocalSubProblem& localSubProblem,
                   const RemoteSubProblem& remoteSubProblem,
                   CouplingOperator& couplingOperator)
    : _localSubProblem(localSubProblem)
    , _remoteSubProblem(remoteSubProblem)
    , _operator(couplingOperator)
  {}

  template<typename LFS>
  struct CouplingLocalFunctionSpace
  {
    typedef typename LFS::template Child<Traits::couplingLFSIndex>::Type Type;
  };

  template<typename LFS>
  static const typename CouplingLocalFunctionSpace<LFS>::Type couplingLocalFunctionSpace(const LFS& lfs)
  {
    return lfs.template getChild<Traits::couplingLFSIndex>();
  }

  const CouplingOperator& couplingOperator() const {
    return _operator;
  }

  const LocalSubProblem& localSubProblem() const {
    return _localSubProblem;
  }

  const RemoteSubProblem& remoteSubProblem() const {
    return _remoteSubProblem;
  }

  template<typename SDS1, typename SDS2>
  bool appliesTo(const SDS1& localSubDomains, const SDS2& remoteSubDomains) const {
    return localSubProblem().appliesTo(localSubDomains) && remoteSubProblem().appliesTo(remoteSubDomains);
  }

private:
  const LocalSubProblem& _localSubProblem;
  const RemoteSubProblem& _remoteSubProblem;

protected:
  CouplingOperator& _operator;
};


template<
  typename LocalSubProblem,
  typename RemoteSubProblem,
  typename CouplingOperator,
  std::size_t couplingLFSIndex
  >
struct is_coupling<EnrichedCoupling<LocalSubProblem,RemoteSubProblem,CouplingOperator,couplingLFSIndex> >
{
  static const bool value = true;
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_COUPLING_HH
