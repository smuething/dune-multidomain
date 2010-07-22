#ifndef DUNE_MULTIDOMAIN_COUPLING_HH
#define DUNE_MULTIDOMAIN_COUPLING_HH

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

template<
  typename LocalSubProblem,
  typename RemoteSubProblem,
  typename CouplingOperator
  >
class Coupling
{

public:

  typedef CouplingTraits<LocalSubProblem,RemoteSubProblem,CouplingOperator> Traits;

  Coupling(const LocalSubProblem& localSubProblem,
           const RemoteSubProblem& remoteSubProblem,
           const CouplingOperator& couplingOperator)
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
  const CouplingOperator& _operator;
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_COUPLING_HH
