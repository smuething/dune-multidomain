#ifndef DUNE_MULTIDOMAIN_SUBDOMAINSET_HH
#define DUNE_MULTIDOMAIN_SUBDOMAINSET_HH

namespace Dune {

namespace PDELab {

namespace MultiDomain {

template<typename SubDomainSet>
void initSubDomainSet(SubDomainSet& subDomainSet)
{
}

template<typename SubDomainSet, typename... T>
void initSubDomainSet(SubDomainSet& subDomainSet, typename SubDomainSet::DomainType subDomain,T... subDomains)
{
  subDomainSet.add(subDomain);
  initSubDomainSet(subDomainSet,subDomains...);
}

template<typename SubDomainSet>
class SubDomainSetHolder
{

public:

  typedef SubDomainSet SubDomainSetType;
  typedef typename SubDomainSet::DomainType SubDomainType;

  template<typename... T>
  SubDomainSetHolder(T... subDomains)
  {
    initSubDomainSet(_subDomainSet,subDomains...);
  }

  const SubDomainSet& subDomainSet() const
  {
    return _subDomainSet;
  }

private:
  SubDomainSet _subDomainSet;
};


template<typename SubDomainSet>
class IncludesSubDomains : public SubDomainSetHolder<SubDomainSet>
{

  typedef SubDomainSetHolder<SubDomainSet> BaseT;

public:

  using BaseT::SubDomainSetType;
  using BaseT::SubDomainType;

  template<typename... T>
  IncludesSubDomains(T... subDomains) :
    BaseT(subDomains...)
  {
  }

  bool operator()(const SubDomainSet& subDomainSet) const {
    return subDomainSet.containsAll(this->subDomainSet());
  }

};


template<typename SubDomainSet>
class EqualsSubDomains : public SubDomainSetHolder<SubDomainSet>
{

  typedef SubDomainSetHolder<SubDomainSet> BaseT;

public:

  using BaseT::SubDomainSetType;
  using BaseT::SubDomainType;

  template<typename... T>
  EqualsSubDomains(T... subDomains) :
    BaseT(subDomains...)
  {
  }

  bool operator()(const SubDomainSet& subDomainSet) const {
    return subDomainSet == this->subDomainSet();
  }

};


} // namespace MultiDomain

} // namespace PDELab

} // namespace Dune

#endif // DUNE_MULTIDOMAIN_SUBDOMAINSET_HH
