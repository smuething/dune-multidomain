#ifndef DUNE_MULTIDOMAIN_SUBDOMAINSET_HH
#define DUNE_MULTIDOMAIN_SUBDOMAINSET_HH

#include <initializer_list>

#include <dune/common/deprecated.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {

template<typename SubDomainSet>
void initSubDomainSet(SubDomainSet& subDomainSet)
{
}

template<typename SubDomainSet, typename... T>
void initSubDomainSet(SubDomainSet& subDomainSet, typename SubDomainSet::SubDomainIndex subDomain,T... subDomains)
{
  subDomainSet.add(subDomain);
  initSubDomainSet(subDomainSet,subDomains...);
}

template<typename MultiDomainGrid>
class SubDomainSetHolder
{

public:

  typedef typename MultiDomainGrid::LeafGridView::IndexSet::SubDomainSet SubDomainSet;
  typedef typename SubDomainSet::SubDomainIndex SubDomainIndex;

  template<typename... T>
  SubDomainSetHolder(T... subDomains)
  {
    initSubDomainSet(_subDomainSet,subDomains...);
  }

  SubDomainSetHolder(std::initializer_list<SubDomainIndex> subdomains)
  {
    for (auto subdomain : subdomains)
      _subDomainSet.add(subdomain);
  }

  const SubDomainSet& subDomainSet() const
  {
    return _subDomainSet;
  }

private:
  SubDomainSet _subDomainSet;
};


template<typename MultiDomainGrid>
class
DUNE_DEPRECATED_MSG("This class has been renamed to SubDomainSubsetCondition")
SubDomainSupersetCondition : public SubDomainSetHolder<MultiDomainGrid>
{

  typedef SubDomainSetHolder<MultiDomainGrid> BaseT;
  typedef typename BaseT::SubDomainSet SubDomainSet;

public:

  template<typename... T>
  SubDomainSupersetCondition(T... subDomains) :
    BaseT(subDomains...)
  {
  }



  bool operator()(const SubDomainSet& subDomainSet) const {
    return subDomainSet.containsAll(this->subDomainSet());
  }

};

template<typename MultiDomainGrid>
class SubDomainSubsetCondition : public SubDomainSetHolder<MultiDomainGrid>
{

  typedef SubDomainSetHolder<MultiDomainGrid> BaseT;
  typedef typename BaseT::SubDomainSet SubDomainSet;
  typedef typename SubDomainSet::SubDomainIndex SubDomainIndex;

public:

  template<typename... T>
  SubDomainSubsetCondition(T... subDomains) :
    BaseT(subDomains...)
  {
  }

  SubDomainSubsetCondition(std::initializer_list<SubDomainIndex> subdomains)
    : BaseT(subdomains)
  {}

  bool operator()(const SubDomainSet& subDomainSet) const {
    return subDomainSet.containsAll(this->subDomainSet());
  }

};


template<typename MultiDomainGrid>
class SubDomainEqualityCondition : public SubDomainSetHolder<MultiDomainGrid>
{

  typedef SubDomainSetHolder<MultiDomainGrid> BaseT;
  typedef typename BaseT::SubDomainSet SubDomainSet;
  typedef typename SubDomainSet::SubDomainIndex SubDomainIndex;

public:

  template<typename... T>
  SubDomainEqualityCondition(T... subDomains) :
    BaseT(subDomains...)
  {
  }

  SubDomainEqualityCondition(std::initializer_list<SubDomainIndex> subdomains)
    : BaseT(subdomains)
  {}

  bool operator()(const SubDomainSet& subDomainSet) const {
    return subDomainSet == this->subDomainSet();
  }

};


} // namespace MultiDomain

} // namespace PDELab

} // namespace Dune

#endif // DUNE_MULTIDOMAIN_SUBDOMAINSET_HH
