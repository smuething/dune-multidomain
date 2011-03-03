// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_ENTITYWRAPPERS_HH
#define DUNE_PDELAB_MULTIDOMAIN_ENTITYWRAPPERS_HH

#include <dune/pdelab/common/geometrywrapper.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


template<typename GV>
class ElementWrapper
  : public ElementGeometry<typename GV::template Codim<0>::Entity>
{

  typedef ElementGeometry<typename GV::template Codim<0>::Entity> BaseT;

public:

  typedef typename GV::template Codim<0>::Entity Entity;
  typedef typename GV::Grid::MDGridTraits::template Codim<0>::SubDomainSet SubDomainSet;

  ElementWrapper(const Entity& e, const SubDomainSet& sds)
    : BaseT(e)
    , _subDomains(sds)
  {}

  const SubDomainSet& subDomains() const
  {
    return _subDomains;
  }

private:

  const SubDomainSet& _subDomains;

};

template<typename GV>
class SkeletonIntersectionWrapper
  : public IntersectionGeometry<typename GV::Intersection>
{

  typedef IntersectionGeometry<typename GV::Intersection> BaseT;
  typedef typename GV::template Codim<0>::EntityPointer EntityPointer;

public:

  typedef typename GV::Intersection Intersection;
  typedef ElementWrapper<GV> EntityWrapper;

  SkeletonIntersectionWrapper(const Intersection& intersection,
                              std::size_t intersectionIndex,
                              const EntityWrapper& insideElement,
                              const typename EntityWrapper::SubDomainSet& outsideSubDomains)
    : BaseT(intersection,intersectionIndex)
    , _inside(insideElement)
    , _outsideEntityPointer(intersection.outside())
    , _outside(*_outsideEntityPointer,outsideSubDomains)
  {}

  const EntityWrapper& insideElement() const
  {
    return _inside;
  }

  const EntityWrapper& outsideElement() const
  {
    return _outside;
  }

private:

  const EntityWrapper& _inside;
  const EntityPointer _outsideEntityPointer;
  const EntityWrapper _outside;

};


template<typename GV>
class BoundaryIntersectionWrapper
  : public IntersectionGeometry<typename GV::Intersection>
{

  typedef IntersectionGeometry<typename GV::Intersection> BaseT;

public:

  typedef typename GV::Intersection Intersection;
  typedef ElementWrapper<GV> EntityWrapper;

  BoundaryIntersectionWrapper(const Intersection& intersection,
                              std::size_t intersectionIndex,
                              const EntityWrapper& insideElement)
    : BaseT(intersection,intersectionIndex)
    , _inside(insideElement)
  {}

  const EntityWrapper& insideElement() const
  {
    return _inside;
  }

private:

  const EntityWrapper& _inside;

};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_ENTITYWRAPPERS_HH
