// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTIDOMAIN_POWERCOUPLINGGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_MULTIDOMAIN_POWERCOUPLINGGRIDFUNCTIONSPACE_HH

#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {

template<typename T, std::size_t k,
         typename OrderingTag = LexicographicOrderingTag>
class PowerCouplingGridFunctionSpace
  : public TypeTree::PowerNode<T,k>
  , public PowerCompositeGridFunctionSpaceBase<PowerCouplingGridFunctionSpace<T, k, OrderingTag>,
                                               typename T::Traits::GridViewType,
                                               typename T::Traits::BackendType,
                                               OrderingTag,
                                               k
                                               >
{
  typedef TypeTree::PowerNode<T,k> BaseT;

  typedef PowerCompositeGridFunctionSpaceBase<
    PowerCouplingGridFunctionSpace,
    typename T::Traits::GridViewType,
    typename T::Traits::BackendType,
    OrderingTag,
    k
    > ImplementationBase;

  friend class PowerCompositeGridFunctionSpaceBase<
    PowerCouplingGridFunctionSpace,
    typename T::Traits::GridViewType,
    typename T::Traits::BackendType,
    OrderingTag,
    k>;

public:
  typedef PowerCouplingGridFunctionSpaceTag ImplementationTag;

  typedef typename TransformPowerGFSToOrdering<OrderingTag>::
  template result<
    typename ImplementationBase::Traits,
    const typename T::Ordering,
    k
    >::type Ordering;

  //! export traits class
  typedef typename ImplementationBase::Traits Traits;
  typedef typename BaseT::ChildType::Intersection Intersection;

  PowerCouplingGridFunctionSpace(T& c)
  : BaseT(c)
  {
    initOrdering();
  }

  PowerCouplingGridFunctionSpace (T& c0,
                                  T& c1)
    : BaseT(c0,c1)
  {
    initOrdering();
  }

  PowerCouplingGridFunctionSpace (T& c0,
                                  T& c1,
                                  T& c2)
    : BaseT(c0,c1,c2)
  {
    initOrdering();
  }

  PowerCouplingGridFunctionSpace (T& c0,
                                  T& c1,
                                  T& c2,
                                  T& c3)
    : BaseT(c0,c1,c2,c3)
  {
    initOrdering();
  }

  PowerCouplingGridFunctionSpace (T& c0,
                                  T& c1,
                                  T& c2,
                                  T& c3,
                                  T& c4)
    : BaseT(c0,c1,c2,c3,c4)
  {
    initOrdering();
  }

  PowerCouplingGridFunctionSpace (T& c0,
                                  T& c1,
                                  T& c2,
                                  T& c3,
                                  T& c4,
                                  T& c5)
    : BaseT(c0,c1,c2,c3,c4,c5)
  {
    initOrdering();
  }

  PowerCouplingGridFunctionSpace (T& c0,
                                  T& c1,
                                  T& c2,
                                  T& c3,
                                  T& c4,
                                  T& c5,
                                  T& c6)
    : BaseT(c0,c1,c2,c3,c4,c5,c6)
  {
    initOrdering();
  }

  PowerCouplingGridFunctionSpace (T& c0,
                                  T& c1,
                                  T& c2,
                                  T& c3,
                                  T& c4,
                                  T& c5,
                                  T& c6,
                                  T& c7)
    : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
  {
    initOrdering();
  }

  PowerCouplingGridFunctionSpace (T& c0,
                                  T& c1,
                                  T& c2,
                                  T& c3,
                                  T& c4,
                                  T& c5,
                                  T& c6,
                                  T& c7,
                                  T& c8)
    : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
  {
    initOrdering();
  }

  PowerCouplingGridFunctionSpace (T& c0,
                                  T& c1,
                                  T& c2,
                                  T& c3,
                                  T& c4,
                                  T& c5,
                                  T& c6,
                                  T& c7,
                                  T& c8,
                                  T& c9)
    : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
  {
    initOrdering();
  }

  //! Direct access to the DOF ordering.
  const Ordering &ordering() const { return *orderingp; }

  //! Direct access to the storage of the DOF ordering.
  shared_ptr<const Ordering> orderingPtr() const { return orderingp; }


  bool contains(const Intersection& is) const
  {
    return this->child(0).contains(is);
  }

private:
  void initOrdering() {
    typename Ordering::NodeStorage transformedChildren;
    for(std::size_t childIndex = 0; childIndex < BaseT::CHILDREN;
        ++childIndex)
      transformedChildren[childIndex] =
        this->child(childIndex).orderingPtr();
    orderingp = make_shared<Ordering>(*this, transformedChildren);
  }

  shared_ptr<Ordering> orderingp;
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_POWERCOUPLINGGRIDFUNCTIONSPACE_HH
