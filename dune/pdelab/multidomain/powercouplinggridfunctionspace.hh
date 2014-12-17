// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTIDOMAIN_POWERCOUPLINGGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_MULTIDOMAIN_POWERCOUPLINGGRIDFUNCTIONSPACE_HH

#include <dune/pdelab/multidomain/couplinglocalfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {

template<typename T, std::size_t k,
         typename Backend,
         typename OrderingTag = LexicographicOrderingTag>
class PowerCouplingGridFunctionSpace
  : public TypeTree::PowerNode<T,k>
  , public PowerCompositeGridFunctionSpaceBase<
      PowerCouplingGridFunctionSpace<T, k, Backend, OrderingTag>,
      typename T::Traits::GridViewType,
      Backend,
      OrderingTag,
      k>
  , public DataHandleProvider<PowerCouplingGridFunctionSpace<T, k, Backend, OrderingTag> >
{
  typedef TypeTree::PowerNode<T,k> BaseT;

  typedef PowerCompositeGridFunctionSpaceBase<
    PowerCouplingGridFunctionSpace,
    typename T::Traits::GridViewType,
    Backend,
    OrderingTag,
    k
    > ImplementationBase;

  friend class PowerCompositeGridFunctionSpaceBase<
    PowerCouplingGridFunctionSpace,
    typename T::Traits::GridViewType,
    Backend,
    OrderingTag,
    k>;

  template<typename,typename>
  friend class GridFunctionSpaceBase;

  typedef TypeTree::TransformTree<
    PowerCouplingGridFunctionSpace,
    gfs_to_ordering<PowerCouplingGridFunctionSpace>
    > ordering_transformation;

public:
  typedef PowerCouplingGridFunctionSpaceTag ImplementationTag;

  typedef typename ordering_transformation::type Ordering;

  //! export traits class
  typedef typename ImplementationBase::Traits Traits;
  typedef typename BaseT::ChildType::Traits::GridView::template Codim<0>::Entity Entity;
  typedef typename BaseT::ChildType::Intersection Intersection;

  PowerCouplingGridFunctionSpace(T& c,
                                 const Backend& backend = Backend(),
                                 const OrderingTag ordering_tag = OrderingTag())
    : BaseT(c)
    , ImplementationBase(backend,ordering_tag)
  {}

  PowerCouplingGridFunctionSpace(T& c0,
                                  T& c1,
                                  const Backend& backend = Backend(),
                                  const OrderingTag ordering_tag = OrderingTag())
    : BaseT(c0,c1)
    , ImplementationBase(backend,ordering_tag)
  {}

  PowerCouplingGridFunctionSpace(T& c0,
                                 T& c1,
                                 T& c2,
                                 const Backend& backend = Backend(),
                                 const OrderingTag ordering_tag = OrderingTag())
    : BaseT(c0,c1,c2)
    , ImplementationBase(backend,ordering_tag)
  {}

  PowerCouplingGridFunctionSpace(T& c0,
                                 T& c1,
                                 T& c2,
                                 T& c3,
                                 const Backend& backend = Backend(),
                                 const OrderingTag ordering_tag = OrderingTag())
    : BaseT(c0,c1,c2,c3)
    , ImplementationBase(backend,ordering_tag)
  {}

  PowerCouplingGridFunctionSpace(T& c0,
                                 T& c1,
                                 T& c2,
                                 T& c3,
                                 T& c4,
                                 const Backend& backend = Backend(),
                                 const OrderingTag ordering_tag = OrderingTag())
    : BaseT(c0,c1,c2,c3,c4)
    , ImplementationBase(backend,ordering_tag)
  {}

  PowerCouplingGridFunctionSpace(T& c0,
                                 T& c1,
                                 T& c2,
                                 T& c3,
                                 T& c4,
                                 T& c5,
                                 const Backend& backend = Backend(),
                                 const OrderingTag ordering_tag = OrderingTag())
    : BaseT(c0,c1,c2,c3,c4,c5)
    , ImplementationBase(backend,ordering_tag)
  {}

  PowerCouplingGridFunctionSpace(T& c0,
                                 T& c1,
                                 T& c2,
                                 T& c3,
                                 T& c4,
                                 T& c5,
                                 T& c6,
                                 const Backend& backend = Backend(),
                                 const OrderingTag ordering_tag = OrderingTag())
    : BaseT(c0,c1,c2,c3,c4,c5,c6)
    , ImplementationBase(backend,ordering_tag)
  {}

  PowerCouplingGridFunctionSpace(T& c0,
                                 T& c1,
                                 T& c2,
                                 T& c3,
                                 T& c4,
                                 T& c5,
                                 T& c6,
                                 T& c7,
                                 const Backend& backend = Backend(),
                                 const OrderingTag ordering_tag = OrderingTag())
    : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
    , ImplementationBase(backend,ordering_tag)
  {}

  PowerCouplingGridFunctionSpace(T& c0,
                                 T& c1,
                                 T& c2,
                                 T& c3,
                                 T& c4,
                                 T& c5,
                                 T& c6,
                                 T& c7,
                                 T& c8,
                                 const Backend& backend = Backend(),
                                 const OrderingTag ordering_tag = OrderingTag())
    : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
    , ImplementationBase(backend,ordering_tag)
  {}

  PowerCouplingGridFunctionSpace(T& c0,
                                 T& c1,
                                 T& c2,
                                 T& c3,
                                 T& c4,
                                 T& c5,
                                 T& c6,
                                 T& c7,
                                 T& c8,
                                 T& c9,
                                 const Backend& backend = Backend(),
                                 const OrderingTag ordering_tag = OrderingTag())
    : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
    , ImplementationBase(backend,ordering_tag)
  {}

  bool contains(const Entity& e) const
  {
    return false;
  }

  bool contains(const Intersection& is) const
  {
    return this->child(0).contains(is);
  }

      //! Direct access to the DOF ordering.
      const Ordering &ordering() const
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return *_ordering;
      }

      //! Direct access to the DOF ordering.
      Ordering &ordering()
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return *_ordering;
      }

      //! Direct access to the storage of the DOF ordering.
      std::shared_ptr<const Ordering> orderingStorage() const
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return _ordering;
      }

      //! Direct access to the storage of the DOF ordering.
      std::shared_ptr<Ordering> orderingStorage()
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return _ordering;
      }

    private:

      // This method here is to avoid a double update of the Ordering when the user calls
      // GFS::update() before GFS::ordering().
      void create_ordering() const
      {
        _ordering = std::make_shared<Ordering>(ordering_transformation::transform(*this));
      }

      mutable std::shared_ptr<Ordering> _ordering;

};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_POWERCOUPLINGGRIDFUNCTIONSPACE_HH
