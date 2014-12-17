#ifndef DUNE_MULTIDOMAIN_COUPLINGGRIDFUNCTIONSPACE_HH
#define DUNE_MULTIDOMAIN_COUPLINGGRIDFUNCTIONSPACE_HH

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/multidomain/couplinglocalfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinggfsordering.hh>
#include <dune/pdelab/multidomain/dofmapper.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {


template<typename GV, typename LeftSubDomainPredicate, typename RightSubDomainPredicate = LeftSubDomainPredicate>
class SubProblemSubProblemInterface
{

public:

  typedef typename GV::Intersection Intersection;

  SubProblemSubProblemInterface(const GV& gv, const LeftSubDomainPredicate& lsp, const RightSubDomainPredicate& rsp)
    : gv_(gv)
    , leftSubDomainPredicate_(lsp)
    , rightSubDomainPredicate_(rsp)
  {}

  bool operator()(const Intersection& is) const
  {
    if (!is.neighbor())
      return false;
    return leftSubDomainPredicate_(gv_.indexSet().subDomains(*(is.inside())))
      && rightSubDomainPredicate_(gv_.indexSet().subDomains(*(is.outside())));
  }

private:
  GV gv_;
  const LeftSubDomainPredicate& leftSubDomainPredicate_;
  const RightSubDomainPredicate& rightSubDomainPredicate_;

};


template<typename GV, typename FEM, typename Predicate_, typename CE=NoConstraints,
         typename B=ISTLVectorBackend<>, typename O=DefaultLeafOrderingTag>
class CouplingGridFunctionSpace
  : public TypeTree::LeafNode
  , public GridFunctionSpaceBase<
             GridFunctionSpace<GV,FEM,CE,B,O>,
             GridFunctionSpaceTraits<GV,FEM,CE,B,O>
             >
  , public GridFunctionOutputParameters
  , public DataHandleProvider<GridFunctionSpace<GV,FEM,CE,B,O> >

{

  typedef TypeTree::TransformTree<CouplingGridFunctionSpace,gfs_to_ordering<CouplingGridFunctionSpace> > ordering_transformation;


  typedef GridFunctionSpaceBase<
    GridFunctionSpace<GV,FEM,CE,B,O>,
    GridFunctionSpaceTraits<GV,FEM,CE,B,O>
    > BaseT;

public:
  //! export Traits class
  typedef GridFunctionSpaceTraits<GV,FEM,CE,B,O> Traits;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::Intersection Intersection;
  typedef typename GV::IntersectionIterator IntersectionIterator;

  typedef Predicate_ Predicate;

  typedef CouplingGridFunctionSpaceTag ImplementationTag;

  typedef O OrderingTag;

  typedef typename ordering_transformation::Type Ordering;

  //! extract type of container storing Ts
  template<typename T>
  struct VectorContainer
  {
    //! \brief define Type as the Type of a container of E's
    typedef typename B::template VectorContainer<CouplingGridFunctionSpace,T> Type;
  private:
    VectorContainer () {}
  };

  //! extract type for storing constraints
  template<typename E>
  struct ConstraintsContainer
  {
    //! \brief define Type as the Type of a container of E's
    typedef ConstraintsTransformation<
      typename Ordering::Traits::DOFIndex,
      typename Ordering::Traits::ContainerIndex,
      E> Type;

  private:
    ConstraintsContainer () {}
  };

  //! constructor
  CouplingGridFunctionSpace (const GV& gridview, const FEM& fem, const Predicate& predicate, const CE& ce_, const typename Traits::Backend& backend = typename Traits::Backend(), const typename Traits::OrderingTag& ordering_tag = typename Traits::OrderingTag())
    : BaseT(backend,ordering_tag)
    , defaultce(ce_)
    , gv(gridview)
    , pfem(stackobject_to_shared_ptr(fem))
    , predicate_(predicate)
    , ce(ce_)
  {
  }

  //! constructor
  CouplingGridFunctionSpace (const GV& gridview, const FEM& fem, const Predicate& predicate, const typename Traits::Backend& backend = typename Traits::Backend(), const typename Traits::OrderingTag& ordering_tag = typename Traits::OrderingTag())
    : BaseT(backend,ordering_tag)
    , gv(gridview)
    , pfem(stackobject_to_shared_ptr(fem))
    , ce(defaultce)
    , predicate_(predicate)
  {
  }

  //! get grid view
  const GV& gridview () const DUNE_DEPRECATED_MSG("Use gridView() instead of gridview()")
  {
    return gv;
  }


  //! get grid view
  const GV& gridView () const
  {
    return gv;
  }

  // get finite element map, I think we dont need it
  const FEM& finiteElementMap () const
  {
    return *pfem;
  }

  //! get dimension of root finite element space
  typename Traits::SizeType globalSize () const
  {
    return ordering().size();
  }

  //! get dimension of this finite element space
  typename Traits::SizeType size () const
  {
    return ordering().size();
  }

  //! get max dimension of shape function space
  //! \todo What are the exact semantics of maxLocalSize?
  typename Traits::SizeType maxLocalSize () const
  {
    return ordering().maxLocalSize();
  }

  //! map index from our index set [0,size()-1] to root index set
  typename Traits::SizeType upMap (typename Traits::SizeType i) const
  {
    return i;
  }

  bool contains(const Intersection& is) const
  {
    return predicate_(is);
  }

  const Predicate& predicate() const
  {
    return predicate_;
  }

  // return constraints engine
  const typename Traits::ConstraintsType& constraints () const
  {
    return ce;
  }


  //! get finite element map
  shared_ptr<const FEM> finiteElementMapStorage () const
  {
    return pfem;
  }

  //! Direct access to the DOF ordering.
  const Ordering &ordering() const
  {
    return *orderingStorage();
  }

  //! Direct access to the DOF ordering.
  Ordering &ordering()
  {
    return *orderingStorage();
  }

  //! Direct access to the storage of the DOF ordering.
  shared_ptr<const Ordering> orderingStorage() const
  {
    if (!_ordering)
      {
        _ordering = make_shared<Ordering>(ordering_transformation::transform(*this));
        _ordering->update();
      }
    return _ordering;
  }

  //! Direct access to the storage of the DOF ordering.
  shared_ptr<Ordering> orderingStorage()
  {
    if (!_ordering)
      {
        _ordering = make_shared<Ordering>(ordering_transformation::transform(*this));
        _ordering->update();
      }
    return _ordering;
  }

  const std::string& name() const
  {
    return _name;
  }

  void name(const std::string& name)
  {
    _name = name;
  }

private:
  CE defaultce;
  const GV& gv;
  shared_ptr<FEM const> pfem;
  typename Traits::SizeType nlocal;
  typename Traits::SizeType nglobal;
  const CE& ce;
  const Predicate_& predicate_;
  mutable shared_ptr<Ordering> _ordering;
  std::string _name;

};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_COUPLINGGRIDFUNCTIONSPACE_HH
