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
{

  typedef TypeTree::TransformTree<CouplingGridFunctionSpace,gfs_to_ordering<CouplingGridFunctionSpace> > ordering_transformation;


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
  CouplingGridFunctionSpace (const GV& gridview, const FEM& fem, const Predicate& predicate, const CE& ce_, const typename Traits::Backend& backend = typename Traits::Backend())
    : defaultce(ce_)
    , gv(gridview)
    , pfem(stackobject_to_shared_ptr(fem))
    , predicate_(predicate)
    , ce(ce_)
    , _backend(backend)
  {
  }

  //! constructor
  CouplingGridFunctionSpace (const GV& gridview, const FEM& fem, const Predicate& predicate, const typename Traits::Backend& backend = typename Traits::Backend())
    : gv(gridview)
    , pfem(stackobject_to_shared_ptr(fem))
    , ce(defaultce)
    , predicate_(predicate)
    , _backend(backend)
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

  //------------------------------
  // generic data handle interface
  //------------------------------

  //! returns true if data for this codim should be communicated
  bool dataHandleContains (int dim, int codim) const
  {
    return (codimUsed.find(codim)!=codimUsed.end());
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool dataHandleFixedSize (int dim, int codim) const
  {
    return false;
  }

  /*! how many objects of type DataType have to be sent for a given entity

    Note: Only the sender side needs to know this size.
  */
  template<class EntityType>
  size_t dataHandleSize (const EntityType& e) const
  {
    Dune::GeometryType gt=e.type();
    typename GV::IndexSet::IndexType index = gtoffset.find(gt)->second + gv.indexSet().index(e);
    return offset[index+1]-offset[index];
  }

  //! return vector of global indices associated with the given entity
  template<class EntityType>
  void dataHandleGlobalIndices (const EntityType& e,
                                std::vector<typename Traits::SizeType>& global) const
  {
    Dune::GeometryType gt=e.type();
    typename GV::IndexSet::IndexType index = gtoffset.find(gt)->second + gv.indexSet().index(e);
    unsigned int n = offset[index+1]-offset[index];
    global.resize(n);
    for (unsigned i=0; i<n; i++)
      global[i] = offset[index]+i;
  }

  //------------------------------


  /*
  // update information, e.g. when grid has changed
  void update ()
  {
    Dune::dinfo << "CouplingGridFunctionSpace:" << std::endl;

    // determine which geometry types are used
    // needs one traversal of the grid
    typedef std::set<Dune::GeometryType> GtUsedSetType;
    GtUsedSetType gtused;
    codimUsed.clear();
    for (ElementIterator it = gv.template begin<0>();
         it!=gv.template end<0>(); ++it)
      {
        for (IntersectionIterator iit = gv.ibegin(*it); iit != gv.iend(*it); ++iit)
          {
            if (!predicate_(*iit))
              continue;

            //EntityWrapper ew = buildWrapper(*iit);

            // check geometry type
            if ((pfem->find(*iit)).type()!=iit->type())
              DUNE_THROW(Exception, "geometry type mismatch in GridFunctionSpace");

            // get local coefficients for this entity
            const typename Traits::FiniteElementType::Traits::LocalCoefficientsType&
              lc = (pfem->find(*iit)).localCoefficients();

            // insert geometry type of all subentities into set
            for (std::size_t i=0; i<lc.size(); ++i)
              {
                Dune::GeometryType gt=Dune::ReferenceElements<double,GV::Grid::dimension - 1>
                  ::general(iit->type()).type(lc.localKey(i).subEntity(),lc.localKey(i).codim());
                gtused.insert(gt);
                codimUsed.insert(GV::Grid::dimension-gt.dim());
              }
          }
      }

    // now we can allocate one number per entity that holds degrees of freedom
    typename Traits::SizeType nentities = 0;
    gtoffset.clear();
    const typename GV::IndexSet& is=gv.indexSet();
    for (typename GtUsedSetType::iterator i=gtused.begin(); i!=gtused.end(); ++i)
      {
        gtoffset[*i] = nentities;
        Dune::dinfo << *i << ": " << is.size(*i)
                    << " entries in offset vector at " << nentities
                    << std::endl;
        nentities += is.size(*i);
      }
    nentities++; // add one additional dummy entry; this allows to compute size of last entity.
    offset.resize(nentities);
    std::fill(offset.begin(),offset.end(),0);
    Dune::dinfo << "allocated offset vector with size " << offset.size()
                << std::endl;

    // now compute the number of entries for each entity
    // requires second grid traversal
    nlocal = 0;
    for (ElementIterator it = gv.template begin<0>();
         it!=gv.template end<0>(); ++it)
      {

        for (IntersectionIterator iit = gv.ibegin(*it); iit != gv.iend(*it); ++iit)
          {

            if (!predicate_(*iit))
              continue;

            //EntityWrapper ew = buildWrapper(*iit);

            // get local coefficients for this entity
            const typename Traits::FiniteElementType::Traits::LocalCoefficientsType&
              lc = (pfem->find(*iit)).localCoefficients();

            // compute maximum number of degrees of freedom per element
            nlocal = std::max(nlocal,static_cast<typename Traits::SizeType>(lc.size()));

            DOFMapper<GV> dm(*iit);

            // compute maximum size for each subentity
            for (std::size_t i=0; i<lc.size(); ++i)
              {
                Dune::GeometryType gt=Dune::ReferenceElements<double,GV::Grid::dimension - 1>
                  ::general(iit->type()).type(lc.localKey(i).subEntity(),lc.localKey(i).codim());
                unsigned int index = gtoffset[gt] +
                  is.subIndex(dm.element(),
                              dm.mapSubIndex(lc.localKey(i).subEntity(),lc.localKey(i).codim()),
                              lc.localKey(i).codim()+1);
                offset[index] = std::max(offset[index],
                                         typename Traits::SizeType(lc.localKey(i).index()+1));
              }
          }
      }

    // now count global number of dofs and compute offset
    nglobal = 0;
    for (typename std::vector<typename Traits::SizeType>::iterator i=offset.begin();
         i!=offset.end(); ++i)
      {
        typename Traits::SizeType size = *i;
        *i = nglobal;
        nglobal += size;
      }
    Dune::dinfo << "total number of dofs is " << nglobal << std::endl;
  }
  */

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

  B& backend()
  {
    return _backend;
  }

  const B& backend() const
  {
    return _backend;
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
  B _backend;
  std::string _name;

  std::map<Dune::GeometryType,typename Traits::SizeType> gtoffset; // offset in vector for given geometry type
  std::vector<typename Traits::SizeType> offset; // offset into big vector for each entity;
  std::set<unsigned int> codimUsed;
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_COUPLINGGRIDFUNCTIONSPACE_HH
