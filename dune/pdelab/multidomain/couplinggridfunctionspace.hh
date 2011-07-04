#ifndef DUNE_MULTIDOMAIN_COUPLINGGRIDFUNCTIONSPACE_HH
#define DUNE_MULTIDOMAIN_COUPLINGGRIDFUNCTIONSPACE_HH

#include <dune/pdelab/multidomain/couplinglocalfunctionspace.hh>
#include <dune/pdelab/common/typetree.hh>

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

template<typename GV>
class DOFMapper
{

public:
  typedef typename GV::template Codim<0>::EntityPointer ElementPointer;
  typedef typename GV::template Codim<0>::Entity Element;
  typedef typename GV::Intersection Intersection;
  typedef typename Element::Geometry::LocalCoordinate Coordinate;

  DOFMapper(const Intersection& is)
    : intersection_(is)
    , inside_(DOFStorageOnInside(is))
    , elementPointer_(inside_ ? is.inside() : is.outside())
  {}

  const Element& element() const
  {
    return *elementPointer_;
  }

  int mapSubIndex(int i, int codim) const
  {
    if (codim == 0)
      return inside_ ? intersection_.indexInInside() : intersection_.indexInOutside();
    if (codim == GV::Grid::dimension - 1)
      {
        Coordinate c = (inside_ ? intersection_.geometryInInside() : intersection_.geometryInOutside()).corner(i);
        const VertexList& vl = vertexList(element().type());
        typename VertexList::const_iterator it = std::find_if(vl.begin(),vl.end(),[c](Coordinate rc) { rc -= c; return rc.two_norm() < 1e-10;}); // TODO: don't use lambda expression!
        assert(it != vl.end());
        return it - vl.begin();
      }
    DUNE_THROW(GridError,"currently, only degrees of freedom on the intersection and the vertices are supported!");
  }

private:
  const Intersection& intersection_;
  const bool inside_;
  ElementPointer elementPointer_;

  typedef std::vector<Coordinate> VertexList;
  typedef std::map<Dune::GeometryType,VertexList> VertexListMap;

  static VertexListMap vertexListMap_;

  static bool DOFStorageOnInside(const Intersection& is)
  {
    if (is.conforming())
      return true;
    const int insideLevel = is.inside().level();
    const int outsideLevel = is.outside().level();
    assert(insideLevel != outsideLevel);
    return insideLevel > outsideLevel;
  }

  static const VertexList& vertexList(Dune::GeometryType gt)
  {
    typename VertexListMap::iterator it = vertexListMap_.find(gt);
    if (it != vertexListMap_.end())
      return it->second;
    vertexListMap_.insert({gt,buildVertexList(gt)});
    return vertexListMap_[gt];
  }

  static VertexList buildVertexList(Dune::GeometryType gt)
  {
    const Dune::GenericReferenceElement<typename GV::Grid::ctype,GV::Grid::dimension>& refEl =
      Dune::GenericReferenceElements<typename GV::Grid::ctype,GV::Grid::dimension>::general(gt);
    VertexList vl(refEl.size(GV::Grid::dimension));
    for (std::size_t i = 0; i < vl.size(); ++i)
      vl[i] = refEl.position(i,GV::Grid::dimension);
    return std::move(vl);
  }

};

template<typename GV>
typename DOFMapper<GV>::VertexListMap DOFMapper<GV>::vertexListMap_;

template<typename GV, typename LFEM, typename Predicate_, typename CE=NoConstraints,
         typename B=StdVectorBackend>
class CouplingGridFunctionSpace
  : public Dune::PDELab::TypeTree::LeafNode
{
public:
  //! export Traits class
  typedef GridFunctionSpaceTraits<GV,LFEM,CE,B> Traits;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::Intersection Intersection;
  typedef typename GV::IntersectionIterator IntersectionIterator;

  typedef Predicate_ Predicate;

  typedef LeafOrdering<CouplingGridFunctionSpace> Ordering;

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
    typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;
  private:
    ConstraintsContainer () {}
  };

  //! define local function space parametrized by self
  typedef CouplingLocalFunctionSpaceNode<CouplingGridFunctionSpace> LocalFunctionSpace;

  //! constructor
  CouplingGridFunctionSpace (const GV& gridview, const LFEM& lfem, const Predicate& predicate, const CE& ce_)
    : defaultce(ce_)
    , gv(gridview)
    , plfem(stackobject_to_shared_ptr(lfem))
    , predicate_(predicate)
    , ordering_(make_shared<Ordering>(*this))
    , ce(ce_)
  {
    update();
  }

  //! constructor
  CouplingGridFunctionSpace (const GV& gridview, const LFEM& lfem, const Predicate& predicate)
    : gv(gridview)
    , plfem(stackobject_to_shared_ptr(lfem))
    , ce(defaultce)
    , predicate_(predicate)
    , ordering_(make_shared<Ordering>(*this))
  {
    update();
  }

  //! get grid view
  const GV& gridview () const
  {
    return gv;
  }

  // get finite element map, I think we dont need it
  const LFEM& localFiniteElementMap () const
  {
    return *plfem;
  }

  //! get dimension of root finite element space
  typename Traits::SizeType globalSize () const
  {
    return nglobal;
  }

  //! get dimension of this finite element space
  typename Traits::SizeType size () const
  {
    return nglobal;
  }

  //! get max dimension of shape function space
  //! \todo What are the exact semantics of maxLocalSize?
  typename Traits::SizeType maxLocalSize () const
  {
    return nlocal;
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

  // return constraints engine
  const typename Traits::ConstraintsType& constraints () const
  {
    return ce;
  }

  //! compute global indices for one element
  void globalIndices (const typename Traits::FiniteElementType& lfe,
                      const Intersection& is,
                      std::vector<typename Traits::SizeType>& global) const
  {
    // get layout of entity
    if (!predicate_(is))
      {
        global.resize(0);
        return;
      }
    const typename Traits::FiniteElementType::Traits::LocalCoefficientsType&
      lc = lfe.localCoefficients();
    global.resize(lc.size());

    DOFMapper<GV> dm(is);

    for (std::size_t i=0; i<lc.size(); ++i)
      {
        // get geometry type of subentity
        Dune::GeometryType gt=Dune::GenericReferenceElements<double,GV::Grid::dimension - 1>
          ::general(lfe.type()).type(lc.localKey(i).subEntity(),lc.localKey(i).codim());

        // evaluate consecutive index of subentity
        int index = gv.indexSet().subIndex(dm.element(),
                                           dm.mapSubIndex(lc.localKey(i).subEntity(),lc.localKey(i).codim()),
                                           lc.localKey(i).codim() + 1);

        // now compute
        global[i] = offset[(gtoffset.find(gt)->second)+index]+lc.localKey(i).index();
      }
  }

  // global Indices from element, needs additional finite element lookup
  void globalIndices (const Intersection& is,
                      std::vector<typename Traits::SizeType>& global) const
  {
    globalIndices(plfem->find(is),is,global);
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
            if ((plfem->find(*iit)).type()!=iit->type())
              DUNE_THROW(Exception, "geometry type mismatch in GridFunctionSpace");

            // get local coefficients for this entity
            const typename Traits::FiniteElementType::Traits::LocalCoefficientsType&
              lc = (plfem->find(*iit)).localCoefficients();

            // insert geometry type of all subentities into set
            for (std::size_t i=0; i<lc.size(); ++i)
              {
                Dune::GeometryType gt=Dune::GenericReferenceElements<double,GV::Grid::dimension - 1>
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
              lc = (plfem->find(*iit)).localCoefficients();

            // compute maximum number of degrees of freedom per element
            nlocal = std::max(nlocal,static_cast<typename Traits::SizeType>(lc.size()));

            DOFMapper<GV> dm(*iit);

            // compute maximum size for each subentity
            for (std::size_t i=0; i<lc.size(); ++i)
              {
                Dune::GeometryType gt=Dune::GenericReferenceElements<double,GV::Grid::dimension - 1>
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
    ordering_->update();
  }

  //! Direct access to the DOF ordering.
  const Ordering &ordering() const { return *ordering_; }

  //! Direct access to the storage of the DOF ordering.
  shared_ptr<const Ordering> orderingPtr() const { return ordering_; }


private:
  CE defaultce;
  const GV& gv;
  shared_ptr<LFEM const> plfem;
  typename Traits::SizeType nlocal;
  typename Traits::SizeType nglobal;
  const CE& ce;
  const Predicate_& predicate_;
  shared_ptr<Ordering> ordering_;

  std::map<Dune::GeometryType,typename Traits::SizeType> gtoffset; // offset in vector for given geometry type
  std::vector<typename Traits::SizeType> offset; // offset into big vector for each entity;
  std::set<unsigned int> codimUsed;
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_COUPLINGGRIDFUNCTIONSPACE_HH
