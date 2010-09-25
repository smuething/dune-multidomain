#ifndef DUNE_MULTIDOMAIN_COUPLINGLOCALFUNCTIONSPACE_HH
#define DUNE_MULTIDOMAIN_COUPLINGLOCALFUNCTIONSPACE_HH

namespace Dune {
namespace PDELab {
namespace MultiDomain {


template<typename GFS>
struct CouplingLocalFunctionSpaceTraits
{
  //! \brief the grid view where grid function is defined upon
  typedef GFS GridFunctionSpaceType;

  //! \brief Type to store indices from Backend
  typedef typename GFS::Traits::GridViewType GridViewType;

  //! \brief Type of codim 0 entity in the grid
  typedef typename GridViewType::Traits::template Codim<0>::Entity Element;

  typedef typename GridViewType::Traits::Intersection Intersection;

  //! \brief Type to store indices from Backend
  typedef typename GFS::Traits::SizeType SizeType;

  //! \brief Type of container to store indices
  typedef typename std::vector<SizeType> IndexContainer;

  //! \brief Type of local finite element
  typedef typename GFS::Traits::LocalFiniteElementType LocalFiniteElementType;

  //! \brief Type of constraints engine
  typedef typename GFS::Traits::ConstraintsType ConstraintsType;

};

template <typename GFS>
class CouplingLocalFunctionSpaceNode
  : public LeafNode
{
  typedef typename GFS::Traits::BackendType B;

  template<typename T, bool b, typename E, typename It, typename Int>
  friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
  template<typename T, typename E, typename It, typename Int, int n, int i>
  friend struct LocalFunctionSpaceBaseVisitChildMetaProgram;

public:
  typedef CouplingLocalFunctionSpaceTraits<GFS> Traits;

  //! \brief construct without associating a global function space
  CouplingLocalFunctionSpaceNode () {}

  //! \brief construct from global function space
  CouplingLocalFunctionSpaceNode (const GFS& gfs) :
    pgfs(&gfs), global(gfs.maxLocalSize())
  {
  }

  //! \brief initialize with grid function space
  void setup (const GFS& gfs)
  {
    pgfs = &gfs;
  }

  //! \brief get current size
  typename Traits::IndexContainer::size_type size () const
  {
    return n;
  }

  //! \brief get maximum possible size (which is maxLocalSize from grid function space)
  typename Traits::IndexContainer::size_type maxSize () const
  {
    return pgfs->maxLocalSize();
  }

  //! \brief get size of an appropriate local vector object
  /**
     this is the number of dofs of the complete local function
     space tree, i.e. the size() of the root node. The local
     vector objects must always have this size and the localIndex
     method maps into the range [0,localVectorSize()[
  */
  typename Traits::IndexContainer::size_type localVectorSize () const
  {
    return lvsize;
  }

  //! \brief map index in this local function space to root local function space
  typename Traits::IndexContainer::size_type localIndex (typename Traits::IndexContainer::size_type index) const
  {
    return offset+index;
  }

  //! \brief map index in this local function space to global index space
  typename Traits::SizeType globalIndex (typename Traits::IndexContainer::size_type index) const
  {
    return i[index];
  }

  //! \brief extract coefficients for one element from container
  template<typename GC, typename LC>
  void vread (const GC& globalcontainer, LC& localcontainer) const
  {
    localcontainer.resize(n);
    for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
      localcontainer[typename LC::size_type(k)] = B::access(globalcontainer,i[k]);
  }

  //! \brief write back coefficients for one element to container
  template<typename GC, typename LC>
  void vwrite (const LC& localcontainer, GC& globalcontainer) const
  {
    for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
      B::access(globalcontainer,i[k]) = localcontainer[typename LC::size_type(k)];
  }

  //! \brief add coefficients for one element to container
  template<typename GC, typename LC>
  void vadd (const LC& localcontainer, GC& globalcontainer) const
  {
    for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
      B::access(globalcontainer,i[k]) += localcontainer[typename LC::size_type(k)];
  }

  //! \brief print debug information about this local function space
  void debug () const
  {
    std::cout << n << " indices = (";
    for (typename Traits::IndexContainer::size_type k=0; k<n; k++)
      std::cout << i[k] << " ";
    std::cout << ")" << std::endl;
  }

  //! \brief bind local function space to entity
  /**

     This is a generic implementation of the bind function. It is
     parametrized with the NodeType, which the type of the derived
     LocalFunctionSpaceNode. Handing the NodeType as a parammeter
     avoid the need for the CRTP construct, but all derived
     classes have to add a method bind, which forward to this
     method.

     \param node reference to the derived node, the address must be the same as this
     \param e entity to bind to
  */
  template<typename NodeType>
  bool bind (NodeType& node, const typename Traits::Intersection& is)
  {
    // we should only call bind on out selfs
    assert(&node == this);

    if (!pgfs->contains(is))
      return false;

    // make offset
    typename Traits::IndexContainer::size_type offset=0;

    // compute sizes
    LocalFunctionSpaceBaseVisitNodeMetaProgram<NodeType,NodeType::isLeaf,
                                               typename Traits::Intersection,
                                               typename Traits::IndexContainer::iterator,
                                               typename Traits::IndexContainer::size_type>::
      reserve(node,is,offset);

    // now reserve space in vector
    global.resize(offset);
    lvsize = global.size();

    // initialize iterators and fill indices
    offset = 0;
    LocalFunctionSpaceBaseVisitNodeMetaProgram<NodeType,NodeType::isLeaf,
                                               typename Traits::Intersection,
                                               typename Traits::IndexContainer::iterator,
                                               typename Traits::IndexContainer::size_type>::
      fill_indices(node,is,global.begin(),offset,lvsize);

    // apply upMap
    assert(offset == global.size());
    for (typename Traits::IndexContainer::size_type i=0; i<offset; i++)
      global[i] = pgfs->upMap(global[i]);

    return true;
  }


  //! \brief get local finite element
  const typename Traits::LocalFiniteElementType& localFiniteElement () const
  {
    return *plfem;
  }

  //! \brief get constraints engine
  const typename Traits::ConstraintsType& constraints () const
  {
    return this->pgfs->constraints();
  }

private:
  CountingPointer<GFS const> pgfs;
  typename Traits::IndexContainer global;
  typename Traits::IndexContainer::iterator i;
  typename Traits::IndexContainer::size_type n;
  typename Traits::IndexContainer::size_type offset;
  typename Traits::IndexContainer::size_type lvsize;
  const typename Traits::LocalFiniteElementType* plfem;
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_COUPLINGLOCALFUNCTIONSPACE_HH
