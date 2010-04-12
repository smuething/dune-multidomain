// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_MULTIDOMAIN_MULTIDOMAINLOCALFUNCTIONSPACE_HH
#define DUNE_MULTIDOMAIN_MULTIDOMAINLOCALFUNCTIONSPACE_HH

#include <vector>
#include <dune/pdelab/multidomain/variadiccompositenode.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/grid/multidomaingrid.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //=======================================
    // local function space base: metaprograms
    //=======================================


template<typename T, bool isLeaf, typename E, typename It, typename Int, Dune::mdgrid::MultiDomainGridType t>
struct GuardedVisit
{
};

template<typename T, bool isLeaf, typename E, typename It, typename Int>
struct GuardedVisit<T,isLeaf,E,It,Int,Dune::mdgrid::multiDomainGrid>
{

  static void fill_indices(T& t, const E& e, It begin, Int& offset)
  {
    if (t.gridview().indexSet().contains(e))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,E,It,Int>::fill_indices(t,e,begin,offset);
  }

  static void reserve(T& t, const E& e, Int& offset)
  {
    if (t.gridview().indexSet().contains(e))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,E,It,Int>::reserve(t,e,offset);
  }

};

template<typename T, bool isLeaf, typename E, typename It, typename Int>
struct GuardedVisit<T,isLeaf,E,It,Int,Dune::mdgrid::subDomainGrid>
{

  static void fill_indices(T& t, const E& e, It begin, Int& offset)
  {
    typedef typename T::Traits::GridView::template Codim<0>::EntityPointer SDEP;
    typedef typename SDEP::Entity SDE;
    const SDEP ep = t.gridview().grid().subDomainEntityPointer(e);
    if (t.gridview().indexSet().contains(*ep))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,SDE,It,Int>::fill_indices(t,*ep,begin,offset);
  }

  static void reserve(T& t, const E& e, Int& offset)
  {
    typedef typename T::Traits::GridView::template Codim<0>::EntityPointer SDEP;
    typedef typename SDEP::Entity SDE;
    const SDEP ep = t.gridview().grid().subDomainEntityPointer(e);
    if (t.gridview().indexSet().contains(*ep))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,SDE,It,Int>::reserve(t,*ep,offset);
  }

};


template<typename T, typename E, typename It, typename Int, int n, int i>
struct MultiDomainLocalFunctionSpaceVisitChildMetaProgram // visit child of inner node
{

  typedef MultiDomainLocalFunctionSpaceVisitChildMetaProgram<T,E,It,Int,n,i+1> NextChild;

  template<typename GFS>
  static void setup (T& t, const GFS& gfs)
  {
    //        std::cout << "setting up child " << i << " of " << n << std::endl;
    t.template getChild<i>().setup(gfs.template getChild<i>());
    NextChild::setup(t,gfs);
  }

  static void fill_indices (T& t, const E& e, It begin, Int& offset)
  {
    // vist children of node t in order
    typedef typename T::template Child<i>::Type C;
    t.offset = offset;
    Int initial_offset = offset; // remember initial offset to compute size later
    GuardedVisit<C,C::isLeaf,E,It,Int,Dune::mdgrid::GridType<typename C::Traits::GridViewType::Grid>::v >::
      fill_indices(t.template getChild<i>(),e,begin,offset);
    for (Int j=initial_offset; j<offset; j++)
      begin[j] = t.pgfs->template subMap<i>(begin[j]);
    NextChild::fill_indices(t,e,begin,offset);
  }

  static void reserve (T& t, const E& e, Int& offset)
  {
    // vist children of node t in order
    typedef typename T::template Child<i>::Type C;
    GuardedVisit<C,C::isLeaf,E,It,Int,Dune::mdgrid::GridType<typename C::Traits::GridViewType::Grid>::v >::
      reserve(t.template getChild<i>(),e,offset);
    NextChild::reserve(t,e,offset);
  }
};


template<typename T, typename E, typename It, typename Int, int n>
struct MultiDomainLocalFunctionSpaceVisitChildMetaProgram<T,E,It,Int,n,n> // end of child recursion
{

  template<typename GFS>
  static void setup (T& t, const GFS& gfs)
  {
  }

  static void fill_indices (T& t, const E& e, It begin, Int& offset)
  {
    return;
  }
  static void reserve (T& t, const E& e, Int& offset)
  {
    return;
  }

};


template<typename GFS, typename N>
struct MultiDomainLocalFunctionSpaceTraits
{
  //! \brief the grid view where grid function is defined upon
  typedef GFS GridFunctionSpaceType;

  //! type of local function space node
  typedef N NodeType;

  //! \brief Type to store indices from Backend
  typedef typename GFS::Traits::GridType GridType;

  //! \brief Type of codim 0 entity in the grid
  typedef typename GridType::Traits::template Codim<0>::Entity Element;

  //! \brief Type to store indices from Backend
  typedef typename GFS::Traits::SizeType SizeType;

  //! \brief Type of container to store indices
  typedef typename std::vector<SizeType> IndexContainer;
};

template<typename... Children>
struct BuildMultiDomainLocalFunctionSpaceNodeBase
{

  struct ExtractLocalFunctionSpaceNode
  {
    template<typename GFS>
    struct transform {
      typedef typename GFS::LocalFunctionSpace::Traits::NodeType type;
    };

    template<typename... Args>
    struct container
    {
      typedef VariadicCompositeNode<CopyStoragePolicy,Args...> type;
    };

  };

  typedef typename transform<ExtractLocalFunctionSpaceNode,Children...>::type type;

};

// local function space for a power grid function space
template<typename GFS, typename... Children>
class MultiDomainLocalFunctionSpaceNode
  : public BuildMultiDomainLocalFunctionSpaceNodeBase<Children...>::type
{
  template<typename T, bool b, typename E, typename It, typename Int>
  friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
  template<typename T, typename E, typename It, typename Int, int n, int i>
  friend struct MultiDomainFunctionSpaceVisitChildMetaProgram;

  typedef typename GFS::Traits::BackendType B;
  typedef typename GFS::Traits::GridViewType::Traits::template Codim<0>::Entity Element;

  typedef typename BuildMultiDomainLocalFunctionSpaceNodeBase<Children...>::type BaseT;

public:
  typedef MultiDomainLocalFunctionSpaceTraits<GFS,MultiDomainLocalFunctionSpaceNode> Traits;

protected:
  typedef MultiDomainLocalFunctionSpaceVisitChildMetaProgram<MultiDomainLocalFunctionSpaceNode,
                                                             typename Traits::Element,
                                                             typename Traits::IndexContainer::iterator,
                                                             typename Traits::IndexContainer::size_type,
                                                             BaseT::CHILDREN,
                                                             0> VisitChildTMP;

public:

  //! \brief empty constructor
  MultiDomainLocalFunctionSpaceNode ()
  {
  }

  //! \brief initialize with grid function space
  MultiDomainLocalFunctionSpaceNode (const GFS& gfs)  : pgfs(&gfs)
  {
    setup(gfs);
  }

  //! \brief initialize with grid function space
  void setup (const GFS& gfs)
  {
    pgfs = &gfs;
    VisitChildTMP::setup(*this,gfs);
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

  // map index in this local function space to root local function space
  typename Traits::IndexContainer::size_type localIndex (typename Traits::IndexContainer::size_type index) const
  {
    return offset+index;
  }

  // map index in this local function space to global index space
  typename Traits::SizeType globalIndex (typename Traits::IndexContainer::size_type index) const
  {
    return i[index];
  }

  /** \brief extract coefficients for one element from container */
  template<typename GC, typename LC>
  void vread (const GC& globalcontainer, LC& localcontainer) const
  {
    localcontainer.resize(n);
    for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
      localcontainer[k] = B::access(globalcontainer,i[k]);
  }

  /** \brief write back coefficients for one element to container */
  template<typename GC, typename LC>
  void vwrite (const LC& localcontainer, GC& globalcontainer) const
  {
    for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
      B::access(globalcontainer,i[k]) = localcontainer[k];
  }

  /** \brief add coefficients for one element to container */
  template<typename GC, typename LC>
  void vadd (const LC& localcontainer, GC& globalcontainer) const
  {
    for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
      B::access(globalcontainer,i[k]) += localcontainer[k];
  }

  void debug () const
  {
    std::cout << n << " indices = (";
    for (typename Traits::IndexContainer::size_type k=0; k<n; k++)
      std::cout << i[k] << " ";
    std::cout << ")" << std::endl;
  }

protected:

  const GFS& gfs() const {
    return *pgfs;
  }

private:
  const GFS* pgfs;
  typename Traits::IndexContainer::iterator i;
  typename Traits::IndexContainer::size_type n;
  typename Traits::IndexContainer::size_type offset;

};



// local function space description that can be bound to an element
// depends on a grid function space
template<typename GFS,typename... Children>
class MultiDomainLocalFunctionSpace : public MultiDomainLocalFunctionSpaceNode<GFS,Children...>
{
  typedef MultiDomainLocalFunctionSpaceNode<GFS,Children...> BaseT;

  typedef typename BaseT::VisitChildTMP VisitChildTMP;

public:
  typedef typename BaseT::Traits Traits;

  MultiDomainLocalFunctionSpace (const GFS& gfs)
    : BaseT(gfs), global(gfs.maxLocalSize())
  {}

  //! \brief bind local function space to entity
  void bind (const typename Traits::Element& e)
  {
    // make offset
    typename Traits::IndexContainer::size_type offset=0;

    // compute sizes
    VisitChildTMP::reserve(*this,e,offset);

    this->n = offset;

    // now reserve space in vector
    global.resize(offset);

    // initialize iterators and fill indices
    offset = 0;
    this->offset = 0;
    this->i = global.begin();
    VisitChildTMP::fill_indices(*this,e,global.begin(),offset);

    // apply upMap
    for (typename BaseT::Traits::IndexContainer::size_type i=0; i<offset; i++)
      global[i] = this->gfs().upMap(global[i]);
  }

private:
  typename BaseT::Traits::IndexContainer global;
};




    //! \} group GridFunctionSpace

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_MULTIDOMAINLOCALFUNCTIONSPACE_HH
