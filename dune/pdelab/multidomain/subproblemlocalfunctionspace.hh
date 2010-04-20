// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_MULTIDOMAIN_SUBPROBLEMLOCALFUNCTIONSPACE_HH
#define DUNE_MULTIDOMAIN_SUBPROBLEMLOCALFUNCTIONSPACE_HH

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

/*
template<typename T, bool isLeaf, typename GV, typename E, typename It, typename Int, Dune::mdgrid::MultiDomainGridType t>
struct GuardedVisit
{
};

template<typename T, bool isLeaf, typename GV, typename E, typename It, typename Int>
struct GuardedVisit<T,isLeaf,GV,E,It,Int,Dune::mdgrid::multiDomainGrid>
{

  static void fill_indices(T& t, GV gv, const E& e, It begin, Int& offset)
  {
    if (gv.indexSet().contains(e))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,E,It,Int>::fill_indices(t,e,begin,offset);
  }

  static void reserve(T& t, GV gv, const E& e, Int& offset)
  {
    if (gv.indexSet().contains(e))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,E,It,Int>::reserve(t,e,offset);
  }

};

template<typename T, bool isLeaf, typename GV, typename E, typename It, typename Int>
struct GuardedVisit<T,isLeaf,GV,E,It,Int,Dune::mdgrid::subDomainGrid>
{

  static void fill_indices(T& t, GV gv, const E& e, It begin, Int& offset)
  {
    typedef typename T::Traits::GridViewType::template Codim<0>::EntityPointer SDEP;
    typedef typename SDEP::Entity SDE;
    const SDEP ep = gv.grid().subDomainEntityPointer(e);
    if (gv.indexSet().contains(*ep))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,SDE,It,Int>::fill_indices(t,*ep,begin,offset);
  }

  static void reserve(T& t, GV gv, const E& e, Int& offset)
  {
    typedef typename T::Traits::GridViewType::template Codim<0>::EntityPointer SDEP;
    typedef typename SDEP::Entity SDE;
    const SDEP ep = gv.grid().subDomainEntityPointer(e);
    if (gv.indexSet().contains(*ep))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,SDE,It,Int>::reserve(t,*ep,offset);
  }

};

*/

template<typename T, typename E, typename It, typename Int, int n, int i>
struct SubProblemLocalFunctionSpaceVisitChildMetaProgram // visit child of inner node
{

  typedef SubProblemLocalFunctionSpaceVisitChildMetaProgram<T,E,It,Int,n,i+1> NextChild;

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
    GuardedVisit<C,C::isLeaf,typename C::Traits::GridViewType,E,It,Int,Dune::mdgrid::GridType<typename C::Traits::GridViewType::Grid>::v >::
      fill_indices(t.template getChild<i>(),t.gfs().template getChild<i>().gridview(),e,begin,offset);
    for (Int j=initial_offset; j<offset; j++)
      begin[j] = t.pgfs->template subMap<i>(begin[j]);
    NextChild::fill_indices(t,e,begin,offset);
  }

  static void reserve (T& t, const E& e, Int& offset)
  {
    // vist children of node t in order
    typedef typename T::template Child<i>::Type C;
    GuardedVisit<C,C::isLeaf,typename C::Traits::GridViewType,E,It,Int,Dune::mdgrid::GridType<typename C::Traits::GridViewType::Grid>::v >::
      reserve(t.template getChild<i>(),t.gfs().template getChild<i>().gridview(),e,offset);
    NextChild::reserve(t,e,offset);
  }
};


template<typename T, typename E, typename It, typename Int, int n>
struct SubProblemLocalFunctionSpaceVisitChildMetaProgram<T,E,It,Int,n,n> // end of child recursion
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
struct SubProblemLocalFunctionSpaceTraits
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

namespace {

// *************************************************************************************
// TMP for constructing the type of the variadic node the SubProblemLocalFunctionSpace inherits from

template<typename MDLFS, int... ChildIndices>
struct build_splfs_node;

template<typename MDLFS, int FirstIndex, int... ChildIndices>
struct build_splfs_node<MDLFS,FirstIndex,ChildIndices...>
{
  template<typename... Types>
  struct result
  {
    typedef typename build_splfs_node<MDLFS,ChildIndices...>::template result<Types..., typename MDLFS::template Child<FirstIndex>::Type>::type type;
  };
};

template<typename MDLFS>
struct build_splfs_node<MDLFS>
{
  template<typename... Types>
  struct result
  {
    typedef VariadicCompositeNode<CopyStoragePolicy,Types...> type;
  };
};


// *************************************************************************************
// System of intermediate base classes for SubProblemLocalFunctionSpace that extracts the
// selected children from the MultiDomainGridFunctionSpace so they can be passed on to the
// constructor of the VariadicNode

template<typename MDLFS, typename VariadicNode, int i, int... ChildIndices>
struct SubProblemLocalFunctionSpaceBase;

template<typename MDLFS, typename VariadicNode, int i, int FirstIndex, int... ChildIndices>
struct SubProblemLocalFunctionSpaceBase<MDLFS,VariadicNode,i,FirstIndex, ChildIndices...> :
  public SubProblemLocalFunctionSpaceBase<MDLFS,VariadicNode,i+1,FirstIndex,ChildIndices...>
{

  template<typename... Children>
  SubProblemLocalFunctionSpaceBase(MDLFS& mdlfs, Children&&... children) :
    SubProblemLocalFunctionSpaceBase<MDLFS,VariadicNode,i+1,FirstIndex,ChildIndices...>(mdlfs, children...)
  {}

};

template<typename MDLFS, typename VariadicNode, int i, int... ChildIndices>
struct SubProblemLocalFunctionSpaceBase<MDLFS,VariadicNode,i,i,ChildIndices...> :
  public SubProblemLocalFunctionSpaceBase<MDLFS,VariadicNode,i+1,ChildIndices...>
{

  template<typename... Children>
  SubProblemLocalFunctionSpaceBase(MDLFS& mdlfs, Children&&... children) :
    SubProblemLocalFunctionSpaceBase<MDLFS,VariadicNode,i+1,ChildIndices...>(mdlfs, children..., mdlfs.template getChild<i>())
  {}

};

template<typename MDLFS, typename VariadicNode, int i>
struct SubProblemLocalFunctionSpaceBase<MDLFS,VariadicNode,i> :
  public VariadicNode
{

  template<typename... Children>
  SubProblemLocalFunctionSpaceBase(MDLFS& mdlfs, Children&&... children) :
    VariadicNode(children...)
  {}

};

} // anonymous namespace


// ********************************************************************************
// LocalFunctionSpace for subproblems

template<typename MDLFS, typename Condition, int... ChildIndices>
class SubProblemLocalFunctionSpace
  : public SubProblemLocalFunctionSpaceBase<MDLFS,
                                            typename build_splfs_node<MDLFS,ChildIndices...>::template result<>::type,
                                            0,
                                            ChildIndices...>
{

  dune_static_assert((sizeof...(ChildIndices) <= MDLFS::CHILDREN),"SubProblemLocalFunctionSpace cannot have more components than the MultiDomainGridFunctionSpace");

  dune_static_assert((check_indices<MDLFS::CHILDREN, ChildIndices...>::value),"Invalid set of child indices (index out of range or duplicate indices)");

  typedef typename MDLFS::Traits::GridFunctionSpaceType GFS;

  template<typename T, bool b, typename E, typename It, typename Int>
  friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
  template<typename T, typename E, typename It, typename Int, int n, int i>
  friend struct SubProblemLocalFunctionSpaceVisitChildMetaProgram;

  typedef typename GFS::Traits::BackendType B;
  typedef typename GFS::Traits::GridType::Traits::template Codim<0>::Entity Element;

  typedef SubProblemLocalFunctionSpaceBase<MDLFS,
                                           typename build_splfs_node<MDLFS,ChildIndices...>::template result<>::type,
                                           0,
                                           ChildIndices...> BaseT;

public:
  typedef SubProblemLocalFunctionSpaceTraits<GFS,SubProblemLocalFunctionSpace> Traits;

protected:
  typedef SubProblemLocalFunctionSpaceVisitChildMetaProgram<SubProblemLocalFunctionSpace,
                                                            typename Traits::Element,
                                                            typename Traits::IndexContainer::iterator,
                                                            typename Traits::IndexContainer::size_type,
                                                            BaseT::CHILDREN,
                                                            0> VisitChildTMP;

public:

  //! \brief empty constructor
  SubProblemLocalFunctionSpace ()
  {
  }

  //! \brief initialize with grid function space
  SubProblemLocalFunctionSpace (MDLFS& mdlfs, const Condition& condition_) :
    BaseT(mdlfs),
    plfs(&mdlfs),
    pgfs(&(mdlfs.gfs())),
    condition(condition_)
  {
    setup(mdlfs);
  }

  //! \brief initialize with grid function space
  void setup (const MDLFS& lfs)
  {
    plfs = &lfs;
    pgfs = &(lfs.gfs());
    VisitChildTMP::setup(*this,*pgfs);
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
      global[i] = pgfs->upMap(global[i]);
  }

private:
  const MDLFS* plfs;
  const GFS* pgfs;
  const Condition condition;
  typename Traits::IndexContainer::iterator i;
  typename Traits::IndexContainer::size_type n;
  typename Traits::IndexContainer::size_type offset;
  typename Traits::IndexContainer global;

};


    //! \} group GridFunctionSpace

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_SUBPROBLEMLOCALFUNCTIONSPACE_HH
