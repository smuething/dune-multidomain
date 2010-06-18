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
      fill_indices(t.template getChild<i>(),t.pgfs->template getChild<i>().gridview(),e,begin,offset);
    for (Int j=initial_offset; j<offset; j++)
      begin[j] = t.pgfs->template subMap<i>(begin[j]);
    NextChild::fill_indices(t,e,begin,offset);
  }

  static void reserve (T& t, const E& e, Int& offset)
  {
    // vist children of node t in order
    typedef typename T::template Child<i>::Type C;
    GuardedVisit<C,C::isLeaf,typename C::Traits::GridViewType,E,It,Int,Dune::mdgrid::GridType<typename C::Traits::GridViewType::Grid>::v >::
      reserve(t.template getChild<i>(),t.pgfs->template getChild<i>().gridview(),e,offset);
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


template<typename GFS, typename N, typename _SubProblem, typename _Constraints>
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

  typedef _SubProblem SubProblem;

  typedef typename SubProblem::Traits::Condition Condition;

  typedef _Constraints Constraints;

  typedef Constraints ConstraintsType;

};


template<typename GFS, typename N, typename BaseLFS, typename _SubProblem, typename _Constraints>
struct SubProblemLeafLocalFunctionSpaceTraits
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

  //! \brief local finite element
  typedef typename BaseLFS::Traits::LocalFiniteElementType LocalFiniteElementType;

  typedef _SubProblem SubProblem;

  typedef typename SubProblem::Traits::Condition Condition;

  typedef _Constraints Constraints;

  typedef Constraints ConstraintsType;

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
    typedef typename build_splfs_node<MDLFS,ChildIndices...>::template result<Types..., const typename MDLFS::template Child<FirstIndex>::Type>::type type;
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

template<typename MDLFS, typename VariadicNode, int... ChildIndices>
struct SubProblemLocalFunctionSpaceBase;


template<typename MDLFS, typename VariadicNode, int i, int... ChildIndices>
struct SubProblemLocalFunctionSpaceBase<MDLFS,VariadicNode,i,ChildIndices...> :
  public SubProblemLocalFunctionSpaceBase<MDLFS,VariadicNode,ChildIndices...>
{

  template<typename... Children>
  SubProblemLocalFunctionSpaceBase(const MDLFS& mdlfs, Children&... children) :
    SubProblemLocalFunctionSpaceBase<MDLFS,VariadicNode,ChildIndices...>(mdlfs, children..., mdlfs.template getChild<i>())
  {}

};

template<typename MDLFS, typename VariadicNode>
struct SubProblemLocalFunctionSpaceBase<MDLFS,VariadicNode> :
  public VariadicNode
{

  template<typename... Children>
  SubProblemLocalFunctionSpaceBase(const MDLFS& mdlfs, Children&... children) :
    VariadicNode(children...)
  {}

};

} // anonymous namespace


// ********************************************************************************
// LocalFunctionSpace for subproblems

template<typename MDLFS, typename SubProblem, typename Constraints, int... ChildIndices>
class SubProblemLocalFunctionSpace
  : public SubProblemLocalFunctionSpaceBase<MDLFS,
                                            typename build_splfs_node<MDLFS,ChildIndices...>::template result<>::type,
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
                                           ChildIndices...> BaseT;

public:
  typedef typename Dune::SelectType<sizeof...(ChildIndices) == 1,
                                    SubProblemLeafLocalFunctionSpaceTraits<GFS,SubProblemLocalFunctionSpace,typename BaseT::template Child<0>::Type,SubProblem,Constraints>,
                                    SubProblemLocalFunctionSpaceTraits<GFS,SubProblemLocalFunctionSpace,SubProblem,Constraints>
                                    >::Type Traits;

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
  SubProblemLocalFunctionSpace (const MDLFS& mdlfs, const SubProblem& subProblem, const Constraints& constraints) :
    BaseT(mdlfs),
    plfs(&mdlfs),
    pgfs(&(mdlfs.gfs())),
    _subProblem(subProblem),
    _constraints(constraints),
    n(VisitChildTMP::size(*this))
  {
    //setup(mdlfs);
  }

  //! \brief initialize with grid function space
  void setup (const MDLFS& lfs)
  {
    assert(false);
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
    typename Traits::IndexContainer::size_type offset2=0;

    // compute sizes
    VisitChildTMP::reserve(*this,e,offset2);

    this->n = offset2;

    // now reserve space in vector
    global.resize(offset2);

    // initialize iterators and fill indices
    offset2 = 0;
    this->offset = 0;
    this->i = global.begin();
    VisitChildTMP::fill_indices(*this,e,global.begin(),offset2);

    // apply upMap
    for (typename Traits::IndexContainer::size_type i=0; i<offset2; i++)
      global[i] = pgfs->upMap(global[i]);
  }

  const SubProblem& subProblem() const {
    return _subProblem;
  }

  const Constraints& constraints() const {
    return _constraints();
  }

  template<typename SubDomainSet>
  bool appliesTo(const SubDomainSet& sds) const {
    return _subProblem.appliesTo(sds);
  }

  const typename Traits::LocalFiniteElementType& localFiniteElement() const {
    return this->template getChild<0>().localFiniteElement();
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
