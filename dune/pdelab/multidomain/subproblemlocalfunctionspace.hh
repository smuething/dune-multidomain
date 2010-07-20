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

  static std::size_t size(const T& t)
  {
    return t.template getChild<i>().size() + NextChild::size(t);
  }

  static void fill_indices (T& t, It it)
  {
    typedef typename T::template Child<i>::Type C;
    const C& child = t.template getChild<i>();
    for(unsigned int j = 0; j < child.size(); ++j, ++it)
        *it = child.globalIndex(j);
    NextChild::fill_indices(t,it);
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

  static std::size_t size(const T& t)
  {
    return 0;
  }

  static void fill_indices (T& t, It it)
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
// compile-time mapper from subproblem child indices to global child indices

template<bool b, int v, typename next, int i>
struct pick_value
{
  static const int value = v;
};

template<int v, typename next, int i>
struct pick_value<false,v,next,i>
{
  static const int value = next::template access<i>::value;
};

template<int i, int... values>
struct value_mapper_impl;

template<int i, int v, int... values>
struct value_mapper_impl<i,v,values...>
{
  static const int _key = i;
  static const int _value = v;

  template<int j>
  struct access
  {
    static const int value = pick_value<_key == j,_value,value_mapper_impl<_key+1,values...>,j>::value;
  };
};

template<int i>
struct value_mapper_impl<i>
{
  // end of recursion
};

template<int... values>
struct value_mapper
{
  template<int i>
  struct map {
    dune_static_assert(i < sizeof...(values),"index out of range");
    static const int value = value_mapper_impl<0,values...>::template access<i>::value;
  };
};

} // anonymous namespace


// ********************************************************************************
// LocalFunctionSpace for subproblems - multi-component version

template<typename MDLFS, typename SubProblem, typename Constraints, int... ChildIndices>
class SubProblemLocalFunctionSpace
{

  dune_static_assert((sizeof...(ChildIndices) <= MDLFS::CHILDREN),"SubProblemLocalFunctionSpace cannot have more components than the MultiDomainGridFunctionSpace");

  dune_static_assert((check_indices<MDLFS::CHILDREN, ChildIndices...>::value),"Invalid set of child indices (index out of range or duplicate indices)");

  typedef value_mapper<ChildIndices...> index_mapper;

public:
  static const bool isLeaf = false;
  static const bool isPower = false;
  static const bool isComposite = true;
  static const unsigned int CHILDREN = sizeof...(ChildIndices);

  template<int i>
  struct Child
  {
    typedef typename MDLFS::template Child<index_mapper::template map<i>::value>::Type Type;
  };

  template<int i>
  const typename Child<i>::Type& getChild () const
  {
    return plfs->template getChild<index_mapper::template map<i>::value>();
  }

private:

  typedef typename MDLFS::Traits::GridFunctionSpaceType GFS;

  template<typename T, bool b, typename E, typename It, typename Int>
  friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
  template<typename T, typename E, typename It, typename Int, int n, int i>
  friend struct SubProblemLocalFunctionSpaceVisitChildMetaProgram;

  typedef typename GFS::Traits::BackendType B;
  typedef typename GFS::Traits::GridType::Traits::template Codim<0>::Entity Element;

public:
  typedef SubProblemLocalFunctionSpaceTraits<GFS,SubProblemLocalFunctionSpace,SubProblem,Constraints> Traits;

protected:
  typedef SubProblemLocalFunctionSpaceVisitChildMetaProgram<const SubProblemLocalFunctionSpace,
                                                            typename Traits::Element,
                                                            typename Traits::IndexContainer::iterator,
                                                            typename Traits::IndexContainer::size_type,
                                                            CHILDREN,
                                                            0> VisitChildTMP;

public:

  //! \brief empty constructor - TODO:Do we need this???
  /*SubProblemLocalFunctionSpace ()
  {
  }*/

  //! \brief initialize with grid function space
  SubProblemLocalFunctionSpace (const MDLFS& mdlfs, const SubProblem& subProblem, const Constraints& constraints) :
    plfs(&mdlfs),
    pgfs(&(mdlfs.gfs())),
    _subProblem(subProblem),
    _constraints(constraints),
    n(VisitChildTMP::size(*this))
  {
    setup(mdlfs);
  }

  //! \brief variant if no local function space is availabe at construction time
  SubProblemLocalFunctionSpace(const SubProblem& subProblem, const Constraints& constraints) :
    plfs(NULL),
    pgfs(NULL),
    _subProblem(subProblem),
    _constraints(constraints)
  {}

  //! \brief initialize with grid function space
  void setup (const MDLFS& lfs) const
  {
    plfs = &lfs;
    pgfs = &(lfs.gfs());
    /*
    VisitChildTMP::setup(*this,*pgfs);
    */

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
    return index;
  }

  // map index in this local function space to global index space
  typename Traits::SizeType globalIndex (typename Traits::IndexContainer::size_type index) const
  {
    return global[index];
  }

  /** \brief extract coefficients for one element from container */
  template<typename GC, typename LC>
  void vread (const GC& globalcontainer, LC& localcontainer) const
  {
    localcontainer.resize(n);
    for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
      localcontainer[k] = B::access(globalcontainer,global[k]);
  }

  /** \brief write back coefficients for one element to container */
  template<typename GC, typename LC>
  void vwrite (const LC& localcontainer, GC& globalcontainer) const
  {
    for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
      B::access(globalcontainer,global[k]) = localcontainer[k];
  }

  /** \brief add coefficients for one element to container */
  template<typename GC, typename LC>
  void vadd (const LC& localcontainer, GC& globalcontainer) const
  {
    for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
      B::access(globalcontainer,global[k]) += localcontainer[k];
  }

  void debug () const
  {
    std::cout << n << " indices = (";
    for (typename Traits::IndexContainer::size_type k=0; k<n; k++)
      std::cout << global[k] << " ";
    std::cout << ")" << std::endl;
  }

  //! \brief bind local function space to entity
  void bind () const
  {
    // make offset
    n = VisitChildTMP::size(*this);
    global.resize(n);
    VisitChildTMP::fill_indices(*this,global.begin());
  }

  const SubProblem& subProblem() const {
    return _subProblem;
  }

  const Constraints& constraints() const {
    return _constraints;
  }

  template<typename SubDomainSet>
  bool appliesTo(const SubDomainSet& sds) const {
    return _subProblem.appliesTo(sds);
  }

private:
  mutable const MDLFS* plfs;
  mutable const GFS* pgfs;
  const SubProblem& _subProblem;
  const Constraints& _constraints;
  mutable typename Traits::IndexContainer::size_type n;
  mutable typename Traits::IndexContainer global;

};


// single-component version


template<typename MDLFS, typename SubProblem, typename Constraints, int ChildIndex>
class SubProblemLocalFunctionSpace<MDLFS,SubProblem,Constraints,ChildIndex>
  : public LeafNode
{

  dune_static_assert((ChildIndex < MDLFS::CHILDREN),"Child index out of range");

  typedef typename MDLFS::Traits::GridFunctionSpaceType GFS;

  template<typename T, bool b, typename E, typename It, typename Int>
  friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
  template<typename T, typename E, typename It, typename Int, int n, int i>
  friend struct SubProblemLocalFunctionSpaceVisitChildMetaProgram;

  typedef typename GFS::Traits::BackendType B;
  typedef typename GFS::Traits::GridType::Traits::template Codim<0>::Entity Element;

  typedef typename MDLFS::template Child<ChildIndex>::Type BaseLFS;

  const BaseLFS& baseLFS() const {
    return plfs->template getChild<ChildIndex>();
  }

public:
  typedef SubProblemLeafLocalFunctionSpaceTraits<GFS,SubProblemLocalFunctionSpace,BaseLFS,SubProblem,Constraints> Traits;

  //! \brief empty constructor - TODO:Do we need this???
  /*SubProblemLocalFunctionSpace ()
  {
  }*/

  //! \brief initialize with grid function space
  SubProblemLocalFunctionSpace (const MDLFS& mdlfs, const SubProblem& subProblem, const Constraints& constraints) :
    plfs(&mdlfs),
    pgfs(&(mdlfs.gfs())),
    _subProblem(subProblem),
    _constraints(constraints)
  {
    setup(mdlfs);
  }

  //! \brief variant if no local function space is availabe at construction time
  SubProblemLocalFunctionSpace(const SubProblem& subProblem, const Constraints& constraints) :
    plfs(NULL),
    pgfs(NULL),
    _subProblem(subProblem),
    _constraints(constraints)
  {}

  //! \brief initialize with grid function space
  void setup (const MDLFS& lfs) const
  {
    /*
      assert(false);*/
    plfs = &lfs;
    pgfs = &(lfs.gfs());/*
    VisitChildTMP::setup(*this,*pgfs);
    */
  }

  void bind() const
  {
    // nothing to do here...
  }

  //! \brief get current size
  typename Traits::IndexContainer::size_type size () const
  {
    return baseLFS().size();
  }

  //! \brief get maximum possible size (which is maxLocalSize from grid function space)
  typename Traits::IndexContainer::size_type maxSize () const
  {
    return pgfs->maxLocalSize();
  }

  // map index in this local function space to root local function space
  typename Traits::IndexContainer::size_type localIndex (typename Traits::IndexContainer::size_type index) const
  {
    return baseLFS().localIndex(index);
  }

  // map index in this local function space to global index space
  typename Traits::SizeType globalIndex (typename Traits::IndexContainer::size_type index) const
  {
    return baseLFS().globalIndex(index);
  }

  /** \brief extract coefficients for one element from container */
  template<typename GC, typename LC>
  void vread (const GC& globalcontainer, LC& localcontainer) const
  {
    baseLFS().vread(globalcontainer,localcontainer);
  }

  /** \brief write back coefficients for one element to container */
  template<typename GC, typename LC>
  void vwrite (const LC& localcontainer, GC& globalcontainer) const
  {
    baseLFS().vwrite(localcontainer,globalcontainer);
  }

  /** \brief add coefficients for one element to container */
  template<typename GC, typename LC>
  void vadd (const LC& localcontainer, GC& globalcontainer) const
  {
    baseLFS().vadd(localcontainer,globalcontainer);
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
    assert(false);
  }

  template<typename GC, typename LC>
  void mwrite (const LC& lc, GC& gc) const
  {
    baseLFS().mwrite(lc,gc);
  }

  const typename Traits::LocalFiniteElementType& localFiniteElement() const {
    return baseLFS().localFiniteElement();
  }

  const SubProblem& subProblem() const {
    return _subProblem;
  }

  const Constraints& constraints() const {
    return _constraints;
  }

  template<typename SubDomainSet>
  bool appliesTo(const SubDomainSet& sds) const {
    return _subProblem.appliesTo(sds);
  }

private:
  mutable const MDLFS * plfs;
  mutable const GFS * pgfs;
  const SubProblem& _subProblem;
  const Constraints& _constraints;
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
