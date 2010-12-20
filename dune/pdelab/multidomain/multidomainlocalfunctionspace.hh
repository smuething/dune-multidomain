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


struct MultiDomainTag {};
struct SubDomainTag {};
struct CouplingTag {};

template<typename T, bool isLeaf, typename GV, typename GC, typename Int, typename Tag>
struct StandardGFSVisitor
{

  template<typename GFS>
  static void setup(T& t, const GFS& gfs)
  {}

  template<typename E>
  static void fill_indices(T& t, GV gv, const E& e, Int& offset, GC * const global)
  {}

  template<typename E>
  static void compute_size(T& t, GV gv, const E& e, Int& offset)
  {}

};

template<typename T, bool isLeaf, typename GV, typename GC, typename Int>
struct StandardGFSVisitor<T,isLeaf,GV,GC,Int,MultiDomainTag>
{

  template<typename GFS>
  static void setup(T& t, const GFS& gfs)
  {
    t.setup(gfs);
  }

  template<typename E>
  static void fill_indices(T& t, GV gv, const E& e, Int& offset, GC * const global)
  {
    if (gv.indexSet().contains(e))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,E,GC,Int>().fill_indices(t,e,offset,global);
  }

  template<typename E>
  static void compute_size(T& t, GV gv, const E& e, Int& offset)
  {
    if (gv.indexSet().contains(e))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,E,GC,Int>().compute_size(t,e,offset);
  }

};

template<typename T, bool isLeaf, typename GV, typename GC, typename Int>
struct StandardGFSVisitor<T,isLeaf,GV,GC,Int,SubDomainTag>
{

  template<typename GFS>
  static void setup(T& t, const GFS& gfs)
  {
    t.setup(gfs);
  }

  template<typename E>
  static void fill_indices(T& t, GV gv, const E& e, Int& offset, GC * const global)
  {
    typedef typename T::Traits::GridViewType::template Codim<0>::EntityPointer SDEP;
    typedef typename SDEP::Entity SDE;
    const SDEP ep = gv.grid().subDomainEntityPointer(e);
    if (gv.indexSet().contains(*ep))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,SDE,GC,Int>().fill_indices(t,*ep,offset,global);
  }

  template<typename E>
  static void compute_size(T& t, GV gv, const E& e, Int& offset)
  {
    typedef typename T::Traits::GridViewType::template Codim<0>::EntityPointer SDEP;
    typedef typename SDEP::Entity SDE;
    const SDEP ep = gv.grid().subDomainEntityPointer(e);
    if (gv.indexSet().contains(*ep))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,SDE,GC,Int>().compute_size(t,*ep,offset);
  }

};


template<typename T, bool isLeaf, typename GV, typename GC, typename Int, typename Tag>
struct CouplingGFSVisitor
{

  template<typename GFS>
  static void setup(T& t, const GFS& gfs)
  {}

  template<typename Intersection>
  static void fill_indices(T& t, GV gv, const Intersection& is, Int& offset, GC * const global)
  {}

  template<typename Intersection>
  static void compute_size(T& t, GV gv, const Intersection& is, Int& offset)
  {}

};


template<typename T, bool isLeaf, typename GV, typename GC, typename Int>
struct CouplingGFSVisitor<T,isLeaf,GV,GC,Int,CouplingTag>
{

  template<typename GFS>
  static void setup(T& t, const GFS& gfs)
  {
    t.setup(gfs);
  }

  template<typename Intersection>
  static void fill_indices(T& t, GV gv, const Intersection& is, Int& offset, GC * const global)
  {
    if (t.gridFunctionSpace().contains(is))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,Intersection,GC,Int>().fill_indices(t,is,offset,global);
  }

  template<typename Intersection>
  static void compute_size(T& t, GV gv, const Intersection& is, Int& offset)
  {
    if (t.gridFunctionSpace().contains(is))
      LocalFunctionSpaceBaseVisitNodeMetaProgram<T,isLeaf,Intersection,GC,Int>().reserve(t,is,offset);
  }

};


template<typename T,
         typename Container,
         typename Int,
         template<typename T1, bool isLeaf, typename GV, typename C1, typename Int1, typename Tag> class Visitor,
         int n,
         int i
         >
struct MultiDomainLocalFunctionSpaceVisitChildMetaProgram // visit child of inner node
{

  typedef MultiDomainLocalFunctionSpaceVisitChildMetaProgram<T,Container,Int,Visitor,n,i+1> NextChild;
  typedef typename T::Traits::GridFunctionSpaceType GFS;

  template<typename GFS>
  static void setup (T& t, const GFS& gfs)
  {
    //        std::cout << "setting up child " << i << " of " << n << std::endl;
    typedef typename T::template Child<i>::Type C;
    Visitor<C,C::isLeaf,typename C::Traits::GridViewType,Container,Int,typename GFS::template ChildInfo<i>::Type::Tag >::
      setup(t.template getChild<i>(),gfs.template getChild<i>());
    NextChild::setup(t,gfs);
  }

  template<typename E>
  static void fill_indices (T& t, const E& e, Int& offset, Container * const global)
  {
    // vist children of node t in order
    typedef typename T::template Child<i>::Type C;
    Int initial_offset = offset; // remember initial offset to compute size later
    Visitor<C,C::isLeaf,typename C::Traits::GridViewType,Container,Int,typename GFS::template ChildInfo<i>::Type::Tag >::
      fill_indices(t.template getChild<i>(),t.gfs().template getChild<i>().gridview(),e,offset,global);
    for (Int j=initial_offset; j<offset; j++)
      (*global)[initial_offset+j] = t.pgfs->template subMap<i>((*global)[initial_offset+j]);
    NextChild::fill_indices(t,e,offset,global);
  }

  template<typename E>
  static void compute_size (T& t, const E& e, Int& offset)
  {
    // vist children of node t in order
    typedef typename T::template Child<i>::Type C;
    Visitor<C,C::isLeaf,typename C::Traits::GridViewType,Container,Int,typename GFS::template ChildInfo<i>::Type::Tag >::
      compute_size(t.template getChild<i>(),t.gfs().template getChild<i>().gridview(),e,offset);
    NextChild::compute_size(t,e,offset);
  }
};


template<typename T,
         typename Container,
         typename Int,
         template<typename T1, bool isLeaf, typename GV, typename C1, typename Int1, typename Tag> class Visitor,
         int n
         >
struct MultiDomainLocalFunctionSpaceVisitChildMetaProgram<T,Container,Int,Visitor,n,n> // end of child recursion
{

  template<typename GFS>
  static void setup (T& t, const GFS& gfs)
  {
  }

  template<typename E>
  static void fill_indices (T& t, const E& e, Int& offset, Container * const global)
  {
    return;
  }

  template<typename E>
  static void compute_size (T& t, const E& e, Int& offset)
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

  typedef typename GFS::Traits::GridViewType GridViewType;

  //! \brief Type of codim 0 entity in the grid
  typedef typename GridType::Traits::template Codim<0>::Entity Element;

  //! \brief Type of intersection in the grid
  typedef typename GridViewType::Intersection Intersection;

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
template<typename GFS,
         template<typename, bool, typename, typename, typename, typename> class Visitor,
         typename... Children>
class MultiDomainLocalFunctionSpaceNode
  : public BuildMultiDomainLocalFunctionSpaceNodeBase<Children...>::type
{
  template<typename T, bool b, typename E, typename It, typename Int>
  friend struct LocalFunctionSpaceBaseVisitNodeMetaProgram;
  template<typename T,
           typename It,
           typename Int,
           template<typename T1, bool isLeaf, typename GV, typename It1, typename Int1, typename Tag> class Visitor1,
           int n,
           int i
           >
  friend struct MultiDomainLocalFunctionSpaceVisitChildMetaProgram;

  typedef typename GFS::Traits::BackendType B;
  typedef typename GFS::Traits::GridType::Traits::template Codim<0>::Entity Element;

  typedef typename BuildMultiDomainLocalFunctionSpaceNodeBase<Children...>::type BaseT;

public:
  typedef MultiDomainLocalFunctionSpaceTraits<GFS,MultiDomainLocalFunctionSpaceNode> Traits;

protected:
  typedef MultiDomainLocalFunctionSpaceVisitChildMetaProgram<MultiDomainLocalFunctionSpaceNode,
                                                             typename Traits::IndexContainer,
                                                             typename Traits::IndexContainer::size_type,
                                                             Visitor,
                                                             BaseT::CHILDREN,
                                                             0> VisitChildTMP;

public:

  //! \brief empty constructor
  MultiDomainLocalFunctionSpaceNode ()
  {
  }

  //! \brief initialize with grid function space
  MultiDomainLocalFunctionSpaceNode (const GFS& gfs) : pgfs(&gfs)
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

  typename Traits::IndexContainer::size_type localVectorSize() const
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
    return global()[offset + index];
  }

  /** \brief extract coefficients for one element from container */
  template<typename GC, typename LC>
  void vread (const GC& globalcontainer, LC& localcontainer) const
  {
    localcontainer.resize(n);
    for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
      localcontainer[k] = B::access(globalcontainer,global()[offset + k]);
  }

  /** \brief write back coefficients for one element to container */
  template<typename GC, typename LC>
  void vwrite (const LC& localcontainer, GC& globalcontainer) const
  {
    for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
      B::access(globalcontainer,global()[offset+k]) = localcontainer[k];
  }

  /** \brief add coefficients for one element to container */
  template<typename GC, typename LC>
  void vadd (const LC& localcontainer, GC& globalcontainer) const
  {
    for (typename Traits::IndexContainer::size_type k=0; k<n; ++k)
      B::access(globalcontainer,global()[offset+k]) += localcontainer[k];
  }

  void debug () const
  {
    std::cout << n << " indices = (";
    for (typename Traits::IndexContainer::size_type k=0; k<n; k++)
      std::cout << global()[offset+k] << " ";
    std::cout << ")" << std::endl;
  }

public:

  const GFS& gfs() const {
    return *pgfs;
  }

protected:
  const GFS* pgfs;
  typename Traits::IndexContainer::size_type n;
  typename Traits::IndexContainer::size_type offset;
  typename Traits::IndexContainer* _global;

  typename Traits::IndexContainer& global()
  {
    return *_global;
  }

  const typename Traits::IndexContainer& global() const
  {
    return *_global;
  }

};



// local function space description that can be bound to an element
// depends on a grid function space
template<typename GFS,typename... Children>
class MultiDomainLocalFunctionSpace : public MultiDomainLocalFunctionSpaceNode<GFS,StandardGFSVisitor,Children...>
{
  typedef MultiDomainLocalFunctionSpaceNode<GFS,StandardGFSVisitor,Children...> BaseT;

  typedef typename BaseT::VisitChildTMP VisitChildTMP;

public:
  typedef typename BaseT::Traits Traits;

  explicit MultiDomainLocalFunctionSpace (const GFS& gfs)
    : BaseT(gfs), global_container(gfs.maxLocalSize())
  {}

  //! \brief bind local function space to entity
  void bind (const typename Traits::Element& e)
  {
    // make offset
    typename Traits::IndexContainer::size_type offset=0;
    this->_global = &global_container;

    // compute sizes
    VisitChildTMP::compute_size(*this,e,offset);

    this->n = offset;

    // now reserve space in vector
    global_container.resize(offset);

    // initialize iterators and fill indices
    offset = 0;
    this->offset = 0;
    VisitChildTMP::fill_indices(*this,e,offset,this->_global);

    // apply upMap
    for (typename Traits::IndexContainer::size_type i=0; i<offset; i++)
      global_container[i] = this->gfs().upMap(global_container[i]);
  }

private:
  typename Traits::IndexContainer global_container;
};


// local function space description that can be bound to an intersection
// depends on a grid function space
template<typename GFS,typename... Children>
class MultiDomainCouplingLocalFunctionSpace : public MultiDomainLocalFunctionSpaceNode<GFS,CouplingGFSVisitor,Children...>
{
  typedef MultiDomainLocalFunctionSpaceNode<GFS,CouplingGFSVisitor,Children...> BaseT;

  typedef typename BaseT::VisitChildTMP VisitChildTMP;

public:
  typedef typename BaseT::Traits Traits;

  explicit MultiDomainCouplingLocalFunctionSpace (const GFS& gfs)
    : BaseT(gfs), global_container(gfs.maxLocalSize())
  {}

  //! \brief bind local function space to entity
  void bind (const typename Traits::Intersection& is)
  {
    // make offset
    typename Traits::IndexContainer::size_type offset=0;
    this->_global = &global_container;

    // compute sizes
    VisitChildTMP::compute_size(*this,is,offset);

    this->n = offset;

    // now reserve space in vector
    global_container.resize(offset);

    // initialize iterators and fill indices
    offset = 0;
    this->offset = 0;
    VisitChildTMP::fill_indices(*this,is,offset,this->_global);

    // apply upMap
    for (typename Traits::IndexContainer::size_type i=0; i<offset; i++)
      global_container[i] = this->gfs().upMap(global_container[i]);
  }

private:
  typename Traits::IndexContainer global_container;
};



    //! \} group GridFunctionSpace

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_MULTIDOMAINLOCALFUNCTIONSPACE_HH
