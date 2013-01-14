// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_MULTIDOMAIN_SUBPROBLEMLOCALFUNCTIONSPACE_HH
#define DUNE_MULTIDOMAIN_SUBPROBLEMLOCALFUNCTIONSPACE_HH

#include <vector>
#include <dune/pdelab/multidomain/utility.hh>
#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/common/typetree/proxynode.hh>
#include <dune/pdelab/common/typetree/filteredcompositenode.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/grid/multidomaingrid.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

namespace {

  struct AccumulateSize
    : public Dune::PDELab::TypeTree::DirectChildrenVisitor
    , public Dune::PDELab::TypeTree::DynamicTraversal
  {

    template<typename SPLFS, typename ChildLFS, typename TreePath, typename ChildIndex>
    void afterChild(SPLFS& splfs, const ChildLFS& childLFS, TreePath treePath, ChildIndex childIndex)
    {
      size += childLFS.size();
    }

    AccumulateSize()
      : size(0)
    {}

    std::size_t size;

  };


  template<typename Container>
  struct FillIndices
    : public Dune::PDELab::TypeTree::DirectChildrenVisitor
    , public Dune::PDELab::TypeTree::DynamicTraversal
  {

    template<typename SPLFS, typename ChildLFS, typename TreePath, typename ChildIndex>
    void afterChild(SPLFS& splfs, const ChildLFS& childLFS, TreePath treePath, ChildIndex childIndex)
    {
      for(std::size_t i = 0; i < childLFS.size(); ++i, ++offset)
        local_container[offset] = childLFS.localIndex(i);
    }

    FillIndices(Container& local_container_)
      : offset(0)
      , local_container(local_container_)
    {}

    std::size_t offset;
    Container& local_container;

  };

} // anonymous namespace

template<typename GFS, typename N, typename BaseLFS, typename SubProblem_>
struct SubProblemLocalFunctionSpaceTraits
  : public LocalFunctionSpaceBaseTraits<GFS,typename GFS::Ordering::Traits::DOFIndex>
{

  //! type of local function space node
  typedef N NodeType;

  //! \brief Type to store indices from Backend
  typedef typename GFS::Traits::GridType GridType;

  typedef typename GFS::Traits::GridType Grid;

  typedef SubProblem_ SubProblem;

  typedef typename SubProblem::Traits::Condition Condition;

  typedef BaseLFS BaseLocalFunctionSpace;

};


template<typename GFS, typename N, typename BaseLFS, typename SubProblem_>
struct SubProblemLeafLocalFunctionSpaceTraits
  : public SubProblemLocalFunctionSpaceTraits<GFS,N,BaseLFS,SubProblem_>
{

  //! \brief finite element
  typedef typename BaseLFS::Traits::FiniteElement FiniteElementType;

  typedef typename BaseLFS::Traits::FiniteElement FiniteElement;

  typedef typename BaseLFS::Traits::ConstraintsType ConstraintsType;

};


namespace {

  // Test whether the argument pack values... contains t
  template<typename T, T t, T... values>
  struct arg_pack_contains_value
    : public std::false_type
  {};

  template<typename T, T t, T v, T... values>
  struct arg_pack_contains_value<T,t,v,values...>
    : public std::conditional<t == v,
                              std::true_type,
                              arg_pack_contains_value<T,t,values...>
                              >::type
  {};

  // Test whether the argument pack values... contains duplicates
  template<typename T, T... values>
  struct arg_pack_contains_duplicate_values
    : public std::false_type
  {};

  template<typename T, T v, T... values>
  struct arg_pack_contains_duplicate_values<T,v,values...>
    : public std::conditional<arg_pack_contains_value<T,v,values...>::value,
                              std::true_type,
                              arg_pack_contains_duplicate_values<T,values...>
                              >::type
  {};

} // anonymous namespace

// ********************************************************************************
// LocalFunctionSpace for subproblems - multi-component version

template<typename MDLFS, typename SubProblem, std::size_t... ChildIndices>
class SubProblemLocalFunctionSpace
  : public Dune::PDELab::TypeTree::FilteredCompositeNode<const MDLFS,Dune::PDELab::TypeTree::IndexFilter<ChildIndices...> >
  , public Dune::PDELab::LocalFunctionSpaceBaseNode<typename MDLFS::Traits::GridFunctionSpaceType,typename MDLFS::Traits::GridFunctionSpace::Ordering::Traits::DOFIndex>
{

  dune_static_assert((sizeof...(ChildIndices) <= MDLFS::CHILDREN),
                     "SubProblemLocalFunctionSpace cannot have more components than the MultiDomainGridFunctionSpace");

  dune_static_assert((!arg_pack_contains_duplicate_values<std::size_t,ChildIndices...>::value),
                     "All child indices have to be distinct");

  typedef Dune::PDELab::TypeTree::FilteredCompositeNode<const MDLFS,Dune::PDELab::TypeTree::IndexFilter<ChildIndices...> > NodeT;
  typedef Dune::PDELab::LocalFunctionSpaceBaseNode<typename MDLFS::Traits::GridFunctionSpaceType, typename MDLFS::Traits::GridFunctionSpace::Ordering::Traits::DOFIndex> BaseT;

  typedef typename MDLFS::Traits::GridFunctionSpaceType GFS;

  // friend declarations for bind() visitors
  friend struct AccumulateSize;

  template<typename>
  friend struct FillIndices;

  typedef typename GFS::Traits::Backend B;
  typedef typename GFS::Traits::GridType::Traits::template Codim<0>::Entity Element;

public:
  typedef SubProblemLocalFunctionSpaceTraits<GFS,
                                             SubProblemLocalFunctionSpace,
                                             void, // we are not directly based on a node in the original tree
                                             SubProblem> Traits;

public:

  //! \brief initialize with grid function space
  SubProblemLocalFunctionSpace (const MDLFS& mdlfs, const SubProblem& subProblem)
    : NodeT(stackobject_to_shared_ptr(mdlfs))
    , BaseT(mdlfs.gridFunctionSpaceStorage())
    , _subProblem(subProblem)
  {
    bind();
  }

  SubProblemLocalFunctionSpace (shared_ptr<const MDLFS> mdlfs, const SubProblem& subProblem)
    : NodeT(mdlfs)
    , BaseT(mdlfs.gridFunctionSpaceStorage())
    , _subProblem(subProblem)
  {
    bind();
  }

  typename Traits::IndexContainer::size_type localVectorSize() const
  {
    return mdlfs().localVectorSize();
  }

  typename Traits::IndexContainer::size_type maxSize() const
  {
    return mdlfs().maxSize();
  }

  typename Traits::IndexContainer::size_type localIndex (typename Traits::IndexContainer::size_type index) const
  {
    return _local_storage[index];
  }

  typename Traits::DOFIndex& dofIndex(typename Traits::IndexContainer::size_type index) const
  {
    return mdlfs().dofIndex(localIndex(index));
  }

  //! \brief bind local function space to entity
  void bind ()
  {
    // make offset
    offset = 0;
    AccumulateSize accumulateSize;
    Dune::PDELab::TypeTree::applyToTree(*this,accumulateSize);
    n = accumulateSize.size;
    _local_storage.resize(n);
    Dune::PDELab::TypeTree::applyToTree(*this,FillIndices<typename Traits::IndexContainer>(_local_storage));
  }

  const SubProblem& subProblem() const {
    return _subProblem;
  }

  template<typename EG>
  bool appliesTo(const EG& eg) const {
    return _subProblem.appliesTo(eg);
  }

private:

  using BaseT::offset;
  using BaseT::n;
  using NodeT::unfiltered;

  typename Traits::IndexContainer _local_storage;

  const MDLFS mdlfs() const
  {
    return unfiltered();
  }

  const SubProblem& _subProblem;

};


// common base class for single-component versions to avoid code duplication


template<typename Traits_>
class SubProblemLocalFunctionSpaceProxy
{

protected:

  typedef typename Traits_::BaseLocalFunctionSpace BaseLFS;

  const BaseLFS& baseLFS() const {
    return static_cast<const typename Traits_::NodeType&>(*this).proxiedNode();
  }

public:

  typedef Traits_ Traits;
  typedef typename Traits::SubProblem SubProblem;

  //! \brief initialize with grid function space
  SubProblemLocalFunctionSpaceProxy (const SubProblem& subProblem)
    : _subProblem(subProblem)
  {
  }

  //! \brief get current size
  typename Traits::IndexContainer::size_type size () const
  {
    return baseLFS().size();
  }

  //! \brief get maximum possible size (which is maxLocalSize from grid function space)
  typename Traits::IndexContainer::size_type maxSize () const
  {
    return baseLFS().maxLocalSize();
  }

  typename Traits::IndexContainer::size_type localVectorSize() const
  {
    return baseLFS().localVectorSize();
  }

  // map index in this local function space to root local function space
  typename Traits::IndexContainer::size_type localIndex (typename Traits::IndexContainer::size_type index) const
  {
    return baseLFS().localIndex(index);
  }

  // map index in this local function space to global index space
  const typename Traits::DOFIndex& dofIndex (typename Traits::IndexContainer::size_type index) const
  {
    return baseLFS().dofIndex(index);
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

  template<typename GC, typename LC>
  void mwrite (const LC& lc, GC& gc) const
  {
    baseLFS().mwrite(lc,gc);
  }

  const SubProblem& subProblem() const {
    return _subProblem;
  }

  template<typename EG>
  bool appliesTo(const EG& eg) const {
    return _subProblem.appliesTo(eg);
  }

  // This method does nothing - it only exists for compatibility with the multi-child version
  void bind()
  {
  }

private:
  const SubProblem& _subProblem;

};


// single-component version base - this needs to be specialized for each supported base LFS


template<typename MDLFS, typename LFS, typename BaseLFS, typename SubProblem, typename LFSTag>
class SubProblemLocalFunctionSpaceBase
  : public SubProblemLocalFunctionSpaceProxy<SubProblemLocalFunctionSpaceTraits<
                                               typename MDLFS::Traits::GridFunctionSpace,
                                               LFS,
                                               BaseLFS,
                                               SubProblem
                                               >
                                             >
{

  typedef SubProblemLocalFunctionSpaceProxy<SubProblemLocalFunctionSpaceTraits<
                                              typename MDLFS::Traits::GridFunctionSpace,
                                              LFS,
                                              BaseLFS,
                                              SubProblem
                                              >
                                            > BaseT;

protected:

  SubProblemLocalFunctionSpaceBase(const SubProblem& subProblem)
    : BaseT(subProblem)
  {}

};


// single-component version base - specialization for leaf function space

template<typename MDLFS, typename LFS, typename BaseLFS, typename SubProblem>
class SubProblemLocalFunctionSpaceBase<MDLFS,LFS,BaseLFS,SubProblem,LeafLocalFunctionSpaceTag>
  : public SubProblemLocalFunctionSpaceProxy<SubProblemLeafLocalFunctionSpaceTraits<
                                               typename MDLFS::Traits::GridFunctionSpaceType,
                                               LFS,
                                               BaseLFS,
                                               SubProblem
                                               >
                                             >
{

public:

  typedef SubProblemLocalFunctionSpaceProxy<SubProblemLeafLocalFunctionSpaceTraits<
                                              typename MDLFS::Traits::GridFunctionSpaceType,
                                              LFS,
                                              BaseLFS,
                                              SubProblem
                                              >
                                            > BaseT;

  typedef typename BaseT::Traits Traits;

  const typename Traits::FiniteElementType& finiteElement() const {
    return this->baseLFS().finiteElement();
  }

  const typename Traits::ConstraintsType& constraints() const {
    return this->baseLFS().constraints();
  }

  template<typename GC, typename LC>
  void insert_constraints (const LC& lc, GC& gc) const
  {
    // LC and GC are maps of maps
    typedef typename LC::const_iterator local_col_iterator;
    typedef typename LC::value_type::second_type::const_iterator local_row_iterator;
    typedef typename GC::iterator global_col_iterator;
    typedef typename GC::value_type::second_type global_row_type;

    for (local_col_iterator cit=lc.begin(); cit!=lc.end(); ++cit)
      {

        // look up entry in global map, if not found, insert an empty one.
        global_col_iterator gcit = gc.insert(std::make_pair(std::ref(this->dofIndex(cit->first)),global_row_type())).first;

        // copy row to global container with transformed indices
        for (local_row_iterator rit=(cit->second).begin(); rit!=(cit->second).end(); ++rit)
          gcit->second[this->dofIndex(rit->first)] = rit->second;
      }
  }

protected:

  SubProblemLocalFunctionSpaceBase(const SubProblem& subProblem)
    : BaseT(subProblem)
  {}

};




template<typename MDLFS, typename SubProblem, std::size_t ChildIndex>
class SubProblemLocalFunctionSpace<MDLFS,SubProblem,ChildIndex>
  : public Dune::PDELab::TypeTree::ProxyNode<const typename MDLFS::template Child<ChildIndex>::Type>
  , public SubProblemLocalFunctionSpaceBase<MDLFS,
                                            SubProblemLocalFunctionSpace<MDLFS,SubProblem,ChildIndex>,
                                            typename MDLFS::template Child<ChildIndex>::Type,
                                            SubProblem,
                                            typename MDLFS::template Child<ChildIndex>::Type::ImplementationTag
                                            >
{

  typedef Dune::PDELab::TypeTree::ProxyNode<const typename MDLFS::template Child<ChildIndex>::Type> NodeT;

  typedef SubProblemLocalFunctionSpaceBase<MDLFS,
                                           SubProblemLocalFunctionSpace<MDLFS,SubProblem,ChildIndex>,
                                           typename MDLFS::template Child<ChildIndex>::Type,
                                           SubProblem,
                                           typename MDLFS::template Child<ChildIndex>::Type::ImplementationTag
                                           > BaseT;

public:

  SubProblemLocalFunctionSpace (const MDLFS& mdlfs, const SubProblem& subProblem)
    : NodeT(mdlfs.template childStorage<ChildIndex>())
    , BaseT(subProblem)
  {
  }

  using NodeT::proxiedNode;

};


    //! \} group GridFunctionSpace

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_SUBPROBLEMLOCALFUNCTIONSPACE_HH
