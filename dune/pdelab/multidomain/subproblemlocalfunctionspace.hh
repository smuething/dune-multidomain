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
        container[offset] = childLFS.globalIndex(i);
    }

    FillIndices(Container& container_)
      : offset(0)
      , container(container_)
    {}

    std::size_t offset;
    Container& container;

  };

} // anonymous namespace

template<typename GFS, typename N, typename BaseLFS, typename SubProblem_, typename Constraints_>
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

  typedef SubProblem_ SubProblem;

  typedef typename SubProblem::Traits::Condition Condition;

  typedef Constraints_ Constraints;

  typedef Constraints ConstraintsType;

  typedef BaseLFS BaseLocalFunctionSpace;

};


template<typename GFS, typename N, typename BaseLFS, typename SubProblem_, typename Constraints_>
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

  //! \brief finite element
  typedef typename BaseLFS::Traits::FiniteElementType FiniteElementType;

  typedef SubProblem_ SubProblem;

  typedef typename SubProblem::Traits::Condition Condition;

  typedef Constraints_ Constraints;

  typedef Constraints ConstraintsType;

  typedef BaseLFS BaseLocalFunctionSpace;

};


namespace {

  // Test whether the argument pack values... contains t
  template<typename T, T t, T... values>
  struct arg_pack_contains_value
  {
    static const bool value = false;
  };

  template<typename T, T t, T v, T... values>
  struct arg_pack_contains_value<T,t,v,values...>
  {
    static const bool value = (t == v) || arg_pack_contains_value<T,t,values...>::value;
  };

  // Test whether the argument pack values... contains duplicates
  template<typename T, T... values>
  struct arg_pack_contains_duplicate_values
  {
    static const bool value = false;
  };

  template<typename T, T v, T... values>
  struct arg_pack_contains_duplicate_values<T,v,values...>
  {
    static const bool value = arg_pack_contains_value<T,v,values...>::value ||
      arg_pack_contains_duplicate_values<T,values...>::value;
  };

} // anonymous namespace

// ********************************************************************************
// LocalFunctionSpace for subproblems - multi-component version

template<typename MDLFS, typename SubProblem, typename Constraints, std::size_t... ChildIndices>
class SubProblemLocalFunctionSpace
  : public Dune::PDELab::TypeTree::FilteredCompositeNode<const MDLFS,Dune::PDELab::TypeTree::IndexFilter<ChildIndices...> >
  , public Dune::PDELab::LocalFunctionSpaceBaseNode<typename MDLFS::Traits::GridFunctionSpaceType>
{

  dune_static_assert((sizeof...(ChildIndices) <= MDLFS::CHILDREN),
                     "SubProblemLocalFunctionSpace cannot have more components than the MultiDomainGridFunctionSpace");

  dune_static_assert((!arg_pack_contains_duplicate_values<std::size_t,ChildIndices...>::value),
                     "All child indices have to be distinct");

  typedef Dune::PDELab::TypeTree::FilteredCompositeNode<const MDLFS,Dune::PDELab::TypeTree::IndexFilter<ChildIndices...> > NodeT;
  typedef Dune::PDELab::LocalFunctionSpaceBaseNode<typename MDLFS::Traits::GridFunctionSpaceType> BaseT;

  typedef typename MDLFS::Traits::GridFunctionSpaceType GFS;

  // friend declarations for bind() visitors
  friend struct AccumulateSize;
  friend struct FillIndices;

  typedef typename GFS::Traits::BackendType B;
  typedef typename GFS::Traits::GridType::Traits::template Codim<0>::Entity Element;

public:
  typedef SubProblemLocalFunctionSpaceTraits<GFS,
                                             SubProblemLocalFunctionSpace,
                                             void, // we are not directly based on a node in the original tree
                                             SubProblem,
                                             Constraints> Traits;

public:

  //! \brief empty constructor - TODO:Do we need this???
  /*SubProblemLocalFunctionSpace ()
  {
  }*/

  //! \brief initialize with grid function space
  SubProblemLocalFunctionSpace (const MDLFS& mdlfs, const SubProblem& subProblem, const Constraints& constraints)
    : NodeT(stackobject_to_shared_ptr(mdlfs))
    , BaseT(mdlfs.gridFunctionSpaceStorage())
    , _subProblem(subProblem)
    , _constraints(constraints)
  {
    bind();
  }

  SubProblemLocalFunctionSpace (shared_ptr<const MDLFS> mdlfs, const SubProblem& subProblem, const Constraints& constraints)
    : NodeT(mdlfs)
    , BaseT(mdlfs.gridFunctionSpaceStorage())
    , _subProblem(subProblem)
    , _constraints(constraints)
  {
    bind();
  }

  typename Traits::IndexContainer::size_type localVectorSize() const
  {
    return mdlfs().localVectorSize();
  }

  //! This method does not work for SubProblemLocalFunctionSpace.
  template<typename T>
  typename Traits::IndexContainer::size_type localIndex (T index) const
  {
    dune_static_assert((AlwaysFalse<T>::value),"SubProblemLocalFunctionSpace::localIndex() is not implemented," \
                       "use the localIndex() methods of its children instead.");
  }

  //! \brief bind local function space to entity
  void bind ()
  {
    // make offset
    offset = 0;
    AccumulateSize accumulateSize;
    Dune::PDELab::TypeTree::applyToTree(*this,accumulateSize);
    n = accumulateSize.size;
    global_storage.resize(n);
    Dune::PDELab::TypeTree::applyToTree(*this,FillIndices<typename Traits::IndexContainer>(global_storage));
  }

  const SubProblem& subProblem() const {
    return _subProblem;
  }

  const Constraints& constraints() const {
    return _constraints;
  }

  template<typename EG>
  bool appliesTo(const EG& eg) const {
    return _subProblem.appliesTo(eg);
  }

private:

  using BaseT::offset;
  using BaseT::global;
  using BaseT::global_storage;
  using BaseT::n;
  using NodeT::unfiltered;


  const MDLFS mdlfs() const
  {
    return unfiltered();
  }

  const SubProblem& _subProblem;
  const Constraints& _constraints;

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
  typedef typename Traits::Constraints Constraints;

  //! \brief initialize with grid function space
  SubProblemLocalFunctionSpaceProxy (const SubProblem& subProblem, const Constraints& constraints)
    : _subProblem(subProblem)
    , _constraints(constraints)
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

  template<typename GC, typename LC>
  void mwrite (const LC& lc, GC& gc) const
  {
    baseLFS().mwrite(lc,gc);
  }

  const SubProblem& subProblem() const {
    return _subProblem;
  }

  const Constraints& constraints() const {
    return _constraints;
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
  const Constraints& _constraints;

};


// single-component version base - this needs to be specialized for each supported base LFS


  template<typename MDLFS, typename LFS, typename BaseLFS, typename SubProblem, typename Constraints, typename LFSTag>
class SubProblemLocalFunctionSpaceBase
  : public SubProblemLocalFunctionSpaceProxy<SubProblemLocalFunctionSpaceTraits<
                                               typename MDLFS::Traits::GridFunctionSpaceType,
                                               LFS,
                                               BaseLFS,
                                               SubProblem,
                                               Constraints
                                               >
                                             >
{

  typedef SubProblemLocalFunctionSpaceProxy<SubProblemLocalFunctionSpaceTraits<
                                              typename MDLFS::Traits::GridFunctionSpaceType,
                                              LFS,
                                              BaseLFS,
                                              SubProblem,
                                              Constraints
                                              >
                                            > BaseT;

protected:

  SubProblemLocalFunctionSpaceBase(const SubProblem& subProblem, const Constraints& constraints)
    : BaseT(subProblem,constraints)
  {}

};


// single-component version base - specialization for leaf function space

  template<typename MDLFS, typename LFS, typename BaseLFS, typename SubProblem, typename Constraints>
  class SubProblemLocalFunctionSpaceBase<MDLFS,LFS,BaseLFS,SubProblem,Constraints,LeafLocalFunctionSpaceTag>
  : public SubProblemLocalFunctionSpaceProxy<SubProblemLeafLocalFunctionSpaceTraits<
                                               typename MDLFS::Traits::GridFunctionSpaceType,
                                               LFS,
                                               BaseLFS,
                                               SubProblem,
                                               Constraints
                                               >
                                             >
{

public:

  typedef SubProblemLocalFunctionSpaceProxy<SubProblemLeafLocalFunctionSpaceTraits<
                                              typename MDLFS::Traits::GridFunctionSpaceType,
                                              LFS,
                                              BaseLFS,
                                               SubProblem,
                                              Constraints
                                              >
                                            > BaseT;

  typedef typename BaseT::Traits Traits;

  const typename Traits::FiniteElementType& finiteElement() const {
    return this->baseLFS().finiteElement();
  }

protected:

  SubProblemLocalFunctionSpaceBase(const SubProblem& subProblem, const Constraints& constraints)
    : BaseT(subProblem,constraints)
  {}

};




template<typename MDLFS, typename SubProblem, typename Constraints, std::size_t ChildIndex>
class SubProblemLocalFunctionSpace<MDLFS,SubProblem,Constraints,ChildIndex>
  : public Dune::PDELab::TypeTree::ProxyNode<const typename MDLFS::template Child<ChildIndex>::Type>
  , public SubProblemLocalFunctionSpaceBase<MDLFS,
                                            SubProblemLocalFunctionSpace<MDLFS,SubProblem,Constraints,ChildIndex>,
                                            typename MDLFS::template Child<ChildIndex>::Type,
                                            SubProblem,
                                            Constraints,
                                            typename MDLFS::template Child<ChildIndex>::Type::ImplementationTag
                                            >
{

  typedef Dune::PDELab::TypeTree::ProxyNode<const typename MDLFS::template Child<ChildIndex>::Type> NodeT;

  typedef SubProblemLocalFunctionSpaceBase<MDLFS,
                                           SubProblemLocalFunctionSpace<MDLFS,SubProblem,Constraints,ChildIndex>,
                                           typename MDLFS::template Child<ChildIndex>::Type,
                                           SubProblem,
                                           Constraints,
                                           typename MDLFS::template Child<ChildIndex>::Type::ImplementationTag
                                           > BaseT;

public:

  SubProblemLocalFunctionSpace (const MDLFS& mdlfs, const SubProblem& subProblem, const Constraints& constraints)
    : NodeT(mdlfs.template childStorage<ChildIndex>())
    , BaseT(subProblem,constraints)
  {
  }

};


    //! \} group GridFunctionSpace

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_SUBPROBLEMLOCALFUNCTIONSPACE_HH
