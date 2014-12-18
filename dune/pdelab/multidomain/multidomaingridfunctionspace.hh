#ifndef DUNE_MULTIDOMAIN_MULTIDOMAINGRIDFUNCTIONSPACE_HH
#define DUNE_MULTIDOMAIN_MULTIDOMAINGRIDFUNCTIONSPACE_HH

#include <tuple>
#include <type_traits>

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powercompositegridfunctionspacebase.hh>
#include <dune/pdelab/ordering/transformations.hh>
#include <dune/pdelab/multidomain/multidomainlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/typemap.hh>
#include <dune/grid/multidomaingrid.hh>
#include <dune/pdelab/multidomain/utility.hh>
#include <dune/pdelab/multidomain/couplinggridfunctionspace.hh>
#include <dune/pdelab/multidomain/datahandleprovider.hh>
#include <dune/pdelab/multidomain/powercouplinggridfunctionspace.hh>
#include <utility>

namespace Dune {

namespace PDELab {

namespace MultiDomain {

template<typename G, typename B, typename M, bool supportMapAccess, int k>
struct MultiDomainGridFunctionSpaceTraits
{
  enum{
    //! \brief True if this grid function space is composed of others.
    isComposite = 1,
    //! \brief number of child spaces
    noChilds = k
  };

  //! \brief the grid view where grid function is defined upon
  typedef G GridType;

  typedef typename G::LeafGridView GridViewType;


  typedef typename G::LeafGridView GridView;

  //! \brief vector backend
  typedef B Backend;

  //! \brief mapper
  typedef M MapperType;

  typedef M OrderingTag;

  //! \brief short cut for size type exported by Backend
  typedef typename B::size_type SizeType;

  static const bool supportsMapAccess = supportMapAccess;

};

struct VerifyChildren
  : public TypeTree::DirectChildrenVisitor
  , public TypeTree::DynamicTraversal
{

  template<typename T, typename Child, typename TreePath, typename ChildIndex>
  void beforeChild(const T& t, const Child& child, TreePath treePath, ChildIndex childIndex) const
  {
    typedef typename T::Traits::GridType MultiDomainGrid;
    typedef typename MultiDomainGrid::SubDomainGrid SubDomainGrid;
    typedef typename Child::Traits::GridViewType::Grid ChildGrid;
    static_assert((is_same<MultiDomainGrid,ChildGrid>::value || is_same<SubDomainGrid,ChildGrid>::value),
                  "MultiDomainGridFunctionSpace only works with a MultiDomainGrid and its associated SubDomainGrids.");
    doVerify<T>(t.grid(),child.gridView().grid());
  }

  template<typename T>
  static void doVerify(const typename T::Traits::GridType& g, const typename T::Traits::GridType& cg)
  {
    assert(&g == &cg);
  }

  template<typename T>
  static void doVerify(const typename T::Traits::GridType& g, const typename T::Traits::GridType::SubDomainGrid& cg)
  {
    assert(&g == &cg.multiDomainGrid());
  }

  template<typename T, typename ST>
  static void adoVerify(const typename T::Traits::GridType& g, const ST& cg)
  {
    // this is only here to keep the compiler from complaining about a missing function
    // if ST is not a Multi-/SubDomainGrid
  }

};


template<typename GFS>
struct gfs_flavor_tag
{

  static const Dune::mdgrid::MultiDomainGridType gridType = Dune::mdgrid::GridType<typename GFS::Traits::GridViewType::Grid>::v;

  typedef typename conditional<gridType == Dune::mdgrid::multiDomainGrid,
                               MultiDomainGFSTag,
                               SubDomainGFSTag
                               >::type type;

};

template<typename GV, typename LFEM, typename Predicate_, typename CE, typename B>
struct gfs_flavor_tag<CouplingGridFunctionSpace<GV,LFEM,Predicate_,CE,B> >
{
  typedef CouplingGFSTag type;
};

template<typename T, std::size_t k, typename Backend, typename Ordering>
struct gfs_flavor_tag<PowerCouplingGridFunctionSpace<T,k,Backend,Ordering> >
{
  typedef CouplingGFSTag type;
};

template<typename G, typename Backend, typename OrderingTag, typename... Children>
class MultiDomainGridFunctionSpace
  : public TypeTree::CompositeNode<Children...>
  , public PowerCompositeGridFunctionSpaceBase<MultiDomainGridFunctionSpace<G,Backend,OrderingTag,Children...>,
                                               typename TypeTree::CompositeNode<Children...>::template Child<0>::Type::Traits::GridViewType,
                                               Backend,
                                               OrderingTag,
                                               sizeof...(Children)
                                               >
  , public DataHandleProvider<MultiDomainGridFunctionSpace<G,Backend,OrderingTag,Children...> >
{

  typedef TypeTree::CompositeNode<Children...> BaseT;

  friend class
  PowerCompositeGridFunctionSpaceBase<MultiDomainGridFunctionSpace<G,Backend,OrderingTag,Children...>,
                                      typename TypeTree::CompositeNode<Children...>::template Child<0>::Type::Traits::GridViewType,
                                      Backend,
                                      OrderingTag,
                                      sizeof...(Children)
                                      >;

  typedef PowerCompositeGridFunctionSpaceBase<MultiDomainGridFunctionSpace<G,Backend,OrderingTag,Children...>,
                                              typename TypeTree::CompositeNode<Children...>::template Child<0>::Type::Traits::GridViewType,
                                              Backend,
                                              OrderingTag,
                                              sizeof...(Children)
                                              > ImplementationBase;

  typedef TypeTree::TransformTree<MultiDomainGridFunctionSpace,gfs_to_ordering<MultiDomainGridFunctionSpace> > ordering_transformation;

  template<typename,typename>
  friend class Dune::PDELab::GridFunctionSpaceBase;

public:

  typedef MultiDomainGridFunctionSpaceTag ImplementationTag;

  template<typename T, std::size_t i, typename Tag_>
  struct GFSChild
  {
    static const std::size_t index = i;
    typedef T type;
    typedef T key;
    typedef Tag_ Tag;
    static const bool definedOnSubDomain = std::is_same<Tag,SubDomainGFSTag>::value;
    static const bool isCouplingSpace = std::is_same<Tag,CouplingGFSTag>::value;
    static const bool isNormalSpace = std::is_same<Tag,MultiDomainGFSTag>::value || std::is_same<Tag,SubDomainGFSTag>::value;
  };

private:

  template<typename Grid, template<typename...> class Container>
  struct tagger
  {

    template<std::size_t i, typename T>
    struct transform {
      typedef GFSChild<T,
                       i,
                       typename gfs_flavor_tag<T>::type
                       > type;
    };

    template<typename... Args>
    struct container {
      typedef Container<Args...> type;
    };
  };

  typedef typename indexed_transform<tagger<G,std::tuple>,Children...>::type ChildEntries;

  typedef typename indexed_transform<tagger<G,embedded_key_map>,Children...>::type ChildEntryMap;

public:

  template<std::size_t i>
  struct ChildInfo
  {
    static_assert(i < BaseT::CHILDREN,"invalid child index");

    typedef typename std::tuple_element<i,ChildEntries>::type Type;
  };

  //! export traits class
  typedef MultiDomainGridFunctionSpaceTraits<G,
                                             Backend,
                                             void,
                                             !ChildEntryMap::has_duplicate_entries, // duplicate entries would lead to ambiguous type lookups
                                             sizeof...(Children)>
  Traits;

  typedef typename ordering_transformation::Type Ordering;

  template<typename T>
  struct IndexForChild
  {
    static_assert(Traits::supportsMapAccess,"This MultiDomainGridFunctionSpace does not support accessing children by type - " \
                                            "it probably contains two child grid function spaces with identical type. " \
                                            "Select your function spaces by index instead.");
    static const std::size_t value = get_map_entry<T,ChildEntryMap>::type::index;
  };

  template<typename T>
  const T& childByType() const
  {
    static_assert(Traits::supportsMapAccess,"This MultiDomainGridFunctionSpace does not support accessing children by type - " \
                                            "it probably contains two child grid function spaces with identical type. " \
                                            "Select your function spaces by index instead.");
    return this->template getChild<get_map_entry<T,ChildEntryMap>::type::index>();
  }

  template<typename T>
  T& childByType()
  {
    static_assert(Traits::supportsMapAccess,"This MultiDomainGridFunctionSpace does not support accessing children by type - " \
                                            "it probably contains two child grid function spaces with identical type. " \
                                            "Select your function spaces by index instead.");
    return this->template getChild<get_map_entry<T,ChildEntryMap>::type::index>();
  }


  typename G::LeafGridView gridView() const
  {
    return grid().leafGridView();
  }

  typename G::LeafGridView gridview() const DUNE_DEPRECATED_MSG("Use gridView() instead of gridview()")
  {
    return grid().leafGridView();
  }

  MultiDomainGridFunctionSpace (G& g, const Backend& backend, const OrderingTag& ordering_tag, Children&... children)
    : BaseT(children...)
    , ImplementationBase(backend,ordering_tag)
    , _g(g)
  {
    static_assert(Dune::mdgrid::GridType<G>::v == Dune::mdgrid::multiDomainGrid,
                  "MultiDomainGridFunctionSpace only works on a MultiDomainGrid");
    TypeTree::applyToTree(*this,VerifyChildren());
  }

  MultiDomainGridFunctionSpace (G& g, const Backend& backend, Children&... children)
    : BaseT(children...)
    , ImplementationBase(backend,OrderingTag())
    , _g(g)
  {
    static_assert(Dune::mdgrid::GridType<G>::v == Dune::mdgrid::multiDomainGrid,
                  "MultiDomainGridFunctionSpace only works on a MultiDomainGrid");
    TypeTree::applyToTree(*this,VerifyChildren());
  }

  MultiDomainGridFunctionSpace (G& g, const OrderingTag& ordering_tag, Children&... children)
    : BaseT(children...)
    , ImplementationBase(Backend(),ordering_tag)
    , _g(g)
  {
    static_assert(Dune::mdgrid::GridType<G>::v == Dune::mdgrid::multiDomainGrid,
                  "MultiDomainGridFunctionSpace only works on a MultiDomainGrid");
    TypeTree::applyToTree(*this,VerifyChildren());
  }

  MultiDomainGridFunctionSpace (G& g, Children&... children)
    : BaseT(children...)
    , ImplementationBase(Backend(),OrderingTag())
    , _g(g)
  {
    static_assert(Dune::mdgrid::GridType<G>::v == Dune::mdgrid::multiDomainGrid,
                  "MultiDomainGridFunctionSpace only works on a MultiDomainGrid");
    TypeTree::applyToTree(*this,VerifyChildren());
  }


  //------------------------------

  G& grid() {
    return _g;
  }

  const G& grid() const {
    return _g;
  }


  //! Direct access to the DOF ordering.
  const Ordering &ordering() const
  {
    if (!this->isRootSpace())
      {
        DUNE_THROW(GridFunctionSpaceHierarchyError,
                   "Ordering can only be obtained for root space in GridFunctionSpace tree.");
      }
    if (!_ordering)
      {
        create_ordering();
        this->update(*_ordering);
      }
    return *_ordering;
  }

  //! Direct access to the DOF ordering.
  Ordering &ordering()
  {
    if (!this->isRootSpace())
      {
        DUNE_THROW(GridFunctionSpaceHierarchyError,
                   "Ordering can only be obtained for root space in GridFunctionSpace tree.");
      }
    if (!_ordering)
      {
        create_ordering();
        this->update(*_ordering);
      }
    return *_ordering;
  }

  //! Direct access to the storage of the DOF ordering.
  shared_ptr<const Ordering> orderingStorage() const
  {
    if (!this->isRootSpace())
      {
        DUNE_THROW(GridFunctionSpaceHierarchyError,
                   "Ordering can only be obtained for root space in GridFunctionSpace tree.");
      }
    if (!_ordering)
      {
        create_ordering();
        this->update(*_ordering);
      }
    return _ordering;
  }

  //! Direct access to the storage of the DOF ordering.
  shared_ptr<Ordering> orderingStorage()
  {
    if (!this->isRootSpace())
      {
        DUNE_THROW(GridFunctionSpaceHierarchyError,
                   "Ordering can only be obtained for root space in GridFunctionSpace tree.");
      }
    if (!_ordering)
      {
        create_ordering();
        this->update(*_ordering);
      }
    return _ordering;
  }


private:

  // This method here is to avoid a double update of the Ordering when the user calls
  // GFS::update() before GFS::ordering().
  void create_ordering() const
  {
    _ordering = make_shared<Ordering>(ordering_transformation::transform(*this));
  }

  const G& _g;
  mutable shared_ptr<Ordering> _ordering;

};

template<typename MultiGFS, typename ChildGFS>
struct get_subproblem_index
{
  static const std::size_t value = MultiGFS::template IndexForChild<ChildGFS>::value;
};


// ********************************************************************************
// register GFS -> (Local)Ordering transformations
// ********************************************************************************

template<typename GridFunctionSpace, typename Params>
composite_gfs_to_ordering_descriptor<
  GridFunctionSpace,
  gfs_to_ordering<Params>,
  typename GridFunctionSpace::OrderingTag
  >
registerNodeTransformation(GridFunctionSpace* gfs, gfs_to_ordering<Params>* t, MultiDomainGridFunctionSpaceTag* tag);

template<typename GridFunctionSpace, typename Params>
composite_gfs_to_local_ordering_descriptor<
  GridFunctionSpace,
  gfs_to_local_ordering<Params>,
  typename GridFunctionSpace::OrderingTag
  >
registerNodeTransformation(GridFunctionSpace* gfs, gfs_to_local_ordering<Params>* t, MultiDomainGridFunctionSpaceTag* tag);


  /*
template<typename MultiGFS, typename ChildGFS>
class TypeBasedGridFunctionSubSpace
  : public GridFunctionSubSpace<MultiGFS,get_subproblem_index<MultiGFS,ChildGFS>::value>
{

  typedef GridFunctionSubSpace<MultiGFS,get_subproblem_index<MultiGFS,ChildGFS>::value> BaseT;

public:
  TypeBasedGridFunctionSubSpace(const MultiGFS& gfs)
    : BaseT(gfs)
  {}
};
  */

} // namespace MultiDomain

} // namespace PDELab

} // namespace Dune

#endif // DUNE_MULTIDOMAIN_MULTIDOMAINGRIDFUNCTIONSPACE_HH
