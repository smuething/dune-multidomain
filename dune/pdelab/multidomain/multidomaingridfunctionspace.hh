#ifndef DUNE_MULTIDOMAIN_MULTIDOMAINGRIDFUNCTIONSPACE_HH
#define DUNE_MULTIDOMAIN_MULTIDOMAINGRIDFUNCTIONSPACE_HH

#include <tuple>
#include <type_traits>
#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powercompositegridfunctionspacebase.hh>
#include <dune/pdelab/multidomain/multidomainlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/typemap.hh>
#include <dune/grid/multidomaingrid.hh>
#include <dune/pdelab/multidomain/utility.hh>
#include <dune/pdelab/multidomain/couplinggridfunctionspace.hh>
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

  //! \brief vector backend
  typedef B BackendType;

  //! \brief mapper
  typedef M MapperType;

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
    dune_static_assert((is_same<MultiDomainGrid,ChildGrid>::value || is_same<SubDomainGrid,ChildGrid>::value),
                       "MultiDomainGridFunctionSpace only works with a MultiDomainGrid and its associated SubDomainGrids.");
    doVerify<T>(t.grid(),t.template getChild<i>().gridview().grid());
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
  static void doVerify(const typename T::Traits::GridType& g, const ST& cg)
  {
    // this is only here to keep the compiler from complaining about a missing function
    // if ST is not a Multi-/SubDomainGrid
  }

};


template<typename GFS>
struct gfs_flavor_tag
{

  static const Dune::mdgrid::MultiDomainGridType gridType = Dune::mdgrid::GridType<typename GFS::Traits::GridViewType::Grid>::v;

  typedef typename Dune::SelectType<gridType == Dune::mdgrid::multiDomainGrid,
                                    MultiDomainGFSTag,
                                    SubDomainGFSTag
                                    >::Type type;

};

template<typename GV, typename LFEM, typename Predicate_, typename CE, typename B>
struct gfs_flavor_tag<CouplingGridFunctionSpace<GV,LFEM,Predicate_,CE,B> >
{
  typedef CouplingGFSTag type;
};


template<typename G, typename... Children>
class MultiDomainGridFunctionSpace
  : public TypeTree::VariadicCompositeNode<Children...>
  , public PowerCompositeGridFunctionSpaceBase<MultiDomainGridFunctionSpace<G,Children...>,
                                               typename TypeTree::VariadicCompositeNode<Children...>::template Child<0>::Type::Traits::GridViewType,
                                               typename TypeTree::VariadicCompositeNode<Children...>::template Child<0>::Type::Traits::BackendType,
                                               GridFunctionSpaceLexicographicMapper,
                                               sizeof...(Children)
                                               >
{

  typedef TypeTree::VariadicCompositeNode<Children...> BaseT;

  friend class
  PowerCompositeGridFunctionSpaceBase<MultiDomainGridFunctionSpace<Children...>,
                                      typename TypeTree::VariadicCompositeNode<Children...>::template Child<0>::Type::Traits::GridViewType,
                                      typename TypeTree::VariadicCompositeNode<Children...>::template Child<0>::Type::Traits::BackendType,
                                      GridFunctionSpaceLexicographicMapper,
                                      sizeof...(Children)
                                      >;

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
    dune_static_assert(i < BaseT::CHILDREN,"invalid child index");

    typedef typename std::tuple_element<i,ChildEntries>::type Type;
  };

  //! export traits class
  typedef MultiDomainGridFunctionSpaceTraits<G,
                                             typename BaseT::template Child<0>::Type::Traits::BackendType,
                                             CopyStoragePolicy,
                                             !ChildEntryMap::has_duplicate_entries, // duplicate entries would lead to ambiguous type lookups
                                             sizeof...(Children)>
  Traits;

  template<typename T>
  struct IndexForChild
  {
    dune_static_assert(Traits::supportsMapAccess,"This MultiDomainGridFunctionSpace does not support accessing children by type - " \
                                                 "it probably contains two child grid function spaces with identical type. " \
                                                 "Select your function spaces by index instead.");
    static const std::size_t value = get_map_entry<T,ChildEntryMap>::type::index;
  };

  template<typename T>
  const T& childByType() const
  {
    dune_static_assert(Traits::supportsMapAccess,"This MultiDomainGridFunctionSpace does not support accessing children by type - " \
                                                 "it probably contains two child grid function spaces with identical type. " \
                                                 "Select your function spaces by index instead.");
    return this->template getChild<get_map_entry<T,ChildEntryMap>::type::index>();
  }

  template<typename T>
  T& childByType()
  {
    dune_static_assert(Traits::supportsMapAccess,"This MultiDomainGridFunctionSpace does not support accessing children by type - " \
                                                 "it probably contains two child grid function spaces with identical type. " \
                                                 "Select your function spaces by index instead.");
    return this->template getChild<get_map_entry<T,ChildEntryMap>::type::index>();
  }

  typename G::LeafGridView gridview() const { return grid().leafView(); }

  // define local function space parametrized by self
  typedef typename Dune::PDELab::TypeTree::TransformTree<MultiDomainGridFunctionSpace,Dune::PDELab::gfs_to_lfs>::Type LocalFunctionSpace;
  //typedef MultiDomainCouplingLocalFunctionSpace<MultiDomainGridFunctionSpace,Children...> CouplingLocalFunctionSpace;


  MultiDomainGridFunctionSpace (G& g, Children&... children) : BaseT(children...), _g(g)
  {
    dune_static_assert(Dune::mdgrid::GridType<G>::v == Dune::mdgrid::multiDomainGrid,
                       "MultiDomainGridFunctionSpace only works on a MultiDomainGrid");
    TypeTree::applyToTree(*this,VerifyChildren());
    this->setup();
  }

  //! map index from our index set [0,size()-1] to root index set
  typename Traits::SizeType upMap (typename Traits::SizeType i) const
  {
    return i;
  }

  //! map index from child i's index set into our index set
  typename Traits::SizeType subMap (typename Traits::SizeType i, typename Traits::SizeType j) const
  {
    return offset[i]+j;
  }

  template<class EntityType>
  size_t dataHandleSize (const EntityType& e) const
  {
    // FIXME!
    //return VisitChildTMP::dataHandleSize(*this,e);
    return 0;
  }

  //! return vector of global indices associated with the given entity
  template<class EntityType>
  void dataHandleGlobalIndices (const EntityType& e,
                                std::vector<typename Traits::SizeType>& global) const
  {
    // FIXME!
    // size_t n=dataHandleSize(e);
    // global.resize(n);
    // VisitChildTMP::dataHandleGlobalIndices(*this,e,global,0,childglobal);
  }

  //------------------------------

  G& grid() {
    return _g;
  }

  const G& grid() const {
    return _g;
  }

private:

  const G& _g;

  using ImplementationBase::childLocalSize;
  using ImplementationBase::childGlobalSize;
  using ImplementationBase::maxlocalsize;
  using ImplementationBase::offset;

  void calculateSizes ()
  {
    Dune::dinfo << "multi domain grid function space:"
                << std::endl;

    Dune::dinfo << "( ";
    offset[0] = 0;
    maxlocalsize = 0;
    for (std::size_t i=0; i<BaseT::CHILDREN; i++)
      {
        Dune::dinfo << childGlobalSize[i] << " ";
        offset[i+1] = offset[i]+childGlobalSize[i];
        maxlocalsize += childLocalSize[i];
      }
    Dune::dinfo << ") total size = " << offset[BaseT::CHILDREN]
                << " max local size = " << maxlocalsize
                << std::endl;
  }

};


template<typename MultiGFS, typename ChildGFS>
struct get_subproblem_index
{
  static const std::size_t value = MultiGFS::template IndexForChild<ChildGFS>::value;
};


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


} // namespace MultiDomain

} // namespace PDELab

} // namespace Dune

#endif // DUNE_MULTIDOMAIN_MULTIDOMAINGRIDFUNCTIONSPACE_HH
