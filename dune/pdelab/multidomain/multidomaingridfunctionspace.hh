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

namespace {

  template<typename Entity>
  struct DataHandleSize
    : public Dune::PDELab::TypeTree::DirectChildrenVisitor
    , public Dune::PDELab::TypeTree::DynamicTraversal
  {

    template<typename LFS, typename Child, typename TreePath, typename ChildIndex>
    void beforeChild(const LFS& lfs, const Child& child, TreePath treePath, ChildIndex childIndex)
    {
      size += data_handle_size(child,typename gfs_flavor_tag<Child>::type());
    }

    template<typename Child>
    std::size_t data_handle_size(const Child& child, MultiDomainGFSTag tag)
    {
      return child.dataHandleSize(e);
    }

    template<typename Child>
    std::size_t data_handle_size(const Child& child, CouplingGFSTag tag)
    {
      return child.dataHandleSize(e);
    }

    template<typename Child>
    std::size_t data_handle_size(const Child& child, SubDomainGFSTag tag)
    {
      typedef typename Child::Traits::GridViewType::template Codim<Entity::codimension>::EntityPointer SDEP;
      typedef typename SDEP::Entity SDE;
      const SDEP ep = child.gridview().grid().subDomainEntityPointer(e);
      if (child.gridview().indexSet().contains(*ep))
        return child.dataHandleSize(*ep);
      else
        return 0;
    }

    DataHandleSize(const Entity& entity)
      : e(entity)
      , size(0)
    {}

    const Entity& e;
    std::size_t size;

  };


  template<typename Entity, typename Container, std::size_t gfs_depth>
  struct DataHandleGlobalIndices
    : public Dune::PDELab::DataHandleGlobalIndicesVisitor<Entity,Container,gfs_depth>
  {

    typedef Dune::PDELab::DataHandleGlobalIndicesVisitor<Entity,Container,gfs_depth> BaseT;
    using BaseT::pos;
    using BaseT::e;
    using BaseT::g;

    template<typename Child, typename TreePath>
    void leaf(const Child& child, TreePath treePath)
    {
      data_handle_global_indices(child,typename gfs_flavor_tag<Child>::type());
    }

    template<typename Child>
    void data_handle_global_indices(const Child& child, MultiDomainGFSTag tag)
    {
      pos += child.dataHandleGlobalIndices(e,g,pos,false);
    }

    template<typename Child>
    void data_handle_global_indices(const Child& child, CouplingGFSTag tag)
    {
      pos += child.dataHandleGlobalIndices(e,g,pos,false);
    }

    template<typename Child>
    void data_handle_global_indices(const Child& child, SubDomainGFSTag tag)
    {
      typedef typename Child::Traits::GridViewType::template Codim<Entity::codimension>::EntityPointer SDEP;
      typedef typename SDEP::Entity SDE;
      const SDEP ep = child.gridview().grid().subDomainEntityPointer(e);
      if (child.gridview().indexSet().contains(*ep))
        pos += child.dataHandleGlobalIndices(*ep,g,pos,false);
    }

    DataHandleGlobalIndices(const Entity& entity, Container& global)
      : BaseT(entity,global)
    {
    }

  };

}

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
    doVerify<T>(t.grid(),child.gridview().grid());
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
  PowerCompositeGridFunctionSpaceBase<MultiDomainGridFunctionSpace<G,Children...>,
                                      typename TypeTree::VariadicCompositeNode<Children...>::template Child<0>::Type::Traits::GridViewType,
                                      typename TypeTree::VariadicCompositeNode<Children...>::template Child<0>::Type::Traits::BackendType,
                                      GridFunctionSpaceLexicographicMapper,
                                      sizeof...(Children)
                                      >;

  typedef PowerCompositeGridFunctionSpaceBase<MultiDomainGridFunctionSpace<G,Children...>,
                                              typename TypeTree::VariadicCompositeNode<Children...>::template Child<0>::Type::Traits::GridViewType,
                                              typename TypeTree::VariadicCompositeNode<Children...>::template Child<0>::Type::Traits::BackendType,
                                              GridFunctionSpaceLexicographicMapper,
                                              sizeof...(Children)
                                              > ImplementationBase;

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
                                             void,
                                             !ChildEntryMap::has_duplicate_entries, // duplicate entries would lead to ambiguous type lookups
                                             sizeof...(Children)>
  Traits;

  typedef CompositeLexicographicOrdering<typename Traits::SizeType,typename Children::Ordering...> Ordering;

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


  MultiDomainGridFunctionSpace (G& g, Children&... children)
    : BaseT(children...)
    , _g(g)
    , _ordering(make_shared<Ordering>(*this,children.orderingPtr()...))
  {
    dune_static_assert(Dune::mdgrid::GridType<G>::v == Dune::mdgrid::multiDomainGrid,
                       "MultiDomainGridFunctionSpace only works on a MultiDomainGrid");
    TypeTree::applyToTree(*this,VerifyChildren());
  }

  template<typename EntityType>
  size_t dataHandleSize (const EntityType& e) const
  {
    DataHandleSize<EntityType> visitor(e);
    TypeTree::applyToTree(*this,visitor);
    return visitor.size;
  }

  //! return vector of global indices associated with the given entity
  template<typename EntityType, typename SizeType>
  void dataHandleGlobalIndices (const EntityType& e,
                                std::vector<SizeType>& global) const
  {
    global.resize(dataHandleSize(e));
    DataHandleGlobalIndices<EntityType,std::vector<SizeType>,TypeTree::TreeInfo<MultiDomainGridFunctionSpace>::depth> visitor(e,global);
    TypeTree::applyToTree(*this,visitor);
  }

  //------------------------------

  G& grid() {
    return _g;
  }

  const G& grid() const {
    return _g;
  }

  Ordering& ordering()
  {
    return *_ordering;
  }

  const Ordering& ordering() const
  {
    return *_ordering;
  }


private:

  const G& _g;
  shared_ptr<Ordering> _ordering;

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
