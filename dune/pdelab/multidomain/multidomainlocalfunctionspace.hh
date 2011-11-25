// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_MULTIDOMAIN_MULTIDOMAINLOCALFUNCTIONSPACE_HH
#define DUNE_MULTIDOMAIN_MULTIDOMAINLOCALFUNCTIONSPACE_HH

#include <vector>
#include <dune/pdelab/common/typetree.hh>
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


struct MultiDomainGFSTag {};
struct SubDomainGFSTag {};
struct CouplingGFSTag {};
struct StandardLFSTag {};
struct CouplingLFSTag {};

template<typename GFS>
struct gfs_flavor_tag;

template<typename Entity, typename Impl>
struct ComputeSizeVisitorBase
  : public Dune::PDELab::TypeTree::DirectChildrenVisitor
  , public Dune::PDELab::TypeTree::DynamicTraversal
{

  template<typename Node, typename TreePath>
  void pre(Node& node, TreePath treePath)
  {
    node.offset = offset;
  }

  template<typename Node, typename TreePath>
  void post(Node& node, TreePath treePath)
  {
    node.n = offset - node.offset;
  }

  template<typename LFS, typename Child, typename TreePath, typename ChildIndex>
  void beforeChild(LFS& lfs, Child& child, TreePath treePath, ChildIndex childIndex)
  {
    impl().compute_size(child,typename gfs_flavor_tag<typename Child::Traits::GridFunctionSpaceType>::type());
  }

  Impl& impl()
  {
    return static_cast<Impl&>(*this);
  }

  ComputeSizeVisitorBase(const Entity& entity, std::size_t offset_)
    : e(entity)
    , offset(offset_)
  {}

  const Entity& e;
  std::size_t offset;

};


template<typename Entity, typename LFSTag>
struct ComputeSizeVisitor;

template<typename Entity>
struct ComputeSizeVisitor<Entity,StandardLFSTag>
  : public ComputeSizeVisitorBase<Entity, ComputeSizeVisitor<Entity,StandardLFSTag> >
{

  typedef ComputeSizeVisitorBase<Entity, ComputeSizeVisitor<Entity,StandardLFSTag> > BaseT;
  using BaseT::e;
  using BaseT::offset;

  ComputeSizeVisitor(const Entity& entity, std::size_t offset_ = 0)
    : BaseT(entity,offset_)
  {}

  template<typename Child, typename Tag>
  void compute_size(Child& child, Tag tag)
  {}

  template<typename Child>
  void compute_size(Child& child, MultiDomainGFSTag tag)
  {
    Dune::PDELab::ComputeSizeVisitor<Entity> child_visitor(e,offset);
    Dune::PDELab::TypeTree::applyToTree(child,child_visitor);
    offset = child_visitor.offset;
  }

  template<typename Child>
  void compute_size(Child& child, SubDomainGFSTag tag)
  {
    typedef typename Child::Traits::GridViewType::template Codim<0>::EntityPointer SDEP;
    typedef typename SDEP::Entity SDE;
    const SDEP ep = child.gridFunctionSpace().gridview().grid().subDomainEntityPointer(e);
    if (child.gridFunctionSpace().gridview().indexSet().contains(*ep))
      {
        Dune::PDELab::ComputeSizeVisitor<SDE> child_visitor(*ep,offset);
        Dune::PDELab::TypeTree::applyToTree(child,child_visitor);
        offset = child_visitor.offset;
      }
    else
      {
        // make sure all children know they are empty
        Dune::PDELab::ClearSizeVisitor<> child_visitor(offset);
        Dune::PDELab::TypeTree::applyToTree(child,child_visitor);
      }
  }

};


template<typename Intersection>
struct ComputeSizeVisitor<Intersection,CouplingLFSTag>
  : public ComputeSizeVisitorBase<Intersection, ComputeSizeVisitor<Intersection,CouplingLFSTag> >
{

  typedef ComputeSizeVisitorBase<Intersection, ComputeSizeVisitor<Intersection,CouplingLFSTag> > BaseT;
  using BaseT::e;
  using BaseT::offset;

  ComputeSizeVisitor(const Intersection& intersection, std::size_t offset_ = 0)
    : BaseT(intersection,offset_)
  {}

  template<typename Child, typename Tag>
  void compute_size(Child& child, Tag tag)
  {}

  template<typename Child>
  void compute_size(Child& child, CouplingGFSTag tag)
  {
    if (child.gridFunctionSpace().contains(e))
      {
        Dune::PDELab::ComputeSizeVisitor<Intersection> child_visitor(e,offset);
        Dune::PDELab::TypeTree::applyToTree(child,child_visitor);
        offset = child_visitor.offset;
      }
    else
      {
        // make sure all children know they are empty
        Dune::PDELab::ClearSizeVisitor<> child_visitor(offset);
        Dune::PDELab::TypeTree::applyToTree(child,child_visitor);
      }
  }

};


template<typename Impl, typename GFS, typename Container>
struct FillIndicesVisitorBase
  : public Dune::PDELab::TypeTree::DirectChildrenVisitor
  , public Dune::PDELab::TypeTree::DynamicTraversal
{

  template<typename LFS, typename Child, typename TreePath, typename ChildIndex>
  void afterChild(LFS& lfs, Child& child, TreePath treePath, ChildIndex childIndex)
  {
    impl().fill_indices(child,typename gfs_flavor_tag<typename Child::Traits::GridFunctionSpaceType>::type());
    for (std::size_t i = 0; i<child.size(); ++i)
      {
        (*lfs.global)[child.localIndex(i)] = lfs.pgfs->subMap(childIndex,(*lfs.global)[child.localIndex(i)]);
        (*lfs._multi_indices)[child.localIndex(i)].push_back(childIndex);
      }
  }

  Impl& impl()
  {
    return static_cast<Impl&>(*this);
  }

  FillIndicesVisitorBase(const GFS& gfs_, Container& container_)
    : gfs(gfs_)
    , container(container_)
  {}

  const GFS& gfs;
  Container& container;

};


template<typename Entity, typename GFS, typename Container, typename LFSTag>
struct FillIndicesVisitor;

template<typename Entity, typename GFS, typename Container>
struct FillIndicesVisitor<Entity,GFS,Container,StandardLFSTag>
  : public FillIndicesVisitorBase<FillIndicesVisitor<
                                    Entity,
                                    GFS,
                                    Container,
                                    StandardLFSTag>,
                                  GFS,
                                  Container>
{

  typedef FillIndicesVisitorBase<FillIndicesVisitor<
                                   Entity,
                                   GFS,
                                   Container,
                                   StandardLFSTag>,
                                 GFS,
                                 Container> BaseT;

  FillIndicesVisitor(const Entity& entity, const GFS& gfs, Container& container)
    : BaseT(gfs,container)
    , e(entity)
  {}

  const Entity& e;

  template<typename Child, typename Tag>
  void fill_indices(Child& child, Tag tag)
  {}

  template<typename Child>
  void fill_indices(Child& child, MultiDomainGFSTag tag)
  {
    Dune::PDELab::FillIndicesVisitor<Entity> child_visitor(e);
    Dune::PDELab::TypeTree::applyToTree(child,child_visitor);
  }

  template<typename Child>
  void fill_indices(Child& child, SubDomainGFSTag tag)
  {
    typedef typename Child::Traits::GridViewType::template Codim<0>::EntityPointer SDEP;
    typedef typename SDEP::Entity SDE;
    const SDEP ep = child.gridFunctionSpace().gridview().grid().subDomainEntityPointer(e);
    if (child.gridFunctionSpace().gridview().indexSet().contains(*ep))
      {
        Dune::PDELab::FillIndicesVisitor<SDE> child_visitor(*ep);
        Dune::PDELab::TypeTree::applyToTree(child,child_visitor);
      }
  }

};


template<typename Intersection, typename GFS, typename Container>
struct FillIndicesVisitor<Intersection,GFS,Container,CouplingLFSTag>
  : public FillIndicesVisitorBase<FillIndicesVisitor<
                                    Intersection,
                                    GFS,
                                    Container,
                                    CouplingLFSTag>,
                                  GFS,
                                  Container>
{

  typedef FillIndicesVisitorBase<FillIndicesVisitor<
                                   Intersection,
                                   GFS,
                                   Container,
                                   CouplingLFSTag>,
                                 GFS,
                                 Container> BaseT;

  FillIndicesVisitor(const Intersection& is_, const GFS& gfs, Container& container)
    : BaseT(gfs,container)
    , is(is_)
  {}

  const Intersection& is;

  template<typename Child, typename Tag>
  void fill_indices(Child& child, Tag tag)
  {}

  template<typename Child>
  void fill_indices(Child& child, CouplingGFSTag tag)
  {
    if (child.gridFunctionSpace().contains(is))
      {
        // TODO: Fix this - it is definitely wrong!
        Dune::PDELab::FillIndicesVisitor<Intersection> child_visitor(is);
        Dune::PDELab::TypeTree::applyToTree(child,child_visitor);
      }
  }

};





template<typename GFS, typename MI, typename N>
struct MultiDomainLocalFunctionSpaceTraits
  : public LocalFunctionSpaceBaseTraits<GFS,MI>
{

  //! type of local function space node
  typedef N NodeType;

  typedef typename GFS::Traits::GridType GridType;

  typedef typename GFS::Traits::GridType Grid;

  //! \brief Type of intersection in the grid
  typedef typename GFS::Traits::GridViewType::Intersection Intersection;

};

// local function space for a MultiDomainGridFunctionSpace
template<typename GFS,
         typename LFSTag,
         typename MultiIndex,
         typename... Children>
class MultiDomainLocalFunctionSpaceNode
  : public LocalFunctionSpaceBaseNode<GFS,MultiIndex>
  , public TypeTree::VariadicCompositeNode<Children...>
{

  typedef typename GFS::Traits::BackendType B;

  typedef LocalFunctionSpaceBaseNode<GFS,MultiIndex> BaseT;

public:
  typedef MultiDomainLocalFunctionSpaceTraits<GFS,MultiIndex,MultiDomainLocalFunctionSpaceNode> Traits;

public:

  template<typename Transformation>
  MultiDomainLocalFunctionSpaceNode(shared_ptr<const GFS> gfs,
                                    const Transformation& t,
                                    shared_ptr<Children>... children)
    : BaseT(gfs)
    , TypeTree::VariadicCompositeNode<Children...>(children...)
  {}

  template<typename Transformation>
  MultiDomainLocalFunctionSpaceNode(const GFS& gfs,
                                    const Transformation& t,
                                    shared_ptr<Children>... children)
    : BaseT(stackobject_to_shared_ptr(gfs))
    , TypeTree::VariadicCompositeNode<Children...>(children...)
  {}

  using BaseT::global_storage;
  using BaseT::n;
  using BaseT::offset;
  using BaseT::global;
  using BaseT::_multi_indices;
  using BaseT::pgfs;

  shared_ptr<const GFS> gridFunctionSpaceStorage() const
  {
    return pgfs;
  }

  template<typename Element>
  void bind (const Element& e)
  {
    ComputeSizeVisitor<Element,LFSTag> csv(e);
    Dune::PDELab::TypeTree::applyToTree(*this,csv);

    global_storage.resize(n);

    // initialize iterators and fill indices
    FillIndicesVisitor<Element,GFS,typename BaseT::Traits::IndexContainer,LFSTag> fiv(e,*pgfs,global_storage);
    Dune::PDELab::TypeTree::applyToTree(*this,fiv);

    // apply upMap
    for (typename Traits::IndexContainer::size_type i=0; i<offset; ++i)
      global_storage[i] = pgfs->upMap(global_storage[i]);
  }

};

struct MultiDomainGridFunctionSpaceTag {};

template<typename MultiDomainGFS, typename LFSTag, typename Transformation>
struct MultiDomainLocalFunctionSpaceTransformationTemplate
{
  template<typename... TC>
  struct result
  {
    typedef MultiDomainLocalFunctionSpaceNode<MultiDomainGFS,
                                              LFSTag,
                                              typename Transformation::MultiIndex,
                                              TC...> type;
  };
};

template<typename MultiDomainGFS, typename data>
Dune::PDELab::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
  MultiDomainGFS,
  Dune::PDELab::gfs_to_lfs<data>,
  MultiDomainLocalFunctionSpaceTransformationTemplate<
    MultiDomainGFS,
    StandardLFSTag,
    gfs_to_lfs<data>
    >::template result
  >
lookupNodeTransformation(MultiDomainGFS*, Dune::PDELab::gfs_to_lfs<data>*, MultiDomainGridFunctionSpaceTag);


template<typename GFS>
struct gfs_to_coupling_lfs {

  //! The MultiIndex type that will be used in the resulting LocalFunctionSpace tree.
  typedef Dune::PDELab::MultiIndex<std::size_t,TypeTree::TreeInfo<GFS>::depth> MultiIndex;

};


template<typename MultiDomainGFS,typename data>
Dune::PDELab::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
  MultiDomainGFS,
  gfs_to_coupling_lfs<data>,
  MultiDomainLocalFunctionSpaceTransformationTemplate<
    MultiDomainGFS,
    CouplingLFSTag,
    gfs_to_coupling_lfs<data>
    >::template result
  >
    lookupNodeTransformation(MultiDomainGFS*, gfs_to_coupling_lfs<data>*, MultiDomainGridFunctionSpaceTag);



template<typename PowerGridFunctionSpace, typename data>
Dune::PDELab::TypeTree::TemplatizedGenericPowerNodeTransformation<
  PowerGridFunctionSpace,
  gfs_to_coupling_lfs<data>,
  power_gfs_to_lfs_template<
    PowerGridFunctionSpace,
    gfs_to_coupling_lfs<data>
    >::template result
  >
lookupNodeTransformation(PowerGridFunctionSpace* pgfs, gfs_to_coupling_lfs<data>* t, PowerGridFunctionSpaceTag tag);

#if HAVE_VARIADIC_TEMPLATES

template<typename CompositeGridFunctionSpace, typename data>
Dune::PDELab::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
  CompositeGridFunctionSpace,
  gfs_to_coupling_lfs<data>,
  variadic_composite_gfs_to_lfs_template<
    CompositeGridFunctionSpace,
    gfs_to_coupling_lfs<data>
    >::template result
  >
lookupNodeTransformation(CompositeGridFunctionSpace* cgfs, gfs_to_coupling_lfs<data>* t, CompositeGridFunctionSpaceTag tag);

#else
#error "require variadic templates"
#endif

template<typename GridFunctionSpace, typename data>
Dune::PDELab::TypeTree::GenericLeafNodeTransformation<
  GridFunctionSpace,
  gfs_to_coupling_lfs<data>,
  LeafLocalFunctionSpaceNode<
    GridFunctionSpace,
    typename gfs_to_coupling_lfs<data>::MultiIndex
    >
  >
lookupNodeTransformation(GridFunctionSpace* gfs, gfs_to_coupling_lfs<data>* t, LeafGridFunctionSpaceTag tag);




// TODO: Add support for tags
template<typename GFS, typename Tag=AnySpaceTag>
class CouplingLocalFunctionSpace
  : public Dune::PDELab::TypeTree::TransformTree<GFS,gfs_to_coupling_lfs<GFS> >::Type
{
  typedef typename Dune::PDELab::TypeTree::TransformTree<GFS,gfs_to_coupling_lfs<GFS> >::Type BaseT;

public:
  CouplingLocalFunctionSpace(const GFS & gfs) : BaseT(TypeTree::TransformTree<GFS,gfs_to_coupling_lfs<GFS> >::transform(gfs)) { this->setup(*this); }
};

    //! \} group GridFunctionSpace

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_MULTIDOMAINLOCALFUNCTIONSPACE_HH
