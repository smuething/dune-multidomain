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
    impl().compute_size(child,typename gfs_flavor_tag<Child>::type());
  }

  Impl& impl()
  {
    return static_cast<Impl&>(*this);
  }

  template<typename Child, typename Tag>
  void compute_size(Child& child, Tag tag)
  {}

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

  template<typename Child>
  void compute_size(Child& child, CouplingGFSTag tag)
  {
    if (child.gridFunctionSpace().contains(e))
      {
        // TODO: Fix this - it is definitely wrong!
        Dune::PDELab::ComputeSizeVisitor<Intersection> child_visitor(e,offset);
        Dune::PDELab::TypeTree::applyToTree(child,child_visitor);
        offset = child_visitor.offset;
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
    impl().fill_indices(child,typename gfs_flavor_tag<Child>::type());
    for (std::size_t i = 0; i<child.size(); ++i)
      (*lfs.global)[child.localIndex(i)] = lfs.pgfs->subMap(childIndex,(*lfs.global)[child.localIndex(i)]);
  }

  Impl& impl()
  {
    return static_cast<Impl&>(*this);
  }

  template<typename Child, typename Tag>
  void fill_indices(Child& child, Tag tag)
  {}

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

// local function space for a MultiDomainGridFunctionSpace
template<typename GFS,
         typename LFSTag,
         typename... Children>
class MultiDomainLocalFunctionSpaceNode
  : public LocalFunctionSpaceBaseNode<GFS>
  , public TypeTree::VariadicCompositeNode<Children...>
{

  typedef typename GFS::Traits::BackendType B;

  typedef LocalFunctionSpaceBaseNode<GFS> BaseT;

public:
  typedef MultiDomainLocalFunctionSpaceTraits<GFS,MultiDomainLocalFunctionSpaceNode> Traits;

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
  using BaseT::pgfs;

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

template<typename MultiDomainGFS, typename LFSTag>
struct MultiDomainLocalFunctionSpaceTransformationTemplate
{
  template<typename... TC>
  struct result
  {
    typedef MultiDomainLocalFunctionSpaceNode<MultiDomainGFS,LFSTag,TC...> type;
  };
};

template<typename MultiDomainGFS>
Dune::PDELab::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
  MultiDomainGFS,
  Dune::PDELab::gfs_to_lfs,
  MultiDomainLocalFunctionSpaceTransformationTemplate<MultiDomainGFS,StandardLFSTag>::template result
  >
lookupNodeTransformation(MultiDomainGFS*, Dune::PDELab::gfs_to_lfs*, MultiDomainGridFunctionSpaceTag);


struct gfs_to_coupling_lfs {};


template<typename MultiDomainGFS>
Dune::PDELab::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
  MultiDomainGFS,
  gfs_to_coupling_lfs,
  MultiDomainLocalFunctionSpaceTransformationTemplate<MultiDomainGFS,CouplingLFSTag>::template result
  >
lookupNodeTransformation(MultiDomainGFS*, gfs_to_coupling_lfs*, MultiDomainGridFunctionSpaceTag);


#if 0
// local function space description that can be bound to an element
// depends on a grid function space
template<typename GFS,typename... Children>
class MultiDomainLocalFunctionSpaceNode : public MultiDomainLocalFunctionSpaceNodeBase<GFS,StandardLFSTag,Children...>
{
  typedef MultiDomainLocalFunctionSpaceNodeBase<GFS,StandardLFSTag,Children...> BaseT;

public:
  typedef typename BaseT::Traits Traits;

  explicit MultiDomainLocalFunctionSpaceNode(shared_ptr<const GFS> gfs,
                                             shared_ptr<Children>... children)
    : BaseT(gfs,children...)
  {}

  explicit MultiDomainLocalFunctionSpaceNode(const GFS& gfs,
                                             shared_ptr<Children>... children)
    : BaseT(gfs,children)
    , TypeTree::VariadicCompositeNode<Children...>(children...)
  {}

  explicit MultiDomainLocalFunctionSpaceNode (const GFS& gfs)
    : BaseT(gfs), global_container(gfs.maxLocalSize())
  {}

  //! \brief bind local function space to entity


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

#endif

    //! \} group GridFunctionSpace

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_MULTIDOMAINLOCALFUNCTIONSPACE_HH
