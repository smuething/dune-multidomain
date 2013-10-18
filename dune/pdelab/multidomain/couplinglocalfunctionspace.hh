#ifndef DUNE_MULTIDOMAIN_COUPLINGLOCALFUNCTIONSPACE_HH
#define DUNE_MULTIDOMAIN_COUPLINGLOCALFUNCTIONSPACE_HH

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include <dune/pdelab/multidomain/dofmapper.hh>
#include <dune/pdelab/multidomain/multidomainlocalfunctionspace.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {


template <typename GFS, typename DI>
class CouplingLocalFunctionSpaceNode;


template<typename GFS, typename DI>
struct CouplingLocalFunctionSpaceTraits
  : public LocalFunctionSpaceBaseTraits<GFS,DI>
{

  typedef typename GFS::Traits::GridViewType::Intersection Intersection;

};


template<typename GFS, typename DI>
struct LeafCouplingLocalFunctionSpaceTraits
  : public CouplingLocalFunctionSpaceTraits<GFS,DI>
{

  //! \brief Type of local finite element
  typedef typename GFS::Traits::FiniteElementType FiniteElementType;

  //! \brief Type of constraints engine
  typedef typename GFS::Traits::ConstraintsType ConstraintsType;

  typedef CouplingLocalFunctionSpaceNode<GFS,DI> NodeType;

};



struct PowerCouplingLocalFunctionSpaceTag {};
struct PowerCouplingGridFunctionSpaceTag {};

// local function space for a power grid function space
template<typename GFS, typename DI, typename ChildLFS, std::size_t k>
class PowerCouplingLocalFunctionSpaceNode
  : public LocalFunctionSpaceBaseNode<GFS,DI>
  , public TypeTree::PowerNode<ChildLFS,k>
{
  typedef LocalFunctionSpaceBaseNode<GFS,DI> BaseT;
  typedef TypeTree::PowerNode<ChildLFS,k> TreeNode;

public:
  typedef CouplingLocalFunctionSpaceTraits<GFS,DI> Traits;

  typedef PowerCouplingLocalFunctionSpaceTag ImplementationTag;

  //! \brief initialize with grid function space
  template<typename Transformation>
  PowerCouplingLocalFunctionSpaceNode (shared_ptr<const GFS> gfs,
                                       const Transformation& t,
                                       const array<shared_ptr<ChildLFS>,k>& children)
    : BaseT(gfs)
    , TreeNode(children)
  {}

  template<typename Transformation>
  PowerCouplingLocalFunctionSpaceNode (const GFS& gfs,
                                       const Transformation& t,
                                       const array<shared_ptr<ChildLFS>,k>& children)
    : BaseT(stackobject_to_shared_ptr(gfs))
    , TreeNode(children)
  {}

};


template<typename SourceNode, typename Transformation>
struct power_coupling_gfs_to_coupling_lfs_template
{
  template<typename TC>
  struct result
  {
    typedef PowerCouplingLocalFunctionSpaceNode<SourceNode,typename Transformation::DOFIndex,TC,SourceNode::CHILDREN> type;
  };
};

template<typename PowerCouplingGridFunctionSpace, typename data>
TypeTree::TemplatizedGenericPowerNodeTransformation<
  PowerCouplingGridFunctionSpace,
  gfs_to_coupling_lfs<data>,
  power_coupling_gfs_to_coupling_lfs_template<
    PowerCouplingGridFunctionSpace,
    gfs_to_coupling_lfs<data>
    >::template result
  >
registerNodeTransformation(PowerCouplingGridFunctionSpace* gfs, gfs_to_coupling_lfs<data>* t, PowerCouplingGridFunctionSpaceTag* tag);


template<typename SourceNode, typename Transformation>
struct power_coupling_gfs_to_lfs_template
{
  template<typename TC>
  struct result
  {
    typedef PowerCouplingLocalFunctionSpaceNode<SourceNode,typename Transformation::DOFIndex,TC,SourceNode::CHILDREN> type;
  };
};

template<typename PowerCouplingGridFunctionSpace, typename data>
TypeTree::TemplatizedGenericPowerNodeTransformation<
  PowerCouplingGridFunctionSpace,
  gfs_to_lfs<data>,
  power_coupling_gfs_to_lfs_template<
    PowerCouplingGridFunctionSpace,
    gfs_to_lfs<data>
    >::template result
  >
registerNodeTransformation(PowerCouplingGridFunctionSpace* gfs, gfs_to_lfs<data>* t, PowerCouplingGridFunctionSpaceTag* tag);


struct CouplingLocalFunctionSpaceTag {};
struct CouplingGridFunctionSpaceTag {};

template <typename GFS, typename DI>
class CouplingLocalFunctionSpaceNode
  : public LocalFunctionSpaceBaseNode<GFS,DI>
  , public TypeTree::LeafNode
{
  typedef typename GFS::Traits::BackendType B;

  typedef LocalFunctionSpaceBaseNode<GFS,DI> BaseT;

  using BaseT::n;

public:
  typedef LeafCouplingLocalFunctionSpaceTraits<GFS,DI> Traits;
  using BaseT::dofIndex;

private:

  typedef FiniteElementInterfaceSwitch<
    typename Traits::FiniteElementType
    > FESwitch;

public:

  typedef CouplingLocalFunctionSpaceTag ImplementationTag;

  template<typename Transformation>
  CouplingLocalFunctionSpaceNode (shared_ptr<const GFS> gfs, const Transformation& t)
    : BaseT(gfs)
  {}

  template<typename Transformation>
  CouplingLocalFunctionSpaceNode (const GFS& gfs, const Transformation& t)
    : BaseT(stackobject_to_shared_ptr(gfs))
  {}


  /*
  //! Calculates the multiindices associated with the given entity.
  template<typename Intersection, typename MultiIndexIterator>
  void multiIndices(const Intersection& is, MultiIndexIterator it, MultiIndexIterator endit)
  {
    if (!this->pgfs->contains(is))
      {
        assert(it == endit && "Intersection not contained in CouplingGFS, but local size > 0");
        return;
      }

    // get layout of entity
    const typename FESwitch::Coefficients &coeffs =
      FESwitch::coefficients(*pfe);

    typedef typename GFS::Traits::GridViewType GV;
    GV gv = this->gridFunctionSpace().gridView();

    DOFMapper<GV> dm(is);

    const Dune::ReferenceElement<typename GV::ctype,GV::Grid::dimension-1>& refEl =
      Dune::ReferenceElements<typename GV::ctype,GV::Grid::dimension-1>::general(this->pfe->type());

    for (std::size_t i = 0; i < std::size_t(coeffs.size()); ++i, ++it)
      {
        // get geometry type of subentity
        Dune::GeometryType gt = refEl.type(coeffs.localKey(i).subEntity(),
                                           coeffs.localKey(i).codim());

        // evaluate consecutive index of subentity
        typename GV::IndexSet::IndexType index = gv.indexSet().subIndex(dm.element(),
                                                                        dm.mapSubIndex(coeffs.localKey(i).subEntity(),coeffs.localKey(i).codim()),
                                                                        coeffs.localKey(i).codim() + 1);

        GFS::Ordering::Traits::DOFIndexAccessor::store(*it,gt,index,coeffs.localKey(i).index());

        // make sure we don't write past the end of the iterator range
        assert(it != endit);
      }
  }
  */

  //! Calculates the multiindices associated with the given entity.
  template<typename Intersection, typename DOFIndexIterator>
  void dofIndices(const Intersection& is, DOFIndexIterator it, DOFIndexIterator endit)
  {

    if (!this->pgfs->contains(is))
      {
        assert(it == endit && "Intersection not contained in CouplingGFS, but local size > 0");
        return;
      }

    // get layout of entity
    const typename FESwitch::Coefficients &coeffs =
      FESwitch::coefficients(*pfe);

    typedef typename GFS::Traits::GridViewType GV;
    GV gv = this->gridFunctionSpace().gridView();

    DOFMapper<GV> dm(is);

    const Dune::ReferenceElement<typename GV::ctype,GV::Grid::dimension-1>& refEl =
      Dune::ReferenceElements<typename GV::ctype,GV::Grid::dimension-1>::general(this->pfe->type());

    for (std::size_t i = 0; i < std::size_t(coeffs.size()); ++i, ++it)
      {
        // get geometry type of subentity
        Dune::GeometryType gt = refEl.type(coeffs.localKey(i).subEntity(),
                                           coeffs.localKey(i).codim());

        // evaluate consecutive index of subentity
        typename GV::IndexSet::IndexType index = gv.indexSet().subIndex(dm.element(),
                                                                        dm.mapSubIndex(coeffs.localKey(i).subEntity(),coeffs.localKey(i).codim()),
                                                                        coeffs.localKey(i).codim() + 1);

        // store data
        GFS::Ordering::Traits::DOFIndexAccessor::store(*it,gt,index,coeffs.localKey(i).index());

        // make sure we don't write past the end of the iterator range
        assert(it != endit);
      }
  }

  //! \brief bind local function space to entity
  /**

     This is a generic implementation of the bind function. It is
     parametrized with the NodeType, which the type of the derived
     LocalFunctionSpaceNode. Handing the NodeType as a parammeter
     avoid the need for the CRTP construct, but all derived
     classes have to add a method bind, which forward to this
     method.

     \param node reference to the derived node, the address must be the same as this
     \param e entity to bind to
  */
  template<typename NodeType>
  bool bind (NodeType& node, const typename Traits::Intersection& is)
  {
    // we should only call bind on out selfs
    assert(&node == this);

    typedef typename Traits::Intersection Intersection;

    if (!this->gridFunctionSpace().contains(is))
      return false;

    Dune::PDELab::ComputeSizeVisitor<Intersection> csv(is);
    TypeTree::applyToTree(*this,csv);

    this->_dof_index_storage.resize(node.n);

    Dune::PDELab::FillIndicesVisitor<Intersection> fiv(is);
    TypeTree::applyToTree(*this,fiv);

    // apply upMap
    //for (typename Traits::IndexContainer::size_type i=0; i<n; ++i)
    //  global_storage[i] = this->gridFunctionSpace().upMap(global_storage[i]);

    return true;
  }


  //! \brief get local finite element
  const typename Traits::FiniteElementType& finiteElement () const
  {
    return *pfe;
  }

  //! \brief get constraints engine
  const typename Traits::ConstraintsType& constraints () const
  {
    return this->pgfs->constraints();
  }

  /** \brief write back coefficients for one element to container */
  template<typename GC, typename LC>
  void mwrite (const LC& lc, GC& gc) const
  {
    // LC and GC are maps of maps
    typedef typename LC::const_iterator local_col_iterator;
    typedef typename LC::value_type::second_type local_row_type;
    typedef typename local_row_type::const_iterator local_row_iterator;
    typedef typename GC::iterator global_col_iterator;
    typedef typename GC::value_type::second_type global_row_type;

    for (local_col_iterator cit=lc.begin(); cit!=lc.end(); ++cit)
      {
        typename Traits::SizeType i = globalIndex(cit->first);
        // insert empty row in global container if necessary
        global_col_iterator gcit = gc.find(i);
        if (gcit==gc.end())
          gc[i] = global_row_type();

        // copy row to global container with transformed indices
        for (local_row_iterator rit=(cit->second).begin(); rit!=(cit->second).end(); ++rit)
          gc[i][globalIndex(rit->first)] = rit->second;
      }
  }

private:

  typename FESwitch::Store pfe;
};

template<typename GridFunctionSpace, typename data>
TypeTree::GenericLeafNodeTransformation<
  GridFunctionSpace,
  gfs_to_coupling_lfs<data>,
  CouplingLocalFunctionSpaceNode<
    GridFunctionSpace,
    typename gfs_to_coupling_lfs<data>::DOFIndex
    >
  >
registerNodeTransformation(GridFunctionSpace* gfs, gfs_to_coupling_lfs<data>* t, CouplingGridFunctionSpaceTag* tag);

template<typename GridFunctionSpace, typename data>
TypeTree::GenericLeafNodeTransformation<
  GridFunctionSpace,
  gfs_to_lfs<data>,
  CouplingLocalFunctionSpaceNode<
    GridFunctionSpace,
    typename gfs_to_lfs<data>::DOFIndex
    >
  >
registerNodeTransformation(GridFunctionSpace* gfs, gfs_to_lfs<data>* t, CouplingGridFunctionSpaceTag* tag);

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_COUPLINGLOCALFUNCTIONSPACE_HH
