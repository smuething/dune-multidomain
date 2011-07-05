#ifndef DUNE_MULTIDOMAIN_COUPLINGLOCALFUNCTIONSPACE_HH
#define DUNE_MULTIDOMAIN_COUPLINGLOCALFUNCTIONSPACE_HH

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {


template <typename GFS>
class CouplingLocalFunctionSpaceNode;

template<typename GFS>
struct CouplingLocalFunctionSpaceTraits
{
  //! \brief the grid view where grid function is defined upon
  typedef GFS GridFunctionSpaceType;

  //! \brief Type to store indices from Backend
  typedef typename GFS::Traits::GridViewType GridViewType;

  //! \brief Type of codim 0 entity in the grid
  typedef typename GridViewType::Traits::template Codim<0>::Entity Element;

  typedef typename GridViewType::Traits::Intersection Intersection;

  //! \brief Type to store indices from Backend
  typedef typename GFS::Traits::SizeType SizeType;

  //! \brief Type of container to store indices
  typedef typename std::vector<SizeType> IndexContainer;

  //! \brief Type of local finite element
  typedef typename GFS::Traits::FiniteElementType FiniteElementType;

  //! \brief Type of constraints engine
  typedef typename GFS::Traits::ConstraintsType ConstraintsType;

  typedef CouplingLocalFunctionSpaceNode<GFS> NodeType;

};



struct CouplingLocalFunctionSpaceTag {};
struct CouplingGridFunctionSpaceTag {};


template <typename GFS>
class CouplingLocalFunctionSpaceNode
  : public LocalFunctionSpaceBaseNode<GFS>
  , public TypeTree::LeafNode
{
  typedef typename GFS::Traits::BackendType B;

  typedef LocalFunctionSpaceBaseNode<GFS> BaseT;

  template<typename>
  friend struct Dune::PDELab::PropagateGlobalStorageVisitor;

  template<typename>
  friend struct Dune::PDELab::ClearSizeVisitor;

  template<typename>
  friend struct Dune::PDELab::ComputeSizeVisitor;

  template<typename>
  friend struct Dune::PDELab::FillIndicesVisitor;

  using BaseT::n;
  using BaseT::global_storage;
  using BaseT::globalIndex;

public:
  typedef CouplingLocalFunctionSpaceTraits<GFS> Traits;

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

    if (!pgfs->contains(is))
      return false;

    ComputeSizeVisitor<Intersection> csv(is);
    TypeTree::applyToTree(*this,csv);

    global_storage.resize(n);

    FillIndicesVisitor<Intersection> fiv(is);
    TypeTree::applyToTree(*this,fiv);

    // apply upMap
    for (typename Traits::IndexContainer::size_type i=0; i<n; ++i)
      global_storage[i] = pgfs->upMap(global_storage[i]);

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

  typename FESwitch::Store plfem;
};

template<typename GridFunctionSpace>
Dune::PDELab::TypeTree::GenericLeafNodeTransformation<GridFunctionSpace,gfs_to_lfs,CouplingLocalFunctionSpaceNode<GridFunctionSpace> >
lookupNodeTransformation(GridFunctionSpace* gfs, gfs_to_lfs* t, CouplingGridFunctionSpaceTag tag);

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_COUPLINGLOCALFUNCTIONSPACE_HH
