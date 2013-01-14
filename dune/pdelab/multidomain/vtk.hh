#ifndef DUNE_PDELAB_MULTIDOMAIN_VTK_HH
#define DUNE_PDELAB_MULTIDOMAIN_VTK_HH

#include <dune/pdelab/gridfunctionspace/vtk.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


  template<typename SubDomain>
  struct subdomain_predicate
  {

    template<typename LFS>
    bool operator()(const LFS& lfs) const
    {
      return lfs.gridFunctionSpace().gridView().grid().domain() == _subDomain;
    }

    subdomain_predicate(const SubDomain& subDomain)
    : _subDomain(subDomain)
    {}

  private:
    SubDomain _subDomain;

  };


  //! Helper class for common data of a DGFTree.
  /**
   * This is a simple extension of the vanilla PDELab data
   * struct that also allows binding to a SubDomaingrid entity.
   */
  template<typename GFS, typename X, typename Pred>
  class DGFTreeCommonData
    : public Dune::PDELab::vtk::DGFTreeCommonData<GFS,X,Pred>
  {

    template<typename LFS, typename Data>
    friend class DGFTreeLeafFunction;

    template<typename LFS, typename Data>
    friend class DGFTreeVectorFunction;

    template<typename, typename>
    friend struct vtk_output_collector;

    typedef typename GFS::Traits::GridView::Grid MultiDomainGrid;
    typedef typename MultiDomainGrid::SubDomainGrid::template Codim<0>::Entity SubDomainCell;

  public:

    typedef GFS GridFunctionSpace;
    typedef X Vector;
    typedef Pred Predicate;

    DGFTreeCommonData(const GFS& gfs, const X& x)
      : Dune::PDELab::vtk::DGFTreeCommonData<GFS,X,Pred>(gfs,x)
    {}

    using Dune::PDELab::vtk::DGFTreeCommonData<GFS,X,Pred>::bind;

    void bind(const SubDomainCell& cell)
    {
      this->bind(MultiDomainGrid::multiDomainEntity(cell));
    }

  };

  template<typename VTKWriter,
           typename GFS,
           typename X,
           typename NameGenerator = vtk::DefaultFunctionNameGenerator,
           typename Predicate = vtk::DefaultPredicate>
  vtk::OutputCollector<
    VTKWriter,
    DGFTreeCommonData<GFS,X,Predicate>
    >
  addSolutionToVTKWriter(VTKWriter& vtk_writer,
                         const GFS& gfs,
                         const X& x,
                         const Predicate& predicate = Predicate(),
                         const NameGenerator& name_generator = vtk::defaultNameScheme())
  {
    typedef DGFTreeCommonData<GFS,X,Predicate> Data;
    vtk::OutputCollector<VTKWriter,Data> collector(vtk_writer,make_shared<Data>(gfs,x),predicate);
    collector.addSolution(name_generator);
    return std::move(collector);
  }

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_VTK_HH
