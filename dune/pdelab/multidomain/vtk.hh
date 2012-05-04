#ifndef DUNE_PDELAB_MULTIDOMAIN_VTK_HH
#define DUNE_PDELAB_MULTIDOMAIN_VTK_HH

#include <dune/pdelab/gridfunctionspace/vtk.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {

  //! Helper class for common data of a DGFTree.
  /**
   * This is a simple extension of the vanilla PDELab data
   * struct that also allows binding to a SubDomaingrid entity.
   */
  template<typename GFS, typename X, typename Pred>
  class DGFTreeCommonData
    : public Dune::PDELab::DGFTreeCommonData<GFS,X,Pred>
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
      : Dune::PDELab::DGFTreeCommonData<GFS,X,Pred>(gfs,x)
    {}

    using Dune::PDELab::DGFTreeCommonData<GFS,X,Pred>::bind;

    void bind(const SubDomainCell& cell)
    {
      this->bind(MultiDomainGrid::multiDomainEntity(cell));
    }

  };

  template<typename VTKWriter, typename GFS, typename X, typename Predicate = Dune::PDELab::default_predicate>
  vtk_output_collector<
    VTKWriter,
    DGFTreeCommonData<GFS,X,Predicate>
    >
  add_solution_to_vtk_writer(VTKWriter& vtk_writer, const GFS& gfs, const X& x, std::string base_name = "", const Predicate& predicate = Predicate())
  {
    typedef DGFTreeCommonData<GFS,X,Predicate> Data;
    vtk_output_collector<VTKWriter,Data> collector(vtk_writer,make_shared<Data>(gfs,x),predicate);
    collector.add_solution(base_name);
    return std::move(collector);
  }

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_VTK_HH
