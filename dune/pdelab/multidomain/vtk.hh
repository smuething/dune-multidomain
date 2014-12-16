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


  namespace detail {

    struct ScalarDiscreteCalculationProviderFlagsBase
    {};

    template<bool enable = false>
    struct ScalarDiscreteCalculationProviderFlags
      : public std::conditional<enable,
                                ScalarDiscreteCalculationProviderFlagsBase,
                                ScalarDiscreteCalculationProviderFlags<true>
                                >::type
    {

      typedef ScalarDiscreteCalculationProviderFlags<true> Computations;

      static constexpr bool evaluateBasis()
      {
        return enable;
      }

      static constexpr bool evaluateSolution()
      {
        return enable;
      }

      static constexpr bool evaluateBasisJacobians()
      {
        return enable;
      }

      static constexpr bool evaluateBasisGradients()
      {
        return enable;
      }

      static constexpr bool evaluateGradient()
      {
        return enable;
      }

    };

  }


  template<
    typename LFS,
    typename RangeField = typename Dune::BasisInterfaceSwitch<
      typename Dune::FiniteElementInterfaceSwitch<
        typename LFS::Traits::FiniteElement
        >::Basis
      >::RangeField
    >
  struct ScalarFunctionPolicy
  {

    typedef Dune::PDELab::GridFunctionTraits<
      typename LFS::Traits::GridView,
      RangeField,
      Dune::BasisInterfaceSwitch<
        typename Dune::FiniteElementInterfaceSwitch<
          typename LFS::Traits::FiniteElement
          >::Basis
        >::dimRange,
      typename Dune::BasisInterfaceSwitch<
        typename Dune::FiniteElementInterfaceSwitch<
          typename LFS::Traits::FiniteElement
          >::Basis
        >::Range
      > Traits;

  };


  template<
    typename LFS,
    typename RangeField = typename Dune::BasisInterfaceSwitch<
      typename Dune::FiniteElementInterfaceSwitch<
        typename LFS::Traits::FiniteElement
        >::Basis
      >::RangeField
    >
  struct CoordinateFunctionPolicy
  {

    typedef Dune::PDELab::GridFunctionTraits<
      typename LFS::Traits::GridView,
      RangeField,
      Dune::FiniteElementInterfaceSwitch<
        typename LFS::Traits::FiniteElement
        >::Basis::Traits::dimDomain,
      Dune::FieldVector<
        RangeField,
        Dune::FiniteElementInterfaceSwitch<
          typename LFS::Traits::FiniteElement
          >::Basis::Traits::dimDomain
        >
      > Traits;

  };


  template<typename LFS, typename Data, typename Impl>
  class ScalarDiscreteCalculationProvider
    : public detail::ScalarDiscreteCalculationProviderFlags<>
  {

  public:

    struct Traits
    {

      typedef typename LFS::Traits::GridView GridView;

      typedef typename GridView::template Codim<0>::Entity Cell;

      typedef typename LFS::Traits::FiniteElement FiniteElement;

      typedef Dune::FiniteElementInterfaceSwitch<
        FiniteElement
        > FESwitch;

      typedef typename FESwitch::Basis Basis;

      typedef Dune::BasisInterfaceSwitch<
        Basis
        > BasisSwitch;

      typedef typename Basis::Traits::DomainFieldType DomainField;
      typedef typename Basis::Traits::DomainType Domain;
      static const std::size_t dimDomain = Basis::Traits::dimDomain;

      typedef typename Basis::Traits::RangeFieldType RangeField;
      typedef typename Basis::Traits::RangeType Range;
      static const std::size_t dimRange = Basis::Traits::dimRange;

      typedef typename Basis::Traits::JacobianType Jacobian;

      typedef Dune::FieldMatrix<
        RangeField,
        GridView::dimensionworld,
        dimDomain
        > GlobalJacobian;

    };

    ScalarDiscreteCalculationProvider(const LFS& lfs, const std::shared_ptr<Data>& data)
      : _lfs(lfs)
      , _data(data)
      , _basis(lfs.maxSize())
      , _basis_jacobians(lfs.maxSize())
      , _basis_global_jacobians(lfs.maxSize())
    {}

    void bind(const typename Traits::Cell& e,
              const typename Traits::Domain& x) const
    {

      _data->bind(e);

      if (impl().evaluateBasis() || impl().evaluateSolution())
        Traits::FESwitch::basis(_lfs.finiteElement()).evaluateFunction(x,_basis);

      if (impl().evaluateSolution())
        {
          _solution = typename Traits::Range(0);

          for (std::size_t i = 0; i < _lfs.size(); ++i)
            _solution.axpy(_data->_x_local(_lfs,i),_basis[i]);
        }

      if (impl().evaluateBasisJacobians() ||
          impl().evaluateBasisGradients() ||
          impl().evaluateGradient())
        {
          Traits::FESwitch::basis(_lfs.finiteElement()).evaluateJacobian(x,_basis_jacobians);
        }

      if (impl().evaluateBasisGradients() ||
          impl().evaluateGradient())
        {
          auto geometry = e.geometry();
          auto JgeoIT = geometry.jacobianInverseTransposed(x);

          for (std::size_t i = 0; i < _lfs.size(); ++i)
            JgeoIT.mv(_basis_jacobians[i][0],_basis_global_jacobians[i][0]);

        }

      if (impl().evaluateGradient())
        {
          _gradient = typename Traits::GlobalJacobian(0);

          for (std::size_t i = 0; i < _lfs.size(); ++i)
            _gradient[0].axpy(_data->_x_local(_lfs,i),_basis_global_jacobians[i][0]);

        }

      if (impl().evaluateGradient())
        {
          _gradient = typename Traits::GlobalJacobian(0);

          for (std::size_t i = 0; i < _lfs.size(); ++i)
            _gradient[0].axpy(_data->_x_local(_lfs,i),_basis_global_jacobians[i][0]);
        }

    }

    //! get a reference to the GridView
    const typename Traits::GridView& gridView() const
    {
      return _lfs.gridFunctionSpace().gridView();
    }

    const LFS& localFunctionSpace() const
    {
      return _lfs;
    }

    const std::vector<typename Traits::Range>& basis() const
    {
      return _basis;
    }

    const typename Traits::Range& solution() const
    {
      return _solution;
    }

    const typename Traits::GlobalJacobian& gradient() const
    {
      return _gradient;
    }

  private:

    const Impl& impl() const
    {
      return static_cast<const Impl&>(*this);
    }

    const LFS& _lfs;
    const std::shared_ptr<Data> _data;
    mutable std::vector<typename Traits::Range> _basis;
    mutable typename Traits::Range _solution;
    mutable std::vector<typename Traits::Jacobian> _basis_jacobians;
    mutable std::vector<typename Traits::GlobalJacobian> _basis_global_jacobians;
    mutable typename Traits::GlobalJacobian _gradient;

  };






  template<typename LFS, typename Data, typename Impl, typename Policy>
  class VTKDiscreteFunctionBase
    : public Dune::TypeTree::LeafNode
    , public ScalarDiscreteCalculationProvider<LFS,
                                               Data,
                                               Impl
                                               >
    , public Dune::PDELab::GridFunctionInterface<typename Policy::Traits,
                                                 Impl
                                                 >
  {

    typedef ScalarDiscreteCalculationProvider<
      LFS,
      Data,
      Impl
      > CalculationProviderBase;

    typedef GridFunctionInterface<
      typename Policy::Traits,
      Impl
      > GridFunctionBase;

  public:

    struct Traits
      : public CalculationProviderBase::Traits
      , public GridFunctionBase::Traits
    {

      using CalculationProviderBase::Traits::dimDomain;
      using GridFunctionBase::Traits::dimRange;

    };

    VTKDiscreteFunctionBase(const LFS& lfs, const std::shared_ptr<Data>& data)
      : CalculationProviderBase(lfs,data)
      , GridFunctionBase(lfs.gridFunctionSpace().dataSetType())
    {}

  };

  /*
  template<typename LFS, typename Data, typename Impl>
  class VTKScalarDiscreteFunctionBase:
    public VTKDiscreteFunctionBase<LFS,
                                   Data,
                                   Impl,
                                   ScalarFunctionPolicy<LFS>
                                   >
  {

    typedef VTKDiscreteFunctionBase<
      LFS,
      Data,
      Impl,
      ScalarFunctionPolicy<
        LFS
        >
      > Base;

  public:

    using Base::Base;

  };

  */

  template<typename LFS, typename Data, typename Impl>
  using VTKScalarDiscreteFunctionBase = VTKDiscreteFunctionBase<LFS,
                                                                Data,
                                                                Impl,
                                                                ScalarFunctionPolicy<
                                                                  LFS
                                                                  >
                                                                >;

  template<typename LFS, typename Data, typename Impl>
  using VTKCoordinateDiscreteFunctionBase = VTKDiscreteFunctionBase<LFS,
                                                                    Data,
                                                                    Impl,
                                                                    CoordinateFunctionPolicy<
                                                                      LFS
                                                                      >
                                                                    >;


} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_VTK_HH
