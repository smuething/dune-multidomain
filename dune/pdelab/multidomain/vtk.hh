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

    struct VTKDiscreteFunctionBaseFlagsBase
    {};

    template<bool enable = false>
    struct VTKDiscreteFunctionBaseFlags
      : public std::conditional<enable,
                                VTKDiscreteFunctionBaseFlagsBase,
                                VTKDiscreteFunctionBaseFlags<true>
                                >::type
    {

      typedef VTKDiscreteFunctionBaseFlags<true> Computations;

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

  template<typename LFS, typename Data, typename Impl>
  class VTKDiscreteFunctionBase
    : public Dune::TypeTree::LeafNode
    , public detail::VTKDiscreteFunctionBaseFlags<>
    , public Dune::PDELab::GridFunctionInterface<Dune::PDELab::GridFunctionTraits<
                                                   typename LFS::Traits::GridView,
                                                   typename Dune::BasisInterfaceSwitch<
                                                     typename Dune::FiniteElementInterfaceSwitch<
                                                       typename LFS::Traits::FiniteElement
                                                       >::Basis
                                                     >::RangeField,
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
                                                   >,
                                                 Impl
                                                 >
  {

    typedef Dune::BasisInterfaceSwitch<
      typename Dune::FiniteElementInterfaceSwitch<
        typename LFS::Traits::FiniteElement
        >::Basis
      > BasisSwitch;

    typedef Dune::PDELab::GridFunctionInterface<
      Dune::PDELab::GridFunctionTraits<
        typename LFS::Traits::GridView,
        typename BasisSwitch::RangeField,
        BasisSwitch::dimRange,
        typename BasisSwitch::Range
        >,
      Impl
      > Base;

  public:

    struct Traits
      : public Base::Traits
    {

      typedef typename Base::Traits::GridViewType GridView;

      typedef typename Dune::FiniteElementInterfaceSwitch<
        typename LFS::Traits::FiniteElement
        >::Basis Basis;

      typedef Dune::FieldVector<
        typename Basis::Traits::RangeFieldType,
        GridView::dimensionworld
        > Gradient;

    };

    typedef detail::VTKDiscreteFunctionBaseFlags<>::Computations Computations;

    VTKDiscreteFunctionBase(const LFS& lfs, const std::shared_ptr<Data>& data)
      : Base(lfs.gridFunctionSpace().dataSetType())
      , _lfs(lfs)
      , _data(data)
      , _basis(lfs.maxSize())
      , _basis_jacobians(lfs.maxSize())
      , _basis_gradients(lfs.maxSize())
    {}

    void bind(const typename Traits::ElementType& e,
              const typename Traits::DomainType& x) const
    {

      _data->bind(e);

      typedef Dune::FiniteElementInterfaceSwitch<
        typename LFS::Traits::FiniteElement
        > FESwitch;

      if (impl().evaluateBasis() || impl().evaluateSolution())
        FESwitch::basis(_lfs.finiteElement()).evaluateFunction(x,_basis);

      if (impl().evaluateSolution())
        {
          _solution = typename Traits::RangeType(0);

          for (std::size_t i = 0; i < _lfs.size(); ++i)
            _solution.axpy(_data->_x_local(_lfs,i),_basis[i]);
        }

      if (impl().evaluateBasisJacobians() ||
          impl().evaluateBasisGradients() ||
          impl().evaluateGradient())
        {
          FESwitch::basis(_lfs.finiteElement()).evaluateJacobian(x,_basis_jacobians);
        }

      if (impl().evaluateBasisGradients() ||
          impl().evaluateGradient())
        {
          auto geometry = e.geometry();
          auto JgeoIT = geometry.jacobianInverseTransposed(x);

          for (std::size_t i = 0; i < _lfs.size(); ++i)
            JgeoIT.mv(_basis_jacobians[i][0],_basis_gradients[i]);

        }

      if (impl().evaluateGradient())
        {
          _gradient = typename Traits::Gradient(0);

          for (std::size_t i = 0; i < _lfs.size(); ++i)
            _gradient.axpy(_data->_x_local(_lfs,i),_basis_gradients[i]);
        }

    }

    //! get a reference to the GridView
    const typename Traits::GridViewType& gridView() const
    {
      return _lfs.gridFunctionSpace().gridView();
    }

    const LFS& localFunctionSpace() const
    {
      return _lfs;
    }

    const std::vector<typename Traits::RangeType>& basis() const
    {
      return _basis;
    }

    const typename Traits::RangeType& solution() const
    {
      return _solution;
    }

    const typename Traits::Gradient& gradient() const
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
    mutable std::vector<typename Traits::RangeType> _basis;
    mutable typename Traits::RangeType _solution;
    mutable std::vector<typename Traits::Basis::Traits::JacobianType> _basis_jacobians;
    mutable std::vector<typename Traits::Gradient> _basis_gradients;
    mutable typename Traits::Gradient _gradient;

  };


} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_VTK_HH
