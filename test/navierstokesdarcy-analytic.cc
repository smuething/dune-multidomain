#include "config.h"

#include <dune/common/parametertreeparser.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/p1fem.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/finiteelementmap/q12dfem.hh>
#include <dune/pdelab/finiteelementmap/pk2dfem.hh>
#include <dune/pdelab/multidomain/istlhelpers.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include <dune/pdelab/localoperator/cg_stokes.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/constraints/constraintsparameters.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/multidomain/vtk.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <typeinfo>

#include <dune/pdelab/finiteelementmap/opbfem.hh>
#include "adrwip.hh"

#include "stokesdarcycouplingoperator.hh"

template<typename GV, typename RF>
class ZeroScalarFunction :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
  ZeroScalarFunction<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, ZeroScalarFunction<GV,RF> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  ZeroScalarFunction(const GV & gv) : BaseT(gv) {}

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    y=0;
  }
};

template<typename GV, typename RF>
class PressureFunction :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
  PressureFunction<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, PressureFunction<GV,RF> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  PressureFunction(const GV & gv) : BaseT(gv) {}

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    y=std::exp(x[0]) * std::sin(x[0]+x[1]);
  }
};


template<typename GV, typename RF>
class NavierStokesParameters
{
public:

  static const bool assemble_navier = false;
  static const bool assemble_full_tensor = true;

  typedef Dune::PDELab::NavierStokesParameterTraits<GV,RF> Traits;

  typedef Dune::FieldVector<typename GV::Grid::ctype,GV::Grid::dimensionworld> Vector;

  //! source term
  template<typename EG>
  typename Traits::VelocityRange
  f(const EG& eg, const typename Traits::Domain& coord) const
  {
    auto global_coord = eg.geometry().global(coord);
    auto x = global_coord[0];
    auto y = global_coord[1];
    typename Traits::VelocityRange r;
    /*
    r[0] = _rho*(-y*std::sin(x*y)*std::cos(x*y) - x*std::sin(x*y)*std::exp(x+y))
      + std::exp(x)*std::sin(x+y) + std::cos(x+y)
      - _mu*(-(x*x + y*y)*std::cos(x*y) -y*y*std::cos(x*y) + std::exp(x+y));
    r[1] = _rho*(std::exp(x+y)*std::cos(x*y) + std::exp(x+y)*std::exp(x+y))
      + std::exp(x)*std::cos(x+y)
      - _mu*(2*std::exp(x+y) -std::sin(x*y) - x*y*std::cos(x*y) + std::exp(x+y));
    */
    r[0] = exp(x)*sin(x+y) + exp(x)*cos(x+y) - _mu*(-2*y*y*cos(x*y)-x*x*cos(x*y)+exp(x+y));
    r[1] = exp(x)*cos(x+y) - _mu*(3*exp(x+y) - sin(x*y) - x*y*cos(x*y));
    return r;
  }

  //! source term
  template<typename EG>
  typename Traits::PressureRange
  g2(const EG& eg, const typename Traits::Domain& coord) const
  {
    auto global_coord = eg.geometry().global(coord);
    auto x = global_coord[0];
    auto y = global_coord[1];
    return -y * std::sin(x*y) + std::exp(x+y);
  }

  //! boundary condition type from local intersection coordinate
  template<typename IG>
  typename Traits::BoundaryCondition::Type
  bctype(const IG& is,
         const typename Traits::IntersectionDomain& x) const
  {
    if (!is.boundary())
      return Traits::BoundaryCondition::DoNothing;
    auto cg = is.geometry().global(x);
    RF xg = cg[0];
    RF yg = cg[1];
    return (xg > (-0.5+1e-9)) && (yg > (-0.5+1e-9)) ? Traits::BoundaryCondition::StressNeumann : Traits::BoundaryCondition::VelocityDirichlet;
  }

  //! Dynamic viscosity value from local cell coordinate
  template<typename EG>
  typename Traits::RangeField
  mu(const EG& eg, const typename Traits::Domain& x) const
  {
    return _mu;
  }

  //! Dynamic viscosity value from local intersection coordinate
  template<typename IG>
  typename Traits::RangeField
  mu(const IG& ig, const typename Traits::IntersectionDomain& x) const
  {
    return _mu;
  }

  //! Dynamic viscosity value from local cell coordinate
  template<typename EG>
  typename Traits::RangeField
  rho (const EG& eg, const typename Traits::Domain& x) const
  {
    return _rho;
  }

  //! Dynamic viscosity value from local intersection coordinate
  template<typename IG>
  typename Traits::RangeField
  rho (const IG& ig, const typename Traits::IntersectionDomain& x) const
  {
    return _rho;
  }

  //! Dirichlet boundary condition value from local cell coordinate
  template<typename EG>
  typename Traits::VelocityRange
  g(const EG& eg, const typename Traits::Domain& coord) const
  {
    auto global_coord = eg.geometry().global(coord);
    auto x = global_coord[0];
    auto y = global_coord[1];
    typename Traits::VelocityRange r;
    r[0] = std::cos(x*y);
    r[1] = std::exp(x+y);
    return r;
  }

  //! Dirichlet boundary condition value from local intersection coordinate
  template<typename IG>
  typename Traits::VelocityRange
  g(const IG& ig, const typename Traits::IntersectionDomain& coord) const
  {
    auto global_coord = ig.geometry().global(coord);
    auto x = global_coord[0];
    auto y = global_coord[1];
    typename Traits::VelocityRange r;
    r[0] = std::cos(x*y);
    r[1] = std::exp(x+y);
    return r;
  }

  //! Neumann boundary condition (stress) - version for vector-valued function
  template<typename IG>
  typename Traits::VelocityRange
  j(const IG& ig,
    const typename Traits::IntersectionDomain& coord,
    const typename Traits::Domain& normal) const
  {
    Dune::FieldMatrix<RF,Traits::Domain::dimension,Traits::Domain::dimension> s(0.0);
    auto global_coord = ig.geometry().global(coord);
    auto x = global_coord[0];
    auto y = global_coord[1];
    s[0][0] = s[1][1] = std::exp(x) * std::sin(x+y);
    s[0][0] += 2 * _mu * y * std::sin(x*y);
    s[0][1] += _mu * (x * std::sin(x*y) - std::exp(x+y));
    s[1][0] += _mu * (x * std::sin(x*y) - std::exp(x+y));
    s[1][1] -= 2 * _mu * std::exp(x+y);
    typename Traits::VelocityRange r(0);
    s.mv(normal,r);
    return r;
  }

  NavierStokesParameters(const Dune::ParameterTree& params)
    : _mu(params.get<RF>("fluid.mu"))
    , _rho(params.get<RF>("fluid.rho"))
  {}

private:
  const RF _mu;
  const RF _rho;
};


template<typename GV, typename P>
class VelocityFunctionAdapter
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<
        GV,
        typename P::Traits::RangeField,
        P::Traits::dimDomain,
        typename P::Traits::VelocityRange
        >,
      VelocityFunctionAdapter<
        GV,
        P
        >
    >
{
public:

  typedef Dune::PDELab::GridFunctionTraits<
    GV,
    typename P::Traits::RangeField,
    P::Traits::dimDomain,
    typename P::Traits::VelocityRange
    > Traits;

  typedef GV GridViewType;

  VelocityFunctionAdapter(const P& params)
    : _params(params)
  {}

  void evaluate(const typename Traits::ElementType& e,
                const typename Traits::DomainType& x,
                typename Traits::RangeType& y) const
  {
    y = _params.g(e,x);
  }

private:

  const P& _params;

};


template<typename GV, typename RF>
class DarcyParameters
{

public:
  typedef Dune::PM::ADRTraits<GV,RF> Traits;
  typedef Dune::PM::ADRBoundaryConditions BC;


  //! constructor
  DarcyParameters(const Dune::ParameterTree& params)
    : _mu(params.get<RF>("parameters.fluid.mu"))
  {
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        {
          _kabs[i][j] = (i==j) ? params.get<RF>("parameters.soil.permeability") : 0;
        }
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x_) const
  {
    return typename Traits::RangeType(0.0);
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  K (const typename Traits::ElementType& e, const typename Traits::DomainType& x_) const
  {
    return _kabs;
  }

  //! nonlinear scalar diffusion term
  typename Traits::RangeFieldType
  w (const typename Traits::ElementType& e, const typename Traits::DomainType& x, typename Traits::RangeFieldType u) const
  {
    return 1.0;
  }

  typename Traits::RangeFieldType
  porosity (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 1.0;
  }

  //! reaction term
  typename Traits::RangeFieldType
  a (const typename Traits::ElementType& e, const typename Traits::DomainType& x_) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x_) const
  {
    auto global_coord = e.geometry().global(x_);
    auto x = global_coord[0];
    auto y = global_coord[1];
    return _kabs[0][0] / _mu * std::exp(x)*(std::sin(x+y) - 2*std::cos(x+y));
  }

  //! boundary condition type function
  typename BC::Type
  bc (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x_) const
  {
    if (!is.boundary())
      return BC::None;
    else {
      Dune::FieldVector<typename GV::Grid::ctype,GV::dimension>
        x = is.geometry().global(x_);
      if (x[0] < 1e-9)
        return BC::Dirichlet;
      else
        return BC::Dirichlet;
    }
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x_) const
  {
    auto global_coord = is.geometry().global(x_);
    auto x = global_coord[0];
    auto y = global_coord[1];
    return std::exp(x) * std::sin(x+y);
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x_) const
  {
    auto global_coord = is.geometry().global(x_);
    auto x = global_coord[0];
    auto y = global_coord[1];

    auto normal = is.unitOuterNormal(x_);

    return -normal[0] * (y*std::sin(x*y) * normal[0] + x*std::sin(x*y) * normal[1])
      +  normal[1] * (normal[0] + normal[1]) * std::exp(x+y);
  }

  //! set time for subsequent evaluation
  void setTime (typename Traits::RangeFieldType t)
  {
  }

private:
  typename Traits::PermTensorType _kabs;
  const typename Traits::RangeFieldType _mu;
};


template<typename GV, typename RF>
class StokesDarcyCouplingParameters
{

public:

  typedef StokesDarcyCouplingParameterTraits<GV,RF> Traits;

  template<typename IG>
  typename Traits::RangeField h1(const IG& ig,
                                 const typename Traits::IntersectionDomain& coord,
                                 const typename Traits::Domain& normal) const
  {
    auto global_coord = ig.geometry().global(coord);
    auto x = global_coord[0];
    auto y = global_coord[1];
    Dune::FieldMatrix<RF,GV::dimension,GV::dimension> s(0.0);
    s[0][0] = 2 * _mu * y * std::sin(x*y);
    s[1][1] = - 2 * _mu * std::exp(x+y);
    s[0][1] = s[1][0] = _mu * (x * std::sin(x*y) - std::exp(x+y));
    typename Traits::Domain ps(0.0);
    s.mv(normal,ps);
    return ps * normal;
  }

  template<typename IG>
  typename Traits::VelocityRange h2(const IG& ig,
                                    const typename Traits::IntersectionDomain& coord,
                                    const typename Traits::Domain& normal) const
  {
    typedef typename Traits::VelocityRange R;
    auto global_coord = ig.geometry().global(coord);
    auto x = global_coord[0];
    auto y = global_coord[1];
    Dune::FieldMatrix<RF,GV::dimension,GV::dimension> s(0.0);
    s[0][0] = s[1][1] = std::exp(x) * std::sin(x+y);
    s[0][1] = s[1][0] = _mu * (x * std::sin(x*y) - std::exp(x+y));
    s[0][0] += 2 * _mu * y * std::sin(x*y);
    s[1][1] -= 2 * _mu * std::exp(x+y);

    R normal_stress(0.0);
    s.mv(normal,normal_stress);
    R normal_stress_normal(normal);
    normal_stress_normal *= normal * normal_stress;
    R normal_stress_tangent(normal_stress);
    normal_stress_tangent -= normal_stress_normal;

    R u(0.0);
    u[0] = std::cos(x*y);
    u[1] = std::exp(x+y);

    R normal_flow(normal);
    normal_flow *= normal * u;

    R tangent_flow(u);
    tangent_flow -= normal_flow;
    tangent_flow *= _alpha/sqrt((_K[0][0] + _K[1][1])/2.0);

    normal_stress_tangent -= tangent_flow;

    return normal_stress_tangent;
  }

  template<typename IG>
  typename Traits::RangeField h3(const IG& ig,
                                 const typename Traits::IntersectionDomain& coord,
                                 const typename Traits::Domain& normal) const
  {
    auto global_coord = ig.geometry().global(coord);
    auto x = global_coord[0];
    auto y = global_coord[1];

    typename Traits::Domain u(0.0);
    u[0] = _K[0][0] * std::exp(x) * (std::sin(x+y) + std::cos(x+y)) - std::cos(x*y);
    u[1] = _K[0][0] * std::exp(x) * std::cos(x+y) - std::exp(x+y);

    return (u * normal);
  }

  template<typename IG>
  typename Traits::RangeField alpha(const IG& ig, const typename Traits::IntersectionDomain& coord) const
  {
    return _alpha;
  }

  template<typename IG>
  typename Traits::RangeField viscosity(const IG& ig, const typename Traits::IntersectionDomain& coord) const
  {
    return 1.0;
  }

  template<typename IG>
  typename Traits::RangeField porosity(const IG& ig, const typename Traits::IntersectionDomain& coord) const
  {
    return 1.0;
  }

  template<typename IG>
  typename Traits::RangeField gamma(const IG& ig, const typename Traits::IntersectionDomain& coord) const
  {
    return 1.0;
  }

  template<typename IG>
  typename Traits::RangeField density(const IG& ig, const typename Traits::IntersectionDomain& coord) const
  {
    return 1.0;
  }

  template<typename IG>
  typename Traits::PermeabilityTensor
  K (const IG& ig, const typename Traits::IntersectionDomain& x_) const
  {
    return _K;
  }

  StokesDarcyCouplingParameters(const Dune::ParameterTree& params)
    : _mu(params.get<RF>("fluid.mu"))
    , _alpha(params.get<RF>("interface.alpha"))
    , _epsilon(params.get("epsilon",1e-8))
  {
    for (std::size_t i = 0; i < Traits::dimDomain; ++i)
      for (std::size_t j = 0; j < Traits::dimDomain; ++j)
        _K[i][j] = (i == j ? params.get<RF>("soil.permeability") : 0);
  }

  double epsilon() const
  {
    return _epsilon;
  }

private:

  const RF _mu;
  typename Traits::PermeabilityTensor _K;
  const RF _alpha;
  const double _epsilon;
};

template<typename DarcyParameters>
struct DarcyBoundaryTypeAdapter :
  public Dune::PDELab::DirichletConstraintsParameters
{

  template<typename I>
  bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    return parameters.bc(ig,x) == DarcyParameters::BC::Dirichlet;
  }

  template<typename I>
  bool isNeumann(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    typename DarcyParameters::BC::Type bc = parameters.bc(ig,x);
    return
      bc == DarcyParameters::BC::Flux ||
      bc == DarcyParameters::BC::Outflow;
  }

  template<typename R>
  void setTime(const R& r)
  {
    parameters.setTime(r);
  }

  DarcyBoundaryTypeAdapter(DarcyParameters& params)
    : parameters(params)
  {}

  DarcyParameters& parameters;
};



template<typename DarcyParameter, typename GV = typename DarcyParameter::Traits::GridViewType>
struct DarcyFlowFromPressure
{

  template<typename LFS, typename Data>
  class Function
    : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<
                                              GV,
                                              typename DarcyParameter::Traits::RangeFieldType,
                                              DarcyParameter::Traits::dimDomain,
                                              typename DarcyParameter::Traits::RangeType
                                              >,
                                            Function<LFS,Data> >
  {

    const DarcyParameter& p;
    const LFS& lfs;
    std::shared_ptr<Data> data;
    mutable std::vector<
      Dune::FieldMatrix<
        typename DarcyParameter::Traits::RangeFieldType,
        1,
        DarcyParameter::Traits::dimDomain
        >
      > _gradient;

  public:

    Function(const DarcyParameter& p_,
             const LFS& lfs_,
             std::shared_ptr<Data> data_)
      : p(p_)
      , lfs(lfs_)
      , data(data_)
      , _gradient(lfs_.maxSize())
    {}

    typedef Dune::PDELab::GridFunctionTraits<
      GV,
      typename DarcyParameter::Traits::RangeFieldType,
      DarcyParameter::Traits::dimDomain,
      typename DarcyParameter::Traits::RangeType
      > Traits;

    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& x,
                          typename Traits::RangeType& y) const
    {

      data->bind(e);

      typedef Dune::FiniteElementInterfaceSwitch<
        typename LFS::Traits::FiniteElement
        > FESwitch;

      typedef Dune::BasisInterfaceSwitch<
        typename FESwitch::Basis
        > BasisSwitch;

      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;

      BasisSwitch::gradient(
        FESwitch::basis(lfs.finiteElement()),
        e.geometry(),
        x,
        _gradient
      );

      typename Traits::RangeType gradp(0.0);

      for (std::size_t i = 0; i < lfs.size(); ++i)
        gradp.axpy(data->_x_local(lfs,i),_gradient[i][0]);

      p.K(Traits::GridViewType::Grid::multiDomainEntity(e),x).mv(gradp,y);

      y *= -1;
    }

    const typename Traits::GridViewType& getGridView ()
    {
      return lfs.gridFunctionSpace().gridview();
    }
  };

  template<typename LFS, typename Data>
  struct create_type
  {
    typedef Function<LFS,Data> type;
  };

  template<typename LFS, typename Data>
  std::shared_ptr<typename create_type<LFS,Data>::type> create(const LFS& lfs, std::shared_ptr<Data> data)
  {
    typedef typename create_type<LFS,Data>::type F;
    return std::make_shared<F>(p,lfs,data);
  }

  const DarcyParameter& p;

  DarcyFlowFromPressure(const DarcyParameter& param)
    : p(param)
  {}

};

namespace Dune {
namespace PDELab {

template<typename T, typename Params, typename X>
class DiscretePressureGridFunction
  : public GridFunctionInterface<
  GridFunctionTraits<
    typename T::Traits::GridViewType,
    typename T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
    T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimRange,
    typename T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
    >,
  DiscretePressureGridFunction<T,Params,X>
  >
{
  typedef T GFS;

  typedef GridFunctionInterface<
    GridFunctionTraits<
      typename T::Traits::GridViewType,
      typename T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
      T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimRange,
      typename T::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
      >,
    DiscretePressureGridFunction<T,Params,X>
    > BaseT;

public:
  typedef typename BaseT::Traits Traits;

  /** \brief Construct a DiscreteGridFunction
   *
   * \param gfs The GridFunctionsSpace
   * \param x_  The coefficients vector
   */
  DiscretePressureGridFunction (const GFS& gfs, const Params& params, const X& x_, const typename Traits::RangeFieldType& factor_)
    : pgfs(Dune::stackobject_to_shared_ptr(gfs)), parameters(params), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), yb(gfs.maxLocalSize()), factor(factor_)
  {
  }

  // Evaluate
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    lfs.bind(e);
    lfs.vread(xg,xl);
    lfs.finiteElement().localBasis().evaluateFunction(x,yb);
    typename Traits::RangeType t = 0;
    for (unsigned int i=0; i<yb.size(); i++)
      t.axpy(xl[i],yb[i]);
    typename Traits::ElementType::Geometry::GlobalCoordinate xg = e.geometry().global(x);
    t -= xg[Traits::dimDomain-1];
    y = t;
  }

  //! get a reference to the GridView
  inline const typename Traits::GridViewType& getGridView () const
  {
    return pgfs->gridView();
  }

private:
  shared_ptr<GFS const> pgfs;
  const Params& parameters;
  const X& xg;
  mutable Dune::PDELab::LocalFunctionSpace<GFS> lfs;
  mutable std::vector<typename Traits::RangeFieldType> xl;
  mutable std::vector<typename Traits::RangeType> yb;
  const typename Traits::RangeFieldType factor;
};

template<typename T, typename Params, typename X>
class DarcyFlowFromPotential
  : public GridFunctionInterface<
  GridFunctionTraits<
    typename T::Traits::GridViewType,
    typename T::Traits::FiniteElementType::Traits::LocalBasisType
    ::Traits::RangeFieldType,
    T::Traits::FiniteElementType::Traits::LocalBasisType::Traits
    ::dimDomain,
    FieldVector<
      typename T::Traits::FiniteElementType::Traits
      ::LocalBasisType::Traits::RangeFieldType,
      T::Traits::FiniteElementType::Traits::LocalBasisType::Traits
      ::dimDomain> >,
  DarcyFlowFromPotential<T,Params,X> >
{
  typedef T GFS;
  typedef typename GFS::Traits::FiniteElementType::Traits::
  LocalBasisType::Traits LBTraits;

public:
  typedef GridFunctionTraits<
  typename GFS::Traits::GridViewType,
  typename LBTraits::RangeFieldType,
  LBTraits::dimDomain,
  FieldVector<
    typename LBTraits::RangeFieldType,
    LBTraits::dimDomain> > Traits;

private:
  typedef GridFunctionInterface<
  Traits,
  DarcyFlowFromPotential<T,Params,X> > BaseT;

public:
  /** \brief Construct a DiscreteGridFunctionGradient
   *
   * \param gfs The GridFunctionsSpace
   * \param x_  The coefficients vector
   */
  DarcyFlowFromPotential (const GFS& gfs, const Params& params, const X& x_)
    : pgfs(Dune::stackobject_to_shared_ptr(gfs)), parameters(params), xg(x_)
  { }

  // Evaluate
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    // get and bind local functions space
    Dune::PDELab::LocalFunctionSpace<GFS> lfs(*pgfs);
    lfs.bind(e);

    // get local coefficients
    std::vector<typename Traits::RangeFieldType> xl(lfs.size());
    lfs.vread(xg,xl);

    // get Jacobian of geometry
    const typename Traits::ElementType::Geometry::Jacobian&
      JgeoIT = e.geometry().jacobianInverseTransposed(x);

    // get local Jacobians/gradients of the shape functions
    std::vector<typename LBTraits::JacobianType> J(lfs.size());
    lfs.finiteElement().localBasis().evaluateJacobian(x,J);

    typename Traits::RangeType gradphi, t;
    t = 0;
    for(unsigned i = 0; i < lfs.size(); ++i) {
      // compute global gradient of shape function i
      gradphi = 0;
      JgeoIT.umv(J[i][0], gradphi);

      // sum up global gradients, weighting them with the appropriate coeff
      t.axpy(xl[i], gradphi);
    }
    t /= -parameters.porosity(getGridView().grid().multiDomainEntity(e),x);
    parameters.K(getGridView().grid().multiDomainEntity(e),x).mv(t,y);
  }

  //! get a reference to the GridView
  inline const typename Traits::GridViewType& getGridView () const
  {
    return pgfs->gridView();
  }

private:
  shared_ptr<GFS const> pgfs;
  const Params& parameters;
  const X& xg;
};

} // namespace PDELab
} // namespace Dune

int main(int argc, char** argv) {

  try {

    Dune::MPIHelper::instance(argc,argv);

    Dune::ParameterTree parameters;
    Dune::ParameterTreeParser::readINITree(argv[1],parameters);

    const int dim = 2;

    std::vector<int> boundaryIndexToPhysicalGroup, elementIndexToPhysicalGroup;

    typedef Dune::UGGrid<dim> BaseGrid;

    std::shared_ptr<BaseGrid> baseGridPtr(
      Dune::GmshReader<BaseGrid>::read(
        parameters["mesh.filename"],
        boundaryIndexToPhysicalGroup,
        elementIndexToPhysicalGroup,
        true,
        false
      )
    );

    BaseGrid& baseGrid = *baseGridPtr;

    typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,4> > Grid;
    Grid grid(baseGrid,false);
    typedef Grid::SubDomainGrid SubDomainGrid;
    SubDomainGrid& stokesGrid = grid.subDomain(0);
    SubDomainGrid& darcyGrid = grid.subDomain(1);
    typedef Grid::ctype ctype;
    typedef Grid::LeafGridView MDGV;
    typedef SubDomainGrid::LeafGridView SDGV;
    MDGV mdgv = grid.leafView();
    SDGV stokesGV = stokesGrid.leafView();
    SDGV darcyGV = darcyGrid.leafView();

    typedef MDGV::Grid::ctype DF;
    typedef double RF;

    const bool interface_by_subdomain = parameters.get<bool>("mesh.interface-by-subdomain");
    const DF interface_position = parameters.get("mesh.interface-position",0.0);
    const int stokes_domain_id = parameters.get("mesh.stokes-domain-id",0);
    const int darcy_domain_id = parameters.get("mesh.darcy-domain-id",1);
    assert(stokes_domain_id != darcy_domain_id && "Stokes and Darcy domain must have distinct ids");

    grid.startSubDomainMarking();
    for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
      {
        if (interface_by_subdomain)
          {
            int domain = elementIndexToPhysicalGroup[mdgv.indexSet().index(*it)];
            if (domain == stokes_domain_id)
              grid.addToSubDomain(0,*it);
            if (domain == darcy_domain_id)
              grid.addToSubDomain(1,*it);
          }
        else
          {
            auto x = it->geometry().center()[0];
            grid.addToSubDomain(x > interface_position ? 0 : 1,*it);
          }
      }
    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();

    grid.globalRefine(parameters.get<int>("mesh.refine"));


    const int stokes_k = 2;
    const int stokes_q = 2*stokes_k;
    const int darcy_k = 3;

    typedef Dune::PDELab::Pk2DLocalFiniteElementMap<SDGV,DF,RF,stokes_k> V_FEM;
    typedef Dune::PDELab::Pk2DLocalFiniteElementMap<SDGV,DF,RF,stokes_k-1> P_FEM;
    typedef Dune::PDELab::OPBLocalFiniteElementMap<DF,RF,darcy_k,dim,Dune::GeometryType::simplex,Dune::GMPField<512> > DarcyFEM;
    typedef Dune::PDELab::P1LocalFiniteElementMap<DF,RF,dim-1> CouplingFEM;

    V_FEM vfem(stokesGV);
    P_FEM pfem(stokesGV);
    DarcyFEM darcyfem;
    CouplingFEM couplingfem;

    typedef Dune::PDELab::NoConstraints NOCON;
    typedef Dune::PDELab::ConformingDirichletConstraints DCON;

    typedef Dune::PDELab::ISTLVectorBackend<> VectorBackend;
    typedef Dune::PDELab::ISTLVectorBackend<> VelocityVectorBackend;
    typedef Dune::PDELab::ISTLVectorBackend<> TaylorHoodVectorBackend;

    typedef Dune::PDELab::LexicographicOrderingTag VelocityOrderingTag;
    typedef Dune::PDELab::LexicographicOrderingTag TaylorHoodOrderingTag;

    typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> EC;
    typedef Dune::PDELab::MultiDomain::SubDomainSupersetCondition<Grid> SC;

    EC c0(0);
    EC c1(1);

    typedef Dune::PDELab::VectorGridFunctionSpace<
      SDGV,V_FEM,dim,
      VelocityVectorBackend,
      VectorBackend,
      DCON,
      VelocityOrderingTag
      > PGFS_V_GFS;
    PGFS_V_GFS powervgfs(stokesGV,vfem);
    powervgfs.name("v");

    typedef Dune::PDELab::GridFunctionSpace<
      SDGV,P_FEM,
      DCON,
      VectorBackend
      > P_GFS;
    P_GFS pgfs(stokesGV,pfem);
    pgfs.name("p");

    typedef Dune::PDELab::CompositeGridFunctionSpace<
      TaylorHoodVectorBackend,
      TaylorHoodOrderingTag,
      PGFS_V_GFS,
      P_GFS
      > StokesGFS;
    StokesGFS stokesgfs(powervgfs,pgfs);

    typedef Dune::PDELab::GridFunctionSpace<
      SDGV,DarcyFEM,
      NOCON,
      VectorBackend
      > DarcyGFS;
    DarcyGFS darcygfs(darcyGV,darcyfem);
    darcygfs.name("p");

    typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<
      Grid,
      VectorBackend,
      StokesGFS,
      DarcyGFS
      > MultiGFS;
    MultiGFS multigfs(grid,stokesgfs,darcygfs);

    typedef MultiGFS::ConstraintsContainer<RF>::Type C;
    C cg;

    typedef NavierStokesParameters<MDGV,RF> NavierStokesParams;
    NavierStokesParams navierStokesParams(parameters.sub("parameters"));

    typedef Dune::PDELab::TaylorHoodNavierStokes<NavierStokesParams> StokesOperator;
    StokesOperator stokesOperator(navierStokesParams,stokes_q);

    typedef DarcyParameters<MDGV,RF> DarcyParams;
    DarcyParams darcyParams(parameters);

    typedef DarcyBoundaryTypeAdapter<DarcyParams> DarcyBoundaryFunction;
    DarcyBoundaryFunction darcyBoundaryFunction(darcyParams);

    typedef Dune::PM::ADRWIP<DarcyParams> DarcyOperator;
    DarcyOperator darcyOperator(darcyParams,Dune::PM::ADRWIPMethod::OBB,Dune::PM::ADRWIPWeights::weightsOn,0.0);

    typedef StokesDarcyCouplingParameters<MDGV,RF> CouplingParams;
    CouplingParams couplingParams(parameters.sub("parameters"));

    typedef StokesDarcyCouplingOperator<CouplingParams> CouplingOperator;
    CouplingOperator couplingOperator(couplingParams);

    typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,MultiGFS,StokesOperator,EC,StokesGFS> StokesSubProblem;

    typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,MultiGFS,DarcyOperator,EC,DarcyGFS> DarcySubProblem;

    StokesSubProblem stokesSubProblem(stokesOperator,c0);
    DarcySubProblem darcySubProblem(darcyOperator,c1);

    //typedef Dune::PDELab::MultiDomain::EnrichedCoupling<StokesSubProblem,DarcySubProblem,CouplingOperator,2> Coupling;
    typedef Dune::PDELab::MultiDomain::Coupling<StokesSubProblem,DarcySubProblem,CouplingOperator> Coupling;
    Coupling coupling(stokesSubProblem,darcySubProblem,couplingOperator);


    typedef VelocityFunctionAdapter<MDGV,NavierStokesParams> VelocityInitialFunction;
    VelocityInitialFunction velocityInitialFunction(navierStokesParams);

    //typedef ZeroScalarFunction<MDGV,RF> PressureInitialFunction;
    typedef PressureFunction<MDGV,RF> PressureInitialFunction;
    PressureInitialFunction pressureInitialFunction(mdgv);

    typedef Dune::PDELab::CompositeGridFunction<
      VelocityInitialFunction,
      PressureInitialFunction
      > StokesInitialFunction;
    StokesInitialFunction stokesInitialFunction(velocityInitialFunction,pressureInitialFunction);

    typedef Dune::PDELab::StokesVelocityDirichletConstraints<NavierStokesParams> VelocityDirichletPredicate;
    VelocityDirichletPredicate velocityDirichletPredicate(navierStokesParams);

    typedef Dune::PDELab::NoDirichletConstraintsParameters NoConstraints;
    NoConstraints noConstraints;

    typedef Dune::PDELab::CompositeConstraintsParameters<
      VelocityDirichletPredicate,
      NoConstraints
      > StokesDirichletPredicate;
    StokesDirichletPredicate stokesDirichletPredicate(velocityDirichletPredicate,noConstraints);

    auto constraints = Dune::PDELab::MultiDomain::constraints<RF>(
      multigfs,
      Dune::PDELab::MultiDomain::constrainSubProblem(
        stokesSubProblem,
        velocityDirichletPredicate
      )
    );

    constraints.assemble(cg);

    std::cout << multigfs.size() << " DOF, " << cg.size() << " restricted" << std::endl;

    typedef Dune::PDELab::ISTLMatrixBackend MBE;

    typedef Dune::PDELab::MultiDomain::GridOperator<
      MultiGFS,MultiGFS,
      MBE,RF,RF,RF,C,C,
      StokesSubProblem,
      DarcySubProblem,
      Coupling
      > GridOperator;

    GridOperator gridoperator(multigfs,multigfs,cg,cg,stokesSubProblem,darcySubProblem,coupling);

    typedef GridOperator::Traits::Domain V;

    V u(multigfs,0.0);

    Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,u,stokesInitialFunction,stokesSubProblem);
    Dune::PDELab::set_nonconstrained_dofs(multigfs,cg,0,u);

    V u0(u);

    //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
    //LS ls(5000,true);
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
    LS ls(true);

    //typedef Dune::PDELab::StationaryLinearProblemSolver<GridOperator,LS,V> PDESOLVER;
    //PDESOLVER pdesolver(gridoperator,ls,1e-10);


    typedef Dune::PDELab::Newton<GridOperator,LS,V> PDESOLVER;
    PDESOLVER pdesolver(gridoperator,ls);
    pdesolver.setReassembleThreshold(0.0);
    pdesolver.setVerbosityLevel(2);
    pdesolver.setReduction(1e-10);
    pdesolver.setMinLinearReduction(1e-4);
    pdesolver.setMaxIterations(25);
    pdesolver.setLineSearchMaxIterations(10);


    pdesolver.apply(u);

    PGFS_V_GFS darcy_v_gfs(darcyGV,vfem);
    darcy_v_gfs.name("v_cg");
    typename Dune::PDELab::BackendVectorSelector<PGFS_V_GFS,RF>::Type darcy_v(darcy_v_gfs,0.0);

    {
      auto data = Dune::make_shared<Dune::PDELab::MultiDomain::DGFTreeCommonData<MultiGFS,V,Dune::PDELab::default_predicate> >(multigfs,u);

      auto darcy_flow_function = DarcyFlowFromPressure<DarcyParams,SDGV>(darcyParams).create(data->_lfs.child<1>(),data);

      Dune::PDELab::interpolate(*darcy_flow_function,darcy_v_gfs,darcy_v);
    }

    /*
     * Output initial guess and solution
     */

    {
      Dune::SubsamplingVTKWriter<SDGV> vtkwriter(stokesGV,parameters.get("mesh.refineoutput",0));
      Dune::PDELab::MultiDomain::add_solution_to_vtk_writer(
        vtkwriter,
        multigfs,
        u,
        Dune::PDELab::MultiDomain::subdomain_predicate<Grid::SubDomainIndex>(stokesGV.grid().domain())
      );
      Dune::PDELab::MultiDomain::add_solution_to_vtk_writer(
        vtkwriter,
        multigfs,
        u0,
        Dune::PDELab::MultiDomain::subdomain_predicate<Grid::SubDomainIndex>(stokesGV.grid().domain()),
        Dune::PDELab::default_vtk_name_scheme().prefix("initial-")
      );

      vtkwriter.write(parameters["io.stokesfile"],Dune::VTKOptions::binaryappended);
    }

    {
      Dune::SubsamplingVTKWriter<SDGV> vtkwriter(darcyGV,parameters.get("mesh.refineoutput",0));
      Dune::PDELab::MultiDomain::add_solution_to_vtk_writer(
        vtkwriter,
        multigfs,
        u,
        Dune::PDELab::MultiDomain::subdomain_predicate<Grid::SubDomainIndex>(darcyGV.grid().domain())
      ).add_vertex_function(
                            DarcyFlowFromPressure<DarcyParams,SDGV>(darcyParams),
        Dune::PDELab::TypeTree::TreePath<1>(),
        "v"
      );
      Dune::PDELab::add_solution_to_vtk_writer(
        vtkwriter,
        darcy_v_gfs,
        darcy_v
      );
      Dune::PDELab::MultiDomain::add_solution_to_vtk_writer(
        vtkwriter,
        multigfs,
        u0,
        Dune::PDELab::MultiDomain::subdomain_predicate<Grid::SubDomainIndex>(darcyGV.grid().domain()),
        Dune::PDELab::default_vtk_name_scheme().prefix("initial-")
      );

      vtkwriter.write(parameters["io.darcyfile"],Dune::VTKOptions::binaryappended);
    }


  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
	return 1;
  }
  catch (std::exception &e){
    std::cerr << "STL reported error: " << e.what() << std::endl;
	return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
	return 1;
  }

}
