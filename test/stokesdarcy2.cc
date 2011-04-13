#include "config.h"

#include <dune/common/parametertreeparser.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/p1fem.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/finiteelementmap/q12dfem.hh>
#include <dune/pdelab/finiteelementmap/pk2dfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinggridfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinglocalfunctionspace.hh>
#include <dune/pdelab/multidomain/multidomaingridoperatorspace.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include <dune/pdelab/localoperator/cg_stokes.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include<typeinfo>

#include "../../dune-pm/dune/pm/finiteelements/opbfem.hh"
#include "../../dune-pm/dune/pm/models/adrwip.hh"

#include "stokesdarcycouplingoperator.hh"
#include "functionmacros.hh"

#define UGGRID

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
class RightScalarFunction :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
  RightScalarFunction<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, RightScalarFunction<GV,RF> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  RightScalarFunction(const GV & gv, double v_) : BaseT(gv), v(v_) {}

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    if (x[0] > 1-1e-6)
      y=v;
    else
      y=v;
  }

private:
  const RangeType v;

};


template<typename GV>
class ScalarNeumannBoundaryFunction :
  public Dune::PDELab::BoundaryGridFunctionBase<
  Dune::PDELab::BoundaryGridFunctionTraits<
    GV,int,1,Dune::FieldVector<int,1> >,
  ScalarNeumannBoundaryFunction<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,ScalarNeumannBoundaryFunction<GV> > BaseT;

    ScalarNeumannBoundaryFunction (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0; // Neumann
  }

  //! get a reference to the GridView
  inline const GV& getGridView ()
  {
    return gv;
  }
};

template<typename GV>
class PressureDropBoundaryFunction :
  public Dune::PDELab::BoundaryGridFunctionBase<
  Dune::PDELab::BoundaryGridFunctionTraits<
    GV,int,1,Dune::FieldVector<int,1> >,
  PressureDropBoundaryFunction<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,PressureDropBoundaryFunction<GV> > BaseT;

    PressureDropBoundaryFunction (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    if (!ig.boundary())
      y = -1; // no bc
    else {
      auto xg = ig.geometry().global(x);
      if (xg[0] < 1e-6 || xg[0] > 100-1e-6)
        y = 0; // Neumann
      else
        y = 1; // Dirichlet
    }
  }

  //! get a reference to the GridView
  inline const GV& getGridView ()
  {
    return gv;
  }
};


template<typename GV, typename RF>
class PressureDropFlux
  : public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
  PressureDropFlux<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,PressureDropFlux<GV,RF> > BaseT;

private:
  typedef typename Traits::DomainFieldType DFT;

  const DFT inflow, outflow;

public:
  PressureDropFlux (const GV& gv, const RF inflow_, const RF outflow_)
    : BaseT(gv), inflow(inflow_), outflow(outflow_)
  {
  }

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    if (x[0] < 1e-6)
      y = inflow;//*0.04*(100-x[1])*(x[1]-110); // parabolic flow profile
    else if (x[0] > 100-1e-6)
      y = -outflow;//*0.04*(100-x[1])*(x[1]-110); // parabolic flow profile
    else
      y = 0;//-10*((x[1]-0.5)*(1.5-x[1])+0.02);
  }
};


template<typename GV, typename RF>
class NavierStokesParameters
  : public Dune::PDELab::TaylorHoodNavierStokesParameters<double>
{
public:

  typedef Dune::FieldVector<typename GV::Grid::ctype,GV::Grid::dimensionworld> Vector;

  double rho() const{
    return _viscosity;
  }

  double mu() const{
    return _density;
  }

  Vector f1() const{
    return _gravity;
  }

  NavierStokesParameters(const Dune::ParameterTree& params)
    : _gravity(0.0)
    , _density(params.get<RF>("fluid.density"))
    , _viscosity(params.get<RF>("fluid.viscosity"))
  {
    _gravity[GV::Grid::dimensionworld-1] = - params.get("gravity",9.81);
  }

private:
  Dune::FieldVector<RF,GV::Grid::dimensionworld> _gravity;
  const RF _density;
  const RF _viscosity;
};


template<typename GV, typename RF>
class DarcyParameters
{
  typedef Dune::PM::ADRBoundaryConditions BC;

public:
  typedef Dune::PM::ADRTraits<GV,RF> Traits;

  //! constructor
  DarcyParameters(GV gridview,
                  const std::vector<int>& physicalGroupMap,
                  const Dune::ParameterTree& params)
    : gamma(params.get("parameters.coupling.gamma",1.0))
    , porosity_(params.get<RF>("parameters.soil.porosity"))
    , elementIndexToPhysicalGroup(physicalGroupMap)
    , gv(gridview)
    , bottomPotential(params.get<RF>("boundaries.bottompotential"))
    , bottomflux((params.get<RF>("boundaries.outflow") - params.get<RF>("boundaries.inflow"))*0.1)
  {
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        {
          kabslow[i][j] = (i==j) ? params.get<RF>("parameters.soil.permeability.low") : 0;
          kabs[i][j] = (i==j) ? params.get<RF>("parameters.soil.permeability") : 0;
        }
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x_) const
  {
    typename Traits::RangeType velocity(0.0);
    return velocity;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  K (const typename Traits::ElementType& e, const typename Traits::DomainType& x_) const
  {
    return elementIndexToPhysicalGroup[gv.indexSet().index(e)] == 1 ? kabs : kabslow;
  }

  //! nonlinear scalar diffusion term
  typename Traits::RangeFieldType
  w (const typename Traits::ElementType& e, const typename Traits::DomainType& x, typename Traits::RangeFieldType u) const
  {
    return gamma;
  }

  typename Traits::RangeFieldType
  porosity (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return porosity_;
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
    return 0.0;
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
      if (x[0] < 1e-6 || x[0] > 100-1e-6)
        return BC::Flux;
      else
        return BC::Flux;//BC::Dirichlet;
    }
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x_) const
  {
    return bottomPotential;
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x_) const
  {
    auto x = is.geometry().global(x_);
    if (x[1] < 1e-6)
      return 0.0; //bottomflux;
    else
      return 0.0;
    /*
    Dune::FieldVector<typename GV::Grid::ctype,GV::dimension>
      x = is.geometry().global(x_);
    if (x[0] > 1.0-1e-6 && x[1] < 0.5+1e-6)
      return 2.5;
    else
    return 0.0;*/
  }

  //! set time for subsequent evaluation
  void setTime (typename Traits::RangeFieldType t)
  {
  }

private:
  typename Traits::PermTensorType kabslow, kabs;
  typename Traits::RangeFieldType gamma;
  typename Traits::RangeFieldType porosity_;
  const std::vector<int>& elementIndexToPhysicalGroup;
  const GV gv;
  const typename Traits::RangeFieldType bottomPotential, bottomflux;
};

// boundary condition type
SIMPLE_BOUNDARYTYPE_FUNCTION(DarcyBoundaryFunctionType,ig,x,y)
{
  Dune::FieldVector<typename GV::Grid::ctype,GV::dimension>
    xg = ig.geometry().global(x);

  if (!ig.boundary())
    {
      y = Traits::None; // no bc on subdomain interface
      return;
    }

  if (xg[0]<1E-6 || xg[0]>100-1E-6)
    {
      y = Traits::Neumann; // Neumann
      return;
    }
  y = Traits::Dirichlet; // Dirichlet
}
END_SIMPLE_BOUNDARYTYPE_FUNCTION


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
    : pgfs(&gfs), parameters(params), xg(x_), lfs(gfs), xl(gfs.maxLocalSize()), yb(gfs.maxLocalSize()), factor(factor_)
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
    return pgfs->gridview();
  }

private:
  CountingPointer<GFS const> pgfs;
  const Params& parameters;
  const X& xg;
  mutable typename GFS::LocalFunctionSpace lfs;
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
    : pgfs(&gfs), parameters(params), xg(x_)
  { }

  // Evaluate
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    // get and bind local functions space
    typename GFS::LocalFunctionSpace lfs(*pgfs);
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
    return pgfs->gridview();
  }

private:
  CountingPointer<GFS const> pgfs;
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

    typedef Dune::UGGrid<dim> BaseGrid;
    BaseGrid baseGrid(500);
    std::vector<int> boundaryIndexToPhysicalGroup, elementIndexToPhysicalGroup;
    Dune::GmshReader<BaseGrid> gmshreader;
    gmshreader.read(baseGrid,parameters["mesh.filename"],boundaryIndexToPhysicalGroup,elementIndexToPhysicalGroup,true,false);

    typedef BaseGrid::LeafGridView GV;

    GV gv = baseGrid.leafView();

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
    grid.startSubDomainMarking();
    for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
      {
        grid.addToSubDomain(elementIndexToPhysicalGroup[mdgv.indexSet().index(*it)] > 0 ? 1 : 0,*it);
      }
    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();

    typedef MDGV::Grid::ctype DF;
    typedef double RF;

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
    typedef Dune::PDELab::ISTLVectorBackend<1> VBE;

    NOCON con;
    DCON dcon;

    typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> EC;
    typedef Dune::PDELab::MultiDomain::SubDomainSupersetCondition<Grid> SC;

    EC c0(0);
    EC c1(1);

    typedef Dune::PDELab::GridFunctionSpace<SDGV,V_FEM,NOCON,VBE> V_GFS;
    V_GFS vgfs(stokesGV,vfem);

    typedef Dune::PDELab::PowerGridFunctionSpace<V_GFS,dim,Dune::PDELab::GridFunctionSpaceLexicographicMapper> PGFS_V_GFS;

    PGFS_V_GFS powervgfs(vgfs);

    typedef Dune::PDELab::GridFunctionSpace<SDGV,P_FEM,NOCON,VBE> P_GFS;
    P_GFS pgfs(stokesGV,pfem);

    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,PGFS_V_GFS,P_GFS> StokesGFS;

    StokesGFS stokesgfs(powervgfs,pgfs);

    typedef Dune::PDELab::GridFunctionSpace<SDGV,DarcyFEM,NOCON,VBE> DarcyGFS;
    DarcyGFS darcygfs(darcyGV,darcyfem);

    typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,StokesGFS,DarcyGFS> MultiGFS;
    MultiGFS multigfs(grid,stokesgfs,darcygfs);

    typedef MultiGFS::ConstraintsContainer<RF>::Type C;
    C cg;

    typedef RightScalarFunction<MDGV,RF> ScalarVelocityInitialFunction;
    typedef Dune::PDELab::PowerGridFunction<ScalarVelocityInitialFunction,dim> VelocityInitialFunction;
    typedef ZeroScalarFunction<MDGV,RF> PressureInitialFunction;
    typedef Dune::PDELab::CompositeGridFunction<VelocityInitialFunction,PressureInitialFunction> StokesInitialFunction;

    ScalarVelocityInitialFunction scalarVelocityInitialFunctionX(mdgv,-1.0);
    ScalarVelocityInitialFunction scalarVelocityInitialFunctionY(mdgv,1.0);
    VelocityInitialFunction velocityInitialFunction(scalarVelocityInitialFunctionX,scalarVelocityInitialFunctionY);
    PressureInitialFunction pressureInitialFunction(mdgv);
    StokesInitialFunction stokesInitialFunction(velocityInitialFunction,pressureInitialFunction);

    typedef PressureDropBoundaryFunction<MDGV> StokesScalarVelocityBoundaryFunction;
    typedef Dune::PDELab::PowerGridFunction<StokesScalarVelocityBoundaryFunction,dim> StokesVelocityBoundaryFunction;
    typedef ScalarNeumannBoundaryFunction<MDGV> StokesPressureBoundaryFunction;
    typedef Dune::PDELab::CompositeGridFunction<StokesVelocityBoundaryFunction,StokesPressureBoundaryFunction> StokesBoundaryFunction;

    StokesScalarVelocityBoundaryFunction stokesScalarVelocityBoundaryFunction(mdgv);
    StokesVelocityBoundaryFunction stokesVelocityBoundaryFunction(stokesScalarVelocityBoundaryFunction);
    StokesPressureBoundaryFunction stokesPressureBoundaryFunction(mdgv);
    StokesBoundaryFunction stokesBoundaryFunction(stokesVelocityBoundaryFunction,stokesPressureBoundaryFunction);

    typedef DarcyBoundaryFunctionType<MDGV> DarcyBoundaryFunction;
    DarcyBoundaryFunction darcyBoundaryFunction(mdgv);

    typedef PressureDropFlux<MDGV,RF> NeumannFlux;
    NeumannFlux neumannFlux(mdgv,parameters.get<RF>("boundaries.inflow"),parameters.get<RF>("boundaries.outflow"));

    typedef NavierStokesParameters<MDGV,RF> NavierStokesParams;
    NavierStokesParams navierStokesParams(parameters.sub("parameters"));

    typedef Dune::PDELab::TaylorHoodNavierStokes<
      StokesBoundaryFunction,NeumannFlux,NavierStokesParams,false,stokes_q>
      StokesOperator;
    StokesOperator stokesOperator(stokesBoundaryFunction,neumannFlux,navierStokesParams);


    typedef DarcyParameters<MDGV,RF> DarcyParams;
    DarcyParams darcyParams(mdgv,elementIndexToPhysicalGroup,parameters);

    typedef Dune::PM::ADRWIP<DarcyParams> DarcyOperator;
    DarcyOperator darcyOperator(darcyParams,Dune::PM::ADRWIPMethod::OBB,Dune::PM::ADRWIPWeights::weightsOn,0.0);

    typedef CouplingParameters<MDGV> CouplingParams;
    CouplingParams couplingParams(parameters.sub("parameters"));

    typedef StokesDarcyCouplingOperator<CouplingParams> CouplingOperator;
    CouplingOperator couplingOperator(couplingParams);

    typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,DCON,MultiGFS,DCON,StokesOperator,EC,StokesGFS> StokesSubProblem;

    typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,DarcyOperator,EC,DarcyGFS> DarcySubProblem;

    StokesSubProblem stokesSubProblem(dcon,dcon,stokesOperator,c0);
    DarcySubProblem darcySubProblem(con,con,darcyOperator,c1);

    //typedef Dune::PDELab::MultiDomain::EnrichedCoupling<StokesSubProblem,DarcySubProblem,CouplingOperator,2> Coupling;
    typedef Dune::PDELab::MultiDomain::Coupling<StokesSubProblem,DarcySubProblem,CouplingOperator> Coupling;
    Coupling coupling(stokesSubProblem,darcySubProblem,couplingOperator);

    Dune::PDELab::MultiDomain::trialSpaceConstraints(stokesBoundaryFunction,multigfs,cg,stokesBoundaryFunction,stokesSubProblem);
    std::cout << multigfs.size() << " DOF, " << cg.size() << " restricted" << std::endl;

    typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

    typedef Dune::PDELab::MultiDomain::MultiDomainGridOperatorSpace<MultiGFS,MultiGFS,MBE,StokesSubProblem,DarcySubProblem,Coupling> MultiGOS;

    MultiGOS multigos(multigfs,multigfs,cg,cg,stokesSubProblem,darcySubProblem,coupling);

    typedef MultiGFS::VectorContainer<RF>::Type V;

    V u(multigfs,0.0);

    //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
    //LS ls(5000,true);
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
    LS ls(true);

    typedef Dune::PDELab::StationaryLinearProblemSolver<MultiGOS,LS,V> PDESOLVER;
    PDESOLVER pdesolver(multigos,ls,1e-10);

    /*typedef Dune::PDELab::Newton<MultiGOS,LS,V> PDESOLVER;
    PDESOLVER pdesolver(multigos,ls);
    pdesolver.setReassembleThreshold(0.0);
    pdesolver.setVerbosityLevel(2);
    pdesolver.setReduction(1e-10);
    pdesolver.setMinLinearReduction(1e-4);
    pdesolver.setMaxIterations(25);
    pdesolver.setLineSearchMaxIterations(10);*/

    pdesolver.apply(u);

    V r(u);
    r = 0.0;
    multigos.residual(u,r);

    typedef Dune::PDELab::MultiDomain::TypeBasedGridFunctionSubSpace<MultiGFS,StokesGFS> StokesSubGFS;
    typedef Dune::PDELab::MultiDomain::TypeBasedGridFunctionSubSpace<MultiGFS,DarcyGFS> DarcySubGFS;
    typedef Dune::PDELab::GridFunctionSubSpace<StokesSubGFS,0> StokesVelocitySubGFS;
    typedef Dune::PDELab::GridFunctionSubSpace<StokesSubGFS,1> StokesPressureSubGFS;

    StokesSubGFS stokessubgfs(multigfs);
    DarcySubGFS darcysubgfs(multigfs);
    StokesVelocitySubGFS stokesvelocitysubgfs(stokessubgfs);
    StokesPressureSubGFS stokespressuresubgfs(stokessubgfs);

    /*
     * Output initial guess
     */


    //Dune::PDELab::FilenameHelper fn0("instationaryexplicitlycoupledpoisson-right");
    {
      typedef Dune::PDELab::VectorDiscreteGridFunction<StokesVelocitySubGFS,V> VDGF;
      typedef Dune::PDELab::DiscreteGridFunction<StokesPressureSubGFS,V> PDGF;
      VDGF vdgf(stokesvelocitysubgfs,u);
      PDGF pdgf(stokespressuresubgfs,u);
      VDGF vdgf2(stokesvelocitysubgfs,r);
      PDGF pdgf2(stokespressuresubgfs,r);
      Dune::SubsamplingVTKWriter<SDGV> vtkwriter(stokesGV,parameters.get("mesh.refineoutput",0));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf,"u"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf,"p"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf2,"ru"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf2,"rp"));
      vtkwriter.write("stokes",Dune::VTKOptions::binaryappended);
      //fn0.increment();
    }

    //Dune::PDELab::FilenameHelper fn1("instationaryexplicitlycoupledpoisson-left");
    {
      typedef Dune::PDELab::DarcyFlowFromPotential<DarcySubGFS,DarcyParams,V> VDGF;
      typedef Dune::PDELab::DiscreteGridFunction<DarcySubGFS,V> PhiDGF;
      typedef Dune::PDELab::DiscretePressureGridFunction<DarcySubGFS,DarcyParams,V> PDGF;
      VDGF vdgf(darcysubgfs,darcyParams,u);
      PhiDGF phidgf(darcysubgfs,u);
      PDGF pdgf(darcysubgfs,darcyParams,u,couplingParams.gravity());
      VDGF vdgf2(darcysubgfs,darcyParams,r);
      PhiDGF phidgf2(darcysubgfs,r);
      PDGF pdgf2(darcysubgfs,darcyParams,r,couplingParams.gravity());
      Dune::SubsamplingVTKWriter<SDGV> vtkwriter(darcyGV,parameters.get("mesh.refineoutput",0));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf,"u"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PhiDGF>(phidgf,"phi"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf,"p"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf2,"ru"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PhiDGF>(phidgf2,"rphi"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf2,"rp"));
      vtkwriter.write("darcy",Dune::VTKOptions::binaryappended);
      //fn1.increment();
    }

  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
	return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
	return 1;
  }

}
