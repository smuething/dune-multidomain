#include "config.h"

#include <dune/common/parametertreeparser.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinggridfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinglocalfunctionspace.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include <dune/pdelab/localoperator/cg_stokes.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/multidomain/vtk.hh>

#include<typeinfo>

#include <dune/pdelab/finiteelementmap/opbfem.hh>
// #include "../../dune-pm/dune/pm/models/adrwip.hh"
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>

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


struct ScalarNeumannBoundaryType :
  public Dune::PDELab::DirichletConstraintsParameters
{

  template<typename I>
  bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    return false;
  }

  template<typename I>
  bool isNeumann(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    return true;
  }

  template<typename R>
  void setTime(const R& r)
  {
  }

};

template<typename GV>
struct StokesBoundaryType
{

  typedef Dune::PDELab::BoundaryGridFunctionTraits<
    GV,
    Dune::PDELab::StokesBoundaryCondition::Type,
    1,
    Dune::FieldVector<
      Dune::PDELab::StokesBoundaryCondition::Type,
      1
      >
    > Traits;

  void setTime(double t) const
  {}

  template<typename I>
  void evaluate(const I& ig,
                const typename Traits::DomainType& x,
                typename Traits::RangeType& y) const
  {
    if (!ig.boundary())
      {
        y = Dune::PDELab::StokesBoundaryCondition::DoNothing;
        return;
      }
    auto xg = ig.geometry().global(x);
    if (xg[0] < 1e-6 || xg[0] > 100-1e-6)
      y = Dune::PDELab::StokesBoundaryCondition::StressNeumann;
    else
      y = Dune::PDELab::StokesBoundaryCondition::VelocityDirichlet;
    y = Dune::PDELab::StokesBoundaryCondition::StressNeumann;
  }

};

struct PressureDropBoundaryType :
  public Dune::PDELab::DirichletConstraintsParameters
{

  template<typename I>
  bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    if (!ig.boundary())
      return false; // no bc
    auto xg = ig.geometry().global(x);
    if (xg[0] < 1e-6 || xg[0] > 100-1e-6)
      return true; // Neumann
    else
      return true; // Dirichlet
  }

  template<typename I>
  bool isNeumann(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    if (!ig.boundary())
      return false; // no bc
    auto xg = ig.geometry().global(x);
    if (xg[0] < 1e-6 || xg[0] > 100-1e-6)
      return false; // Neumann
    else
      return false; // Dirichlet
  }

  template<typename R>
  void setTime(const R& r)
  {
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
      y = x[1] > 50 ? inflow : -outflow;//0.04*(100-x[1])*(x[1]-110); // parabolic flow profile
    else if (x[0] > 100-1e-6)
      y = x[1] > 50 ? -outflow : inflow;//*0.04*(100-x[1])*(x[1]-110); // parabolic flow profile
    else
      y = 0;//-10*((x[1]-0.5)*(1.5-x[1])+0.02);
    y *= 100.0 - x[1];
    std::cout << x << std::endl;
  }

};


template<typename GV, typename RF, typename B, typename V, typename J>
class NavierStokesParameters
  : public Dune::PDELab::NavierStokesDefaultParameters<GV,RF,V,B,V,J,false,true>
{
public:

  typedef Dune::PDELab::NavierStokesDefaultParameters<GV,RF,V,B,V,J,false,true> BaseT;
  typedef typename BaseT::Traits Traits;

  NavierStokesParameters(const Dune::ParameterTree& params, B& b, V& v, J& j)
    : BaseT(params.sub("fluid"),v,b,v,j)
    , _gravity(0.0)
  {
    _gravity[GV::Grid::dimensionworld-1] = - params.get("gravity",9.81);
  }

  template<typename EG>
  typename Traits::VelocityRange
  f(const EG& e, const typename Traits::Domain& x) const
  {
    return _gravity;
  }


private:
  Dune::FieldVector<RF,GV::Grid::dimensionworld> _gravity;
};


template<typename GV, typename RF>
class DarcyParameters
{

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions BC;


  //! constructor
  DarcyParameters(GV gridview,
                  const std::vector<int>& physicalGroupMap,
                  const Dune::ParameterTree& params)
    : gamma(params.get("parameters.coupling.gamma",1.0))
    , porosity_(params.get<RF>("parameters.soil.porosity"))
    , elementIndexToPhysicalGroup(physicalGroupMap)
    , gv(gridview)
    , rightTopPotential(params.get<RF>("boundaries.rightpotential.top"))
    , rightBottomPotential(params.get<RF>("boundaries.rightpotential.bottom"))
    , bottomflux(params.get<RF>("boundaries.bottomflux"))
    , _gravity(params.get<RF>("parameters.gravity"))
    , _density(params.get<RF>("parameters.fluid.rho"))
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
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x_) const
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

  typename Traits::RangeFieldType
  gravity (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return _gravity;
  }

  typename Traits::RangeFieldType
  density (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return _density;
  }

  //! reaction term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x_) const
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
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x_) const
  {
    if (!is.boundary())
      return BC::None;
    else {
      Dune::FieldVector<typename GV::Grid::ctype,GV::dimension>
        x = is.geometry().global(x_);
      if (x[1] < 1e-6 && x[0] < 100)
          return BC::Outflow;
      else if (x[0] > 100 - 1e-6 && x[1] < 14.5)
        return x[1] < 13 ? BC::Dirichlet : BC::Neumann;
      else
        return BC::Neumann;
    }
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& is, const typename Traits::DomainType& x_) const
  {
    auto x = is.geometry().global(x_);
    auto scale = x[1] / 15.0;
    return scale * rightTopPotential + (1.0 - scale) * rightBottomPotential;
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x_) const
  {
    auto x = is.geometry().global(x_);
    if (x[1] < 1e-6)
      return bottomflux * x[0] * (100.0 - x[0]);
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

  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x_) const
  {
    return j(is,x_);
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
  const typename Traits::RangeFieldType rightTopPotential;
  const typename Traits::RangeFieldType rightBottomPotential;
  const typename Traits::RangeFieldType bottomflux;
  const typename Traits::RangeFieldType _gravity;
  const typename Traits::RangeFieldType _density;
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


template<typename LFS, typename Data, typename DarcyParameters>
class VTKPressureFromPotential
  : public Dune::PDELab::VTKScalarDiscreteFunctionBase<LFS,
                                                       Data,
                                                       VTKPressureFromPotential<
                                                         LFS,
                                                         Data,
                                                         DarcyParameters
                                                         >
                                                       >
{

  typedef Dune::PDELab::VTKScalarDiscreteFunctionBase<
    LFS,
    Data,
    VTKPressureFromPotential
    > Base;

public:

  using Base::Computations::evaluateSolution;

  typedef typename Base::Traits Traits;

  VTKPressureFromPotential(const LFS& lfs,
                           const std::shared_ptr<Data>& data,
                           const DarcyParameters& params)
    : Base(lfs,data)
    , _params(params)
  {}

  // Evaluate
  void evaluate(const typename Traits::ElementType& e,
                const typename Traits::DomainType& x,
                typename Traits::RangeType& y) const
  {
    this->bind(e,x);

    y = this->solution();
    auto& multidomain_entity = this->gridView().grid().multiDomainEntity(e);
    typename Traits::ElementType::Geometry::GlobalCoordinate xg = e.geometry().global(x);
    auto gravity = _params.gravity(multidomain_entity,x);
    auto density = _params.density(multidomain_entity,x);
    y -= xg[Traits::dimDomain-1];
    y *= gravity * density;
  }

private:

  const DarcyParameters& _params;

};


template<typename Params>
class VTKDarcyFlowFromPotential
{

public:

  template<typename LFS, typename Data>
  class Function
    : public Dune::PDELab::VTKCoordinateDiscreteFunctionBase<LFS,Data,Function<LFS,Data> >
  {

    typedef Dune::PDELab::VTKCoordinateDiscreteFunctionBase<
      LFS,
      Data,
      Function
      > Base;

  public:

    using Base::Computations::evaluateGradient;

    typedef typename Base::Traits Traits;

    Function(const LFS& lfs, const std::shared_ptr<Data>& data, const Params& params)
      : Base(lfs,data)
      , _params(params)
    {}

    // Evaluate
    void evaluate(const typename Traits::ElementType& e,
                  const typename Traits::DomainType& x,
                  typename Traits::RangeType& y) const
    {
      this->bind(e,x);

      auto gradient = this->gradient()[0];
      auto& multidomain_entity = this->gridView().grid().multiDomainEntity(e);

      _params.A(multidomain_entity,x).mv(gradient,y);
      y /= -_params.porosity(multidomain_entity,x);

    }

  private:

    const Params& _params;

  };


  template<typename LFS, typename Data>
  std::shared_ptr<Function<LFS,Data> > create(const LFS& lfs, std::shared_ptr<Data> data)
  {
    return std::make_shared<Function<LFS,Data> >(lfs,data,_params);
  }

  VTKDarcyFlowFromPotential(const Params& params)
    : _params(params)
  {}

private:

  const Params& _params;

};



int main(int argc, char** argv) {

  try {

    Dune::MPIHelper::instance(argc,argv);

    Dune::ParameterTree parameters;
    Dune::ParameterTreeParser::readINITree(argv[1],parameters);

    const int dim = 2;

    typedef Dune::UGGrid<dim> BaseGrid;
    std::vector<int> boundaryIndexToPhysicalGroup, elementIndexToPhysicalGroup;
    Dune::GmshReader<BaseGrid> gmshreader;
    Dune::shared_ptr<BaseGrid> baseGridPtr(gmshreader.read(parameters["mesh.filename"],boundaryIndexToPhysicalGroup,elementIndexToPhysicalGroup,true,false));

    BaseGrid& baseGrid = *baseGridPtr;

    typedef BaseGrid::LeafGridView GV;

    GV gv = baseGrid.leafGridView();

    typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,4> > Grid;
    Grid grid(baseGrid,false);
    typedef Grid::SubDomainGrid SubDomainGrid;
    SubDomainGrid& stokesGrid = grid.subDomain(0);
    SubDomainGrid& darcyGrid = grid.subDomain(1);
    typedef Grid::ctype ctype;
    typedef Grid::LeafGridView MDGV;
    typedef SubDomainGrid::LeafGridView SDGV;
    MDGV mdgv = grid.leafGridView();
    SDGV stokesGV = stokesGrid.leafGridView();
    SDGV darcyGV = darcyGrid.leafGridView();
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

    typedef Dune::PDELab::PkLocalFiniteElementMap<SDGV,DF,RF,stokes_k> V_FEM;
    typedef Dune::PDELab::PkLocalFiniteElementMap<SDGV,DF,RF,stokes_k-1> P_FEM;
    typedef Dune::PDELab::OPBLocalFiniteElementMap<DF,RF,darcy_k,dim,Dune::GeometryType::simplex,Dune::GMPField<512> > DarcyFEM;
    typedef Dune::PDELab::PkLocalFiniteElementMap<MDGV,DF,RF,dim-1> CouplingFEM;

    V_FEM vfem(stokesGV);
    P_FEM pfem(stokesGV);
    DarcyFEM darcyfem;
    CouplingFEM couplingfem(mdgv);

    typedef Dune::PDELab::NoConstraints NOCON;
    typedef Dune::PDELab::ConformingDirichletConstraints DCON;
    typedef Dune::PDELab::ISTLVectorBackend<> VBE;

    NOCON con;
    DCON dcon;

    typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> EC;
    typedef Dune::PDELab::MultiDomain::SubDomainSupersetCondition<Grid> SC;

    EC c0({0});
    EC c1({1});

    typedef Dune::PDELab::VectorGridFunctionSpace<
      SDGV,
      V_FEM,
      dim,
      VBE,
      VBE,
      DCON
      > PGFS_V_GFS;
    PGFS_V_GFS powervgfs(stokesGV,vfem);
    powervgfs.name("v");

    typedef Dune::PDELab::GridFunctionSpace<SDGV,P_FEM,DCON,VBE> P_GFS;
    P_GFS pgfs(stokesGV,pfem);
    pgfs.name("p");

    typedef Dune::PDELab::CompositeGridFunctionSpace<
      VBE,
      Dune::PDELab::LexicographicOrderingTag,
      PGFS_V_GFS,
      P_GFS> StokesGFS;

    StokesGFS stokesgfs(powervgfs,pgfs);

    typedef Dune::PDELab::GridFunctionSpace<SDGV,DarcyFEM,NOCON,VBE> DarcyGFS;
    DarcyGFS darcygfs(darcyGV,darcyfem);
    darcygfs.name("phi");

    typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<
      Grid,
      VBE,
      Dune::PDELab::LexicographicOrderingTag,
      StokesGFS,
      DarcyGFS
      > MultiGFS;
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

    typedef StokesBoundaryType<MDGV> StokesBoundary;
    StokesBoundary stokesBoundary;

    /*
    typedef PressureDropBoundaryType StokesScalarVelocityBoundaryFunction;
    typedef Dune::PDELab::PowerConstraintsParameters<StokesScalarVelocityBoundaryFunction,dim> StokesVelocityBoundaryFunction;
    typedef ScalarNeumannBoundaryType StokesPressureBoundaryFunction;
    typedef Dune::PDELab::CompositeConstraintsParameters<StokesVelocityBoundaryFunction,StokesPressureBoundaryFunction> StokesBoundaryFunction;

    StokesScalarVelocityBoundaryFunction stokesScalarVelocityBoundaryFunction;
    StokesVelocityBoundaryFunction stokesVelocityBoundaryFunction(stokesScalarVelocityBoundaryFunction);
    StokesPressureBoundaryFunction stokesPressureBoundaryFunction;
    StokesBoundaryFunction stokesBoundaryFunction(stokesVelocityBoundaryFunction,stokesPressureBoundaryFunction);

    */

    typedef PressureDropFlux<MDGV,RF> NeumannFlux;
    NeumannFlux neumannFlux(mdgv,parameters.get<RF>("boundaries.inflow"),parameters.get<RF>("boundaries.outflow"));

    typedef NavierStokesParameters<MDGV,RF,StokesBoundary,StokesInitialFunction,NeumannFlux> NavierStokesParams;
    NavierStokesParams navierStokesParams(parameters.sub("parameters"),stokesBoundary,stokesInitialFunction,neumannFlux);

    typedef Dune::PDELab::StokesVelocityDirichletConstraints<NavierStokesParams> StokesVelocityBoundaryFunction;
    typedef Dune::PDELab::StokesPressureDirichletConstraints<NavierStokesParams> StokesPressureBoundaryFunction;
    typedef Dune::PDELab::CompositeConstraintsParameters<
      StokesVelocityBoundaryFunction,
      StokesPressureBoundaryFunction
      > StokesBoundaryFunction;

    StokesVelocityBoundaryFunction stokesVelocityBoundaryFunction(navierStokesParams);
    StokesPressureBoundaryFunction stokesPressureBoundaryFunction(navierStokesParams);
    StokesBoundaryFunction stokesBoundaryFunction(stokesVelocityBoundaryFunction,stokesPressureBoundaryFunction);

    typedef Dune::PDELab::TaylorHoodNavierStokes<NavierStokesParams> StokesOperator;
    StokesOperator stokesOperator(navierStokesParams,stokes_q);

    typedef DarcyParameters<MDGV,RF> DarcyParams;
    DarcyParams darcyParams(mdgv,elementIndexToPhysicalGroup,parameters);

    typedef DarcyBoundaryTypeAdapter<DarcyParams> DarcyBoundaryFunction;
    DarcyBoundaryFunction darcyBoundaryFunction(darcyParams);

    typedef Dune::PDELab::ConvectionDiffusionDG<DarcyParams,DarcyFEM> DarcyOperator;
    DarcyOperator darcyOperator(
      darcyParams,
      Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,
      Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,
      parameters.get("parameters.darcy.alpha",1.0)
    );

    typedef CouplingParameters<MDGV> CouplingParams;
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


    auto constraints = Dune::PDELab::MultiDomain::constraints<RF>(
      multigfs,
      Dune::PDELab::MultiDomain::constrainSubProblem(
        stokesSubProblem,
        stokesBoundaryFunction
      )
    );

    constraints.assemble(cg);

    //Dune::PDELab::MultiDomain::trialSpaceConstraints(stokesBoundaryFunction,multigfs,cg,stokesBoundaryFunction,stokesSubProblem);
    std::cout << multigfs.size() << " DOF, " << cg.size() << " restricted" << std::endl;

    typedef Dune::PDELab::MultiDomain::GridOperator<
      MultiGFS,MultiGFS,
      Dune::PDELab::ISTLMatrixBackend,
      RF,RF,RF,C,C,
      StokesSubProblem,
      DarcySubProblem,
      Coupling
      > GridOperator;

    GridOperator gridoperator(multigfs,multigfs,cg,cg,stokesSubProblem,darcySubProblem,coupling);

    typedef GridOperator::Traits::Domain V;

    V u(multigfs,0.0);

    //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
    //LS ls(5000,true);
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
    LS ls(true);

    typedef Dune::PDELab::StationaryLinearProblemSolver<GridOperator,LS,V> PDESOLVER;
    PDESOLVER pdesolver(gridoperator,ls,1e-15);

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
    gridoperator.residual(u,r);

    /*
    typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,TypeTree::TreePath<0,0> > StokesVelocitySubGFS;
    typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,TypeTree::TreePath<0,1> > StokesPressureSubGFS;
    typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,TypeTree::TreePath<1> > DarcySubGFS;

    StokesSubGFS stokessubgfs(multigfs);
    DarcySubGFS darcysubgfs(multigfs);
    StokesVelocitySubGFS stokesvelocitysubgfs(stokessubgfs);
    StokesPressureSubGFS stokespressuresubgfs(stokessubgfs);
    */

    /*
     * Output initial guess
     */


    //Dune::PDELab::FilenameHelper fn0("instationaryexplicitlycoupledpoisson-right");
    {
      /*
      typedef Dune::PDELab::VectorDiscreteGridFunction<StokesVelocitySubGFS,V> VDGF;
      typedef Dune::PDELab::DiscreteGridFunction<StokesPressureSubGFS,V> PDGF;
      VDGF vdgf(stokesvelocitysubgfs,u);
      PDGF pdgf(stokespressuresubgfs,u);
      VDGF vdgf2(stokesvelocitysubgfs,r);
      PDGF pdgf2(stokespressuresubgfs,r);
      */
      Dune::SubsamplingVTKWriter<SDGV> vtkwriter(stokesGV,parameters.get("mesh.refineoutput",0));
      /*
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf,"u"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf,"p"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf2,"ru"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf2,"rp"));
      */
      Dune::PDELab::MultiDomain::addSolutionToVTKWriter(
        vtkwriter,
        multigfs,
        u,
        Dune::PDELab::MultiDomain::subdomain_predicate<Grid::SubDomainIndex>(stokesGV.grid().domain())
      );
      vtkwriter.write("stokes",Dune::VTK::appendedraw);
      //fn0.increment();
    }

    //Dune::PDELab::FilenameHelper fn1("instationaryexplicitlycoupledpoisson-left");
    {
      /*
      typedef Dune::PDELab::DarcyFlowFromPotential<DarcySubGFS,DarcyParams,V> VDGF;
      typedef Dune::PDELab::DiscreteGridFunction<DarcySubGFS,V> PhiDGF;
      typedef Dune::PDELab::DiscretePressureGridFunction<DarcySubGFS,DarcyParams,V> PDGF;
      VDGF vdgf(darcysubgfs,darcyParams,u);
      PhiDGF phidgf(darcysubgfs,u);
      PDGF pdgf(darcysubgfs,darcyParams,u,couplingParams.gravity());
      VDGF vdgf2(darcysubgfs,darcyParams,r);
      PhiDGF phidgf2(darcysubgfs,r);
      PDGF pdgf2(darcysubgfs,darcyParams,r,couplingParams.gravity());
      */
      Dune::SubsamplingVTKWriter<SDGV> vtkwriter(darcyGV,parameters.get("mesh.refineoutput",0));
      /*
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf,"u"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PhiDGF>(phidgf,"phi"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf,"p"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf2,"ru"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PhiDGF>(phidgf2,"rphi"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf2,"rp"));
      */

      Dune::PDELab::MultiDomain::addSolutionToVTKWriter(
        vtkwriter,
        multigfs,
        u,
        Dune::PDELab::MultiDomain::subdomain_predicate<Grid::SubDomainIndex>(darcyGV.grid().domain())
      ).addVertexFunction<VTKPressureFromPotential>(
        Dune::TypeTree::TreePath<1>(),
        "p",
        darcyParams
      ).addVertexFunction(
        VTKDarcyFlowFromPotential<DarcyParams>(darcyParams),
        Dune::TypeTree::TreePath<1>(),
        "v"
      );
      vtkwriter.write("darcy",Dune::VTK::appendedraw);
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
