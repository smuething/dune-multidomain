#include "config.h"

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

#include <dune/pdelab/finiteelementmap/opbfem.hh>
#include "adrwip.hh"

#include "stokesdarcycouplingoperator.hh"

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
      Dune::FieldVector<typename GV::Grid::ctype,GV::dimension>
      xg = ig.geometry().global(x);
      if (xg[0] < 1e-6)// || (xg[0] > 1.0-1e-6 && xg[1] < 1.5-1e-6))
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

  const DFT pressure;
  const DFT length;
  const DFT origin;
  const int direction;

public:
  PressureDropFlux (const GV& gv, const RF p_, const RF l_, const RF o_, const int d_)
    : BaseT(gv), pressure(p_), length(l_), origin(o_), direction(d_)
  {
    const int dim = GV::dimension;
    assert(direction >=0 && direction <dim);
  }

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    if (x[0] < 1e-6)
      y = x[1]*x[1];
    else
      y = -0.3;
  }
};



class NavierStokesParameters
  : public Dune::PDELab::TaylorHoodNavierStokesParameters<double>
{
public:

  double rho() const{
    return 1.0;
  }

  double mu() const{
    return 1.0;
  }
};


template<typename GV, typename RF>
class DarcyParameters
{
  typedef Dune::PM::ADRBoundaryConditions BC;

public:
  typedef Dune::PM::ADRTraits<GV,RF> Traits;

  //! constructor
  DarcyParameters()
  {
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        kabs[i][j] = (i==j) ? 1.0 : 0;
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
    return kabs;
  }

  //! nonlinear scalar diffusion term
  typename Traits::RangeFieldType
  w (const typename Traits::ElementType& e, const typename Traits::DomainType& x, typename Traits::RangeFieldType u) const
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
      if (x[0] > 1.0-1e-6 && x[1] < 0.5+1e-6)
        return BC::Flux;
      else
        return BC::Flux;
    }
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x_) const
  {
    typename Traits::DomainType x = is.geometry().global(x_);
    if (x[0]<1e-6)
      return 1.0;
    else
      return 0.0;
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x_) const
  {
    Dune::FieldVector<typename GV::Grid::ctype,GV::dimension>
      x = is.geometry().global(x_);
    if (x[0] > 1.0-1e-6 && x[1] < 0.5+1e-6)
      return 2.5;
    else
      return 0.0;
  }

  //! set time for subsequent evaluation
  void setTime (typename Traits::RangeFieldType t)
  {
  }

private:
  typename Traits::PermTensorType kabs;
};

int main(int argc, char** argv) {

  try {

    Dune::MPIHelper::instance(argc,argv);

    const int dim = 2;

#ifdef YASPGRID
    typedef Dune::YaspGrid<dim> BaseGrid;
    const Dune::FieldVector<int,dim> s(1);
    const Dune::FieldVector<double,dim> h(1.0);
    const Dune::FieldVector<bool,dim> p(false);
    BaseGrid baseGrid(h,s,p,0);
    baseGrid.globalRefine(2);
#endif

#ifdef UGGRID
    typedef Dune::UGGrid<dim> BaseGrid;
    BaseGrid baseGrid(500);
    std::vector<int> boundaryIndexToPhysicalGroup, elementIndexToPhysicalGroup;
    Dune::GmshReader<BaseGrid> gmshreader;
    gmshreader.read(baseGrid,"stokesdarcy.msh",boundaryIndexToPhysicalGroup,elementIndexToPhysicalGroup,true,false);
#endif

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
#ifdef UGGRID
        grid.addToSubDomain(elementIndexToPhysicalGroup[mdgv.indexSet().index(*it)],*it);
#else
        // FIXME: !!!!!!!!
#endif
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

    typedef Dune::PDELab::MultiDomain::SubProblemSubProblemInterface<MDGV,EC,EC> StokesDarcyBoundary;

    StokesDarcyBoundary stokesDarcyBoundary(mdgv,c0,c1);

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

    typedef Dune::PDELab::MultiDomain::CouplingGridFunctionSpace<MDGV,CouplingFEM,StokesDarcyBoundary,NOCON,VBE> CouplingGFS;

    CouplingGFS couplinggfs(mdgv,couplingfem,stokesDarcyBoundary);

    //typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,StokesGFS,DarcyGFS,CouplingGFS> MultiGFS;
    //MultiGFS multigfs(grid,stokesgfs,darcygfs,couplinggfs);
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

    typedef PressureDropFlux<MDGV,RF> NeumannFlux;
    NeumannFlux neumannFlux(mdgv,5.0,0,0,0);

    NavierStokesParameters navierStokesParameters;

    typedef Dune::PDELab::TaylorHoodNavierStokes<
      StokesBoundaryFunction,NeumannFlux,NavierStokesParameters,true,stokes_q>
      StokesOperator;
    StokesOperator stokesOperator(stokesBoundaryFunction,neumannFlux,navierStokesParameters);


    typedef DarcyParameters<MDGV,RF> DarcyParams;
    DarcyParams darcyParams;

    typedef Dune::PM::ADRWIP<DarcyParams> DarcyOperator;
    DarcyOperator darcyOperator(darcyParams,Dune::PM::ADRWIPMethod::OBB,Dune::PM::ADRWIPWeights::weightsOn,0.0);

    typedef CouplingParameters<MDGV> CouplingParams;
    CouplingParams couplingParams;

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
    V r(u);

    //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
    //LS ls(5000,true);
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
    LS ls(true);

    typedef Dune::PDELab::StationaryLinearProblemSolver<MultiGOS,LS,V> PDESOLVER;
    PDESOLVER pdesolver(multigos,ls,1e-10);

    pdesolver.apply(u);

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
      Dune::VTKWriter<SDGV> vtkwriter(stokesGV,Dune::VTKOptions::conforming);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf,"u"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf,"p"));
      vtkwriter.write("stokes",Dune::VTKOptions::binaryappended);
      //fn0.increment();
    }

    //Dune::PDELab::FilenameHelper fn1("instationaryexplicitlycoupledpoisson-left");
    {
      typedef Dune::PDELab::DiscreteGridFunctionGradient<DarcySubGFS,V> VDGF;
      typedef Dune::PDELab::DiscreteGridFunction<DarcySubGFS,V> PDGF;
      VDGF vdgf(darcysubgfs,u);
      PDGF pdgf(darcysubgfs,u);
      Dune::VTKWriter<SDGV> vtkwriter(darcyGV,Dune::VTKOptions::conforming);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf,"v"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf,"phi"));
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
