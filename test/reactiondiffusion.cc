#include "config.h"

#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/instationarymultidomaingridoperatorspace.hh>
#include <dune/pdelab/multidomain/instationarysubproblem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include <dune/pdelab/localoperator/diffusion.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>

#include<typeinfo>

#include "functionmacros.hh"
#include "../../dune-pdelab-howto/src/course-examples/example05_operator.hh"
#include "../../dune-pdelab-howto/src/course-examples/example05_toperator.hh"
#include "../../dune-pdelab-howto/src/course-examples/example05_initial.hh"
#include "simpletimeoperator.hh"

template<typename GV, typename RF>
class K
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
                                                  K<GV,RF> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> > Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,K<GV,RF> > BaseT;

  K (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
							  typename Traits::RangeType& y) const
  {
    for (int i = 0; i < GV::dimension; ++i)
      for (int j = 0; j < GV::dimension; ++j)
        y[i][j] = i == j ? _permeability : 0.0;
  }

  K(const GV& gv, double permeability) :
    BaseT(gv),
    _permeability(permeability)
  {}

private:
  const double _permeability;
};

// source term
SIMPLE_ANALYTIC_FUNCTION(A0,x,y)
{
  y=0;
}
END_SIMPLE_ANALYTIC_FUNCTION

// source term
SIMPLE_ANALYTIC_FUNCTION(F,x,y)
{
  if (x[0]>0.75 && x[0]<0.875 && x[1]>0.25 && x[1]<0.375)
    y = 50.0;
  else
    y = 0.0;
  //y=0;
}
END_SIMPLE_ANALYTIC_FUNCTION


// boundary condition type
SIMPLE_BOUNDARYTYPE_FUNCTION(B,ig,x,y)
{
  Dune::FieldVector<typename GV::Grid::ctype,GV::dimension>
    xg = ig.geometry().global(x);

  if (!ig.boundary())
    {
      y = Traits::None; // no bc on subdomain interface
      return;
    }

  y = Traits::Neumann;
/*
  if (xg[1]<1E-6 || xg[1]>1.0-1E-6)
    {
      y = Traits::Neumann; // Neumann
      return;
    }
  if (xg[0]>1.0-1E-6 && xg[1]>0.5+1E-6)
    {
      y = Traits::Neumann; // Neumann
      return;
    }
  y = Traits::Dirichlet; // Dirichlet*/
}
END_SIMPLE_BOUNDARYTYPE_FUNCTION


// dirichlet bc
SIMPLE_ANALYTIC_FUNCTION(G,x,y)
{
  typename Traits::DomainType center;
  for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
  center -= x;
  y = exp(-center.two_norm2());
}
END_SIMPLE_ANALYTIC_FUNCTION


// neumann bc
SIMPLE_ANALYTIC_FUNCTION(J,x,y)
{
  if (x[1]<1E-6 || x[1]>1.0-1E-6)
    {
      y = 0;
      return;
    }
  if (x[0]>1.0-1E-6 && x[1]>0.5+1E-6)
    {
      y = -5.0;
      return;
    }
}
END_SIMPLE_ANALYTIC_FUNCTION



int main(int argc, char** argv) {

  if (argc < 5) {
    std::cerr << "Usage: " << argv[0] << " <refinement level> <dtstart> <dtmax> <tend>" << std::endl;
    return 1;
  }

  Dune::MPIHelper::instance(argc,argv);

  const double dtstart = atof(argv[2]);
  const double dtmax = atof(argv[3]);
  const double tend = atof(argv[4]);

  const int dim = 2;
  typedef Dune::YaspGrid<dim> BaseGrid;
  const Dune::FieldVector<double,dim> h(2.0);
  const Dune::array<int,dim> s = { {1,1} };
  const std::bitset<dim> p(false);
  BaseGrid baseGrid(h,s,p,0);
  baseGrid.globalRefine(atoi(argv[1]));
  typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,4> > Grid;
  Grid grid(baseGrid,false);
  typedef Grid::SubDomainGrid SubDomainGrid;
  SubDomainGrid& sdg0 = grid.subDomain(0);
  SubDomainGrid& sdg1 = grid.subDomain(1);
  typedef Grid::ctype ctype;
  typedef Grid::LeafGridView MDGV;
  typedef SubDomainGrid::LeafGridView SDGV;
  MDGV mdgv = grid.leafGridView();
  SDGV sdgv0 = sdg0.leafGridView();
  SDGV sdgv1 = sdg1.leafGridView();
  sdg0.hostEntityPointer(*sdgv0.begin<0>());
  grid.startSubDomainMarking();
  for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
    {
      if (it->geometry().center()[0] > 1.0+1e-8)
        grid.addToSubDomain(0,*it);
      else
        grid.addToSubDomain(1,*it);
    }
  grid.preUpdateSubDomains();
  grid.updateSubDomains();
  grid.postUpdateSubDomains();

  typedef MDGV::Grid::ctype DF;

  typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM;

  typedef FEM::Traits::LocalFiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  FEM fem0, fem1;
  typedef Dune::PDELab::NoConstraints NOCON;
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;

  typedef Dune::PDELab::GridFunctionSpace<MDGV,FEM,NOCON,VBE> GFS0;
  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM,NOCON,VBE> GFS1;

  typedef GFS0::ConstraintsContainer<R>::Type C;
  C cg;

  GFS0 gfs0(mdgv,fem0);
  GFS1 gfs1(sdgv1,fem1);

  typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,GFS0,GFS1> MultiGFS;

  MultiGFS multigfs(grid,gfs0,gfs1);

  double d_0 = 0.00028, d_1 = 0.005, lambda = 1.0, sigma = 1.0, kappa = -0.05, tau = 0.1;

  typedef B<MDGV> BType;
  BType b(mdgv);

  typedef F<MDGV,R> FType;
  FType f(mdgv);

  typedef G<MDGV,R> GType;
  GType g(mdgv);

  typedef J<MDGV,R> JType;
  JType j(mdgv);

  typedef K<MDGV,R> KType;
  KType k(mdgv,d_0);

  typedef A0<MDGV,R> A0Type;
  A0Type a0(mdgv);

  typedef Dune::PDELab::Diffusion<KType,A0Type,FType,BType,JType> LOP0;
  LOP0 lop0(k,a0,f,b,j,2);

  typedef SimpleTimeOperator TLOP0;
  TLOP0 tlop0(2);

  typedef Example05LocalOperator LOP1;
  LOP1 lop1(d_0,d_1,lambda,sigma,kappa,2);

  typedef Example05TimeLocalOperator TLOP1;
  TLOP1 tlop1(tau,2);

  typedef MDGV::IndexSet::SubDomainSet SDS;
  typedef Dune::PDELab::MultiDomain::EqualsSubDomains<SDS> EC;

  EC ec0(0);
  EC ec1(1);

  NOCON nocon;

  typedef Dune::PDELab::MultiDomain::InstationarySubProblem<double,MultiGFS,NOCON,MultiGFS,NOCON,LOP0,TLOP0,EC,0> SubProblem0;
  typedef Dune::PDELab::MultiDomain::InstationarySubProblem<double,MultiGFS,NOCON,MultiGFS,NOCON,LOP1,TLOP1,EC,0,1> SubProblem1;
  SubProblem0 sp0(nocon,nocon,lop0,tlop0,ec0);
  SubProblem1 sp1(nocon,nocon,lop1,tlop1,ec1);

  typedef U0Initial<MDGV,double> U0InitialType;
  U0InitialType u0initial(mdgv);
  typedef U1Initial<MDGV,double> U1InitialType;
  U1InitialType u1initial(mdgv);
  typedef Dune::PDELab::CompositeGridFunction<U0InitialType,U1InitialType> UInitialType;
  UInitialType uinitial(u0initial,u1initial);

  SubProblem0::Traits::LocalTrialFunctionSpace
    splfs0(sp0,sp0.trialGridFunctionSpaceConstraints());
  SubProblem1::Traits::LocalTrialFunctionSpace
    splfs1(sp1,sp1.trialGridFunctionSpaceConstraints());

  typedef MultiGFS::VectorContainer<R>::Type V;
  V uold(multigfs);
  uold = 0.0;

  Dune::PDELab::MultiDomain::interpolate(multigfs,uold,f,splfs0,uinitial,splfs1);

  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

  typedef Dune::PDELab::MultiDomain::InstationaryMultiDomainGridOperatorSpace<double,V,MultiGFS,MultiGFS,MBE,SubProblem0,SubProblem1> MultiGOS;

  MultiGOS multigos(multigfs,multigfs,cg,cg,sp0,sp1);

  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,false);

  typedef Dune::PDELab::Newton<MultiGOS,LS,V> PDESOLVER;
  PDESOLVER pdesolver(multigos,ls);
  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(2);
  pdesolver.setReduction(1e-10);
  pdesolver.setMinLinearReduction(1e-4);
  pdesolver.setMaxIterations(25);
  pdesolver.setLineSearchMaxIterations(10);

  Dune::PDELab::Alexander2Parameter<double> method;
  Dune::PDELab::OneStepMethod<double,MultiGOS,PDESOLVER,V,V> osm(method,multigos,pdesolver);
  osm.setVerbosityLevel(2);

  typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,1> SGFS0;
  typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,0> SGFS1;
  SGFS0 sgfs0(multigfs);
  SGFS1 sgfs1(multigfs);

  Dune::PDELab::FilenameHelper fn0("reactiondiffusion-c0");
  {
    typedef Dune::PDELab::DiscreteGridFunction<SGFS0,V> U0DGF;
    U0DGF u0dgf(sgfs0,uold);
    Dune::VTKWriter<MDGV> vtkwriter(mdgv,Dune::VTK::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"c0"));
    vtkwriter.write(fn0.getName(),Dune::VTK::appendedraw);
    fn0.increment();
  }
  Dune::PDELab::FilenameHelper fn1("reactiondiffusion-c1");
  {
    typedef Dune::PDELab::DiscreteGridFunction<SGFS1,V> U1DGF;
    U1DGF u1dgf(sgfs1,uold);
    Dune::VTKWriter<SDGV> vtkwriter(sdgv1,Dune::VTK::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"c1"));
    vtkwriter.write(fn1.getName(),Dune::VTK::appendedraw);
    fn1.increment();
  }


  V unew(multigfs,0.0);
  unew = uold;
  double dt = dtstart;
  double time = 0;

  while(time<tend-1e-8)
    {
      osm.apply(time,dt,uold,unew);
      uold = unew;
      time += dt;
      if (dt<dtmax-1e-8) dt = std::min(dt*1.1,dtmax);
      {
        typedef Dune::PDELab::DiscreteGridFunction<SGFS0,V> U0DGF;
        U0DGF u0dgf(sgfs0,uold);
        Dune::VTKWriter<MDGV> vtkwriter(mdgv,Dune::VTK::conforming);
        vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"c0"));
        vtkwriter.write(fn0.getName(),Dune::VTK::appendedraw);
        fn0.increment();
      }
      {
        typedef Dune::PDELab::DiscreteGridFunction<SGFS1,V> U1DGF;
        U1DGF u1dgf(sgfs1,uold);
        Dune::VTKWriter<SDGV> vtkwriter(sdgv1,Dune::VTK::conforming);
        vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"c1"));
        vtkwriter.write(fn1.getName(),Dune::VTK::appendedraw);
        fn1.increment();
      }
    }
}
