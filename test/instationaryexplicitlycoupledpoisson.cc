#include "config.h"

#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/instationarymultidomaingridoperatorspace.hh>
#include <dune/pdelab/multidomain/instationarysubproblem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/multidomain/instationarycoupling.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include<typeinfo>

#include "functionmacros.hh"
#include "proportionalflowcoupling.hh"
#include "simpletimeoperator.hh"
#include "instationarypoissonoperator.hh"

// source term
INSTATIONARY_ANALYTIC_FUNCTION(F,x,y)
{
  if (x[0]>0.4 && x[0]<0.5 && x[1]>0.4 && x[1]<0.55)
    y = sin(time()*0.1)*350.0;
  else if ((x[0] > 0.1 && x[0] < 0.2 && x[1] > 0.75 && x[1] < 0.85) || (x[0] > 0.8 && x[0] < 0.95 && x[1] > 0.05 && x[1] < 0.2))
    y = cos(time()*0.15)*(-100);
  else
    y = 0.0;
}
END_INSTATIONARY_ANALYTIC_FUNCTION


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
  y = Traits::Dirichlet; // Dirichlet
}
END_SIMPLE_BOUNDARYTYPE_FUNCTION


// dirichlet bc
INSTATIONARY_ANALYTIC_FUNCTION(G,x,y)
{
  typename Traits::DomainType center;
  for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
  center -= x;
  y = exp(-center.two_norm2());
}
END_INSTATIONARY_ANALYTIC_FUNCTION


// neumann bc
INSTATIONARY_ANALYTIC_FUNCTION(J,x,y)
{
  if (x[1]<1E-6 || x[1]>1.0-1E-6)
    {
      y = 0;
      return;
    }
  if (x[0]>1.0-1E-6 && x[1]>0.5+1E-6)
    {
      y = -5.0*sin(time()*M_PI/50);
      return;
    }
}
END_INSTATIONARY_ANALYTIC_FUNCTION


int main(int argc, char** argv) {

  try {

  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <refinement level> <coupling intensity>" << std::endl;
    return 1;
  }

  Dune::Timer timer;
  Dune::Timer totalTimer;
  timer.start();
  const int dim = 2;
  //typedef Dune::SGrid<dim,dim> BaseGrid;
  typedef Dune::YaspGrid<dim> BaseGrid;
  const Dune::FieldVector<int,dim> s(1);
  const Dune::FieldVector<double,dim> h(1.0);
  const Dune::FieldVector<bool,dim> p(false);
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
  MDGV mdgv = grid.leafView();
  SDGV sdgv0 = sdg0.leafView();
  SDGV sdgv1 = sdg1.leafView();
  sdg0.hostEntityPointer(*sdgv0.begin<0>());
  grid.startSubDomainMarking();
  for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
    {
      if (it->geometry().center()[0] > 0.5)
        grid.addToSubDomain(0,*it);
      else
        grid.addToSubDomain(1,*it);
    }
  grid.preUpdateSubDomains();
  grid.updateSubDomains();
  grid.postUpdateSubDomains();

  std::cout << "grid setup: " << timer.elapsed() << " sec" << std::endl;
  timer.reset();

  typedef MDGV::Grid::ctype DF;

  typedef Dune::PDELab::Q22DLocalFiniteElementMap<ctype,double> FEM0;
  typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM1;

  typedef FEM0::Traits::LocalFiniteElementType::Traits::
  LocalBasisType::Traits::RangeFieldType R;

  FEM0 fem0;
  FEM1 fem1;
  typedef Dune::PDELab::NoConstraints NOCON;
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;

  CON con;

  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM0,NOCON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS0;

  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM1,NOCON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS1;

  typedef GFS0::ConstraintsContainer<R>::Type C;
  C cg;

  GFS0 gfs0(sdgv0,fem0);
  GFS1 gfs1(sdgv1,fem1);

  typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,GFS0,GFS1> MultiGFS;

  MultiGFS multigfs(grid,gfs0,gfs1);

  std::cout << "function space setup: " << timer.elapsed() << " sec" << std::endl;
  timer.reset();

  typedef B<MDGV> BType;
  BType b(mdgv);

  typedef F<MDGV,R,double> FType;
  FType f(mdgv);

  typedef G<MDGV,R,double> GType;
  GType g(mdgv);

  typedef J<MDGV,R,double> JType;
  JType j(mdgv);

  // different integration order
  typedef InstationaryPoisson<FType,BType,JType,4> LOP0;
  LOP0 lop0(f,b,j);
  typedef InstationaryPoisson<FType,BType,JType,2> LOP1;
  LOP1 lop1(f,b,j);


  typedef SimpleTimeOperator TOP;
  TOP top0(4);
  TOP top1(2);

  typedef MDGV::IndexSet::SubDomainSet SDS;
  typedef Dune::PDELab::MultiDomain::EqualsSubDomains<SDS> EC;

  EC ec0(0);
  EC ec1(1);

  typedef Dune::PDELab::MultiDomain::InstationarySubProblem<double,MultiGFS,CON,MultiGFS,CON,LOP0,TOP,EC,0> SubProblem0;
  typedef Dune::PDELab::MultiDomain::InstationarySubProblem<double,MultiGFS,CON,MultiGFS,CON,LOP1,TOP,EC,1> SubProblem1;
  SubProblem0 sp0(con,con,lop0,top0,ec0);
  SubProblem1 sp1(con,con,lop1,top1,ec1);

  SubProblem0::Traits::LocalTrialFunctionSpace
    splfs0(sp0,sp0.trialGridFunctionSpaceConstraints());
  SubProblem1::Traits::LocalTrialFunctionSpace
    splfs1(sp1,sp1.trialGridFunctionSpaceConstraints());

  ProportionalFlowCoupling proportionalFlowCoupling(atof(argv[2]));

  typedef Dune::PDELab::MultiDomain::InstationaryCoupling<double,SubProblem0,SubProblem1,ProportionalFlowCoupling> Coupling;
  Coupling coupling(sp0,sp1,proportionalFlowCoupling);

  std::cout << "subproblem / coupling setup: " << timer.elapsed() << " sec" << std::endl;
  timer.reset();

  constraints(b,multigfs,cg,b,splfs0,b,splfs1);

  std::cout << "constraints evaluation: " << timer.elapsed() << " sec" << std::endl;
  timer.reset();

  // make coefficent Vector and initialize it from a function
  typedef MultiGFS::VectorContainer<R>::Type V;
  V xold(multigfs);
  xold = 0.0;
  Dune::PDELab::MultiDomain::interpolate(multigfs,xold,g,splfs0,g,splfs1);

  std::cout << "interpolation: " << timer.elapsed() << " sec" << std::endl;
  std::cout << xold.size() << " dof total, " << cg.size() << " dof constrained" << std::endl;
  timer.reset();

  typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,0> SGFS0;
  typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,1> SGFS1;
  SGFS0 sgfs0(multigfs);
  SGFS1 sgfs1(multigfs);

  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

  typedef Dune::PDELab::MultiDomain::InstationaryMultiDomainGridOperatorSpace<double,V,MultiGFS,MultiGFS,MBE,SubProblem0,SubProblem1,Coupling> MultiGOS;

  MultiGOS multigos(multigfs,multigfs,cg,cg,sp0,sp1,coupling);

  std::cout << "operator space setup: " << timer.elapsed() << " sec" << std::endl;
  timer.reset();

  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,false);

  typedef Dune::PDELab::StationaryLinearProblemSolver<MultiGOS,LS,V> PDESOLVER;
  PDESOLVER pdesolver(multigos,ls,1e-10);

  Dune::PDELab::Alexander2Parameter<double> method;
  Dune::PDELab::OneStepMethod<double,MultiGOS,PDESOLVER,V,V> osm(method,multigos,pdesolver);
  osm.setVerbosityLevel(2);

  Dune::PDELab::FilenameHelper fn0("instationaryexplicitlycoupledpoisson-right");
  {
    typedef Dune::PDELab::DiscreteGridFunction<SGFS0,V> U0DGF;
    U0DGF u0dgf(sgfs0,xold);
    Dune::SubsamplingVTKWriter<SDGV> vtkwriter(sdgv0,2);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"c0"));
    vtkwriter.write(fn0.getName(),Dune::VTKOptions::binaryappended);
    fn0.increment();
  }
  Dune::PDELab::FilenameHelper fn1("instationaryexplicitlycoupledpoisson-left");
  {
    typedef Dune::PDELab::DiscreteGridFunction<SGFS1,V> U1DGF;
    U1DGF u1dgf(sgfs1,xold);
    Dune::VTKWriter<SDGV> vtkwriter(sdgv1,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"c1"));
    vtkwriter.write(fn1.getName(),Dune::VTKOptions::binaryappended);
    fn1.increment();
  }

  V xnew(xold);
  double time = 0;
  double dt = 1.0;

  for (int i = 0; i < 100; ++i)
    {
      osm.apply(time,dt,xold,xnew);
      xold = xnew;
      time += dt;

      {
        typedef Dune::PDELab::DiscreteGridFunction<SGFS0,V> U0DGF;
        U0DGF u0dgf(sgfs0,xold);
        Dune::SubsamplingVTKWriter<SDGV> vtkwriter(sdgv0,2);
        vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"c0"));
        vtkwriter.write(fn0.getName(),Dune::VTKOptions::binaryappended);
        fn0.increment();
      }

      {
        typedef Dune::PDELab::DiscreteGridFunction<SGFS1,V> U1DGF;
        U1DGF u1dgf(sgfs1,xold);
        Dune::VTKWriter<SDGV> vtkwriter(sdgv1,Dune::VTKOptions::conforming);
        vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"c1"));
        vtkwriter.write(fn1.getName(),Dune::VTKOptions::binaryappended);
        fn1.increment();
      }
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
