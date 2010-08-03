#include "config.h"

#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemgridfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/multidomaingridoperatorspace.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/multidomain/coupling.hh>

#include<typeinfo>

#include "functionmacros.hh"
#include "proportionalflowcoupling.hh"

// source term
SIMPLE_ANALYTIC_FUNCTION(F,x,y)
{
  if (x[0]>0.4 && x[0]<0.5 && x[1]>0.4 && x[1]<0.55)
    y = 350.0;
  else if ((x[0] > 0.1 && x[0] < 0.2 && x[1] > 0.75 && x[1] < 0.85) || (x[0] > 0.8 && x[0] < 0.95 && x[1] > 0.05 && x[1] < 0.2))
    y = -100;
  else
    y = 0.0;
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

  typedef F<MDGV,R> FType;
  FType f(mdgv);

  typedef G<MDGV,R> GType;
  GType g(mdgv);

  typedef J<MDGV,R> JType;
  JType j(mdgv);

  typedef Dune::PDELab::Poisson<FType,BType,JType,2> LOP;
  LOP lop(f,b,j);

  typedef MDGV::IndexSet::SubDomainSet SDS;
  typedef Dune::PDELab::MultiDomain::EqualsSubDomains<SDS> EC;

  EC ec0(0);
  EC ec1(1);

  typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,CON,MultiGFS,CON,LOP,EC,GFS0> SubProblem0;
  typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,CON,MultiGFS,CON,LOP,EC,GFS1> SubProblem1;
  SubProblem0 sp0(con,con,lop,ec0);
  SubProblem1 sp1(con,con,lop,ec1);

  SubProblem0::Traits::LocalTrialFunctionSpace
    splfs0(sp0,sp0.trialGridFunctionSpaceConstraints());
  SubProblem1::Traits::LocalTrialFunctionSpace
    splfs1(sp1,sp1.trialGridFunctionSpaceConstraints());

  ProportionalFlowCoupling proportionalFlowCoupling(atof(argv[2]));

  typedef Dune::PDELab::MultiDomain::Coupling<SubProblem0,SubProblem1,ProportionalFlowCoupling> Coupling;
  Coupling coupling(sp0,sp1,proportionalFlowCoupling);

  std::cout << "subproblem / coupling setup: " << timer.elapsed() << " sec" << std::endl;
  timer.reset();

  constraints(b,multigfs,cg,b,splfs0,b,splfs1);

  std::cout << "constraints evaluation: " << timer.elapsed() << " sec" << std::endl;
  timer.reset();

  // make coefficent Vector and initialize it from a function
  typedef MultiGFS::VectorContainer<R>::Type V;
  V x0(multigfs);
  x0 = 0.0;
  Dune::PDELab::MultiDomain::interpolate(multigfs,x0,g,splfs0,g,splfs1);

  Dune::PDELab::set_shifted_dofs(cg,0.0,x0);

  std::cout << "interpolation: " << timer.elapsed() << " sec" << std::endl;
  std::cout << x0.size() << " dof total, " << cg.size() << " dof constrained" << std::endl;
  timer.reset();

  typedef Dune::PDELab::MultiDomain::TypeBasedGridFunctionSubSpace<MultiGFS,GFS0> SGFS0;
  typedef Dune::PDELab::MultiDomain::TypeBasedGridFunctionSubSpace<MultiGFS,GFS1> SGFS1;
  SGFS0 sgfs0(multigfs);
  SGFS1 sgfs1(multigfs);

  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

  typedef Dune::PDELab::MultiDomain::MultiDomainGridOperatorSpace<MultiGFS,MultiGFS,MBE,SubProblem0,SubProblem1,Coupling> MultiGOS;

  MultiGOS multigos(multigfs,multigfs,cg,cg,sp0,sp1,coupling);

  std::cout << "operator space setup: " << timer.elapsed() << " sec" << std::endl;
  timer.reset();

  typedef MultiGOS::MatrixContainer<R>::Type M;
  M m(multigos);
  m = 0.0;


  std::cout << "matrix construction: " << timer.elapsed() << " sec" << std::endl;
  timer.reset();

  multigos.jacobian(x0,m);

  std::cout << "jacobian evaluation: " << timer.elapsed() << " sec" << std::endl;
  timer.reset();

  V r(multigfs);

  r = 0.0;

  multigos.residual(x0,r);
  std::cout << "residual evaluation: " << timer.elapsed() << " sec" << std::endl;
  timer.reset();

  Dune::MatrixAdapter<M,V,V> opa(m);
  Dune::SeqSSOR<M,V,V> ssor(m,1,1.0);
  Dune::CGSolver<V> solver(opa,ssor,1e-10,5000,2);
  Dune::InverseOperatorResult stat;

  r *= -1.0;

  V x(multigfs,0.0);
  solver.apply(x,r,stat);

  std::cout << "linear system: " << timer.elapsed() << " sec" << std::endl;
  std::cout << "total time: " << totalTimer.elapsed() << " sec" << std::endl;
  timer.reset();

  x += x0;

  {
    typedef Dune::PDELab::DiscreteGridFunction<SGFS0,V> DGF;
    DGF dgf(sgfs0,x);
    Dune::SubsamplingVTKWriter<SDGV> vtkwriter(sdgv0,2);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.write("poisson-right",Dune::VTKOptions::ascii);
  }

  {
    typedef Dune::PDELab::DiscreteGridFunction<SGFS1,V> DGF;
    DGF dgf(sgfs1,x);
    Dune::VTKWriter<SDGV> vtkwriter(sdgv1,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.write("poisson-left",Dune::VTKOptions::ascii);
  }

  std::cout << "output I/O: " << timer.elapsed() << " sec" << std::endl;

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
