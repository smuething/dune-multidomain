#include "config.h"

#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/multidomaingridoperatorspace.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include<dune/pdelab/common/vtkexport.hh>

#include<typeinfo>

#include "functionmacros.hh"

// source term
SIMPLE_ANALYTIC_FUNCTION(F,x,y)
{
  if (x[0]>0.25 && x[0]<0.375 && x[1]>0.25 && x[1]<0.375)
    y = 50.0;
  else
    y = 0.0;
  y=0;
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

    if (argc < 2) {
      std::cerr << "Usage: " << argv[0] << " <refinement level>" << std::endl;
      return 1;
    }

    Dune::MPIHelper::instance(argc,argv);

    const int dim = 2;
    //typedef Dune::SGrid<dim,dim> BaseGrid;
    typedef Dune::YaspGrid<dim> BaseGrid;
    const Dune::FieldVector<double,dim> h(1.0);
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
        if (it->geometry().center()[0] > 0.5)
          grid.addToSubDomain(0,*it);
        else
          grid.addToSubDomain(1,*it);
      }
    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();

    typedef MDGV::Grid::ctype DF;

    typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM;

    typedef FEM::Traits::FiniteElementType::Traits::
      LocalBasisType::Traits::RangeFieldType R;

    FEM fem;
    typedef Dune::PDELab::NoConstraints NOCON;
    typedef Dune::PDELab::ConformingDirichletConstraints CON;
    typedef Dune::PDELab::ISTLVectorBackend<> VBE;

    CON con;

    typedef Dune::PDELab::GridFunctionSpace<MDGV,FEM,NOCON,VBE> GFS;

    typedef GFS::ConstraintsContainer<R>::Type C;
    C cg;

    GFS gfs(mdgv,fem);

    typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,GFS> MultiGFS;

    MultiGFS multigfs(grid,gfs);

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

    typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> Condition;

    Condition c0(0);
    Condition c1(1);

    typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,CON,MultiGFS,CON,LOP,Condition,0> SubProblem;
    SubProblem sp0(con,con,lop,c0);
    SubProblem sp1(con,con,lop,c1);

    SubProblem::Traits::LocalTrialFunctionSpace
      splfs0(sp0,sp0.trialGridFunctionSpaceConstraints()),
      splfs1(sp1,sp1.trialGridFunctionSpaceConstraints());

    constraints(b,multigfs,cg,b,splfs0,b,splfs1);

    // make coefficent Vector and initialize it from a function
    typedef MultiGFS::VectorContainer<R>::Type V;
    V x0(multigfs);
    x0 = 0.0;
    Dune::PDELab::MultiDomain::interpolate(multigfs,x0,g,splfs0,g,splfs1);
    Dune::PDELab::set_shifted_dofs(cg,0.0,x0);

    typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

    typedef Dune::PDELab::MultiDomain::MultiDomainGridOperatorSpace<MultiGFS,MultiGFS,MBE,SubProblem,SubProblem> MultiGOS;

    MultiGOS multigos(multigfs,multigfs,cg,cg,sp0,sp1);

    typedef MultiGOS::MatrixContainer<R>::Type M;
    M m(multigos);
    m = 0.0;

    multigos.jacobian(x0,m);

    V r(multigfs);

    r = 0.0;

    multigos.residual(x0,r);

    Dune::MatrixAdapter<M,V,V> opa(m);
    Dune::SeqSSOR<M,V,V> ssor(m,1,1.0);
    Dune::CGSolver<V> solver(opa,ssor,1e-10,5000,2);
    Dune::InverseOperatorResult stat;

    r *= -1.0;

    V x(multigfs,0.0);
    solver.apply(x,r,stat);

    x += x0;

    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF dgf(gfs,x);

    Dune::VTKWriter<MDGV> vtkwriter(mdgv,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.write("poisson.vtu",Dune::VTKOptions::ascii);

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
