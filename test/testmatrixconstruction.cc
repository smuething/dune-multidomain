#include "config.h"

#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/finiteelementmap/q12dfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/multidomaingridoperatorspace.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/localoperator/poisson.hh>
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

  Dune::MPIHelper::instance(argc,argv);

  const int dim = 2;
  typedef Dune::YaspGrid<dim> BaseGrid;
  const Dune::FieldVector<int,dim> s(1);
  const Dune::FieldVector<double,dim> h(1.0);
  const Dune::FieldVector<bool,dim> p(false);
  BaseGrid baseGrid(h,s,p,0);
  baseGrid.globalRefine(6);
  typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,4> > Grid;
  Grid grid(baseGrid,false);
  typedef Grid::SubDomainGrid SubDomainGrid;
  SubDomainGrid& sdg0 = grid.subDomain(0);
  SubDomainGrid& sdg1 = grid.subDomain(1);
  SubDomainGrid& sdg2 = grid.subDomain(2);
  typedef Grid::ctype ctype;
  typedef Grid::LeafGridView MDGV;
  typedef SubDomainGrid::LeafGridView SDGV;
  MDGV mdgv = grid.leafView();
  SDGV sdgv0 = sdg0.leafView();
  SDGV sdgv1 = sdg1.leafView();
  SDGV sdgv2 = sdg2.leafView();
  grid.startSubDomainMarking();
  for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
    {
      Dune::FieldVector<ctype,dim> center = it->geometry().center();
      if (center[0] > 0.4)
        grid.addToSubDomain(0,*it);
      if (center[0] < 0.6)
        grid.addToSubDomain(1,*it);
      if (center[1] > 0.3 && center[1] < 0.7)
        grid.addToSubDomain(2,*it);
    }
  grid.preUpdateSubDomains();
  grid.updateSubDomains();
  grid.postUpdateSubDomains();

  typedef MDGV::Grid::ctype DF;

  typedef Dune::PDELab::Q22DLocalFiniteElementMap<ctype,double> FEM;
  typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM0;
  typedef Dune::PDELab::Q12DLocalFiniteElementMap<ctype,double> FEM1;
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<ctype,double> FEM2;

  typedef FEM0::Traits::LocalFiniteElementType::Traits::
  LocalBasisType::Traits::RangeFieldType R;

  FEM fem;
  FEM0 fem0;
  FEM1 fem1;
  FEM2 fem2;
  typedef Dune::PDELab::NoConstraints NOCON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;

  NOCON con;

  typedef Dune::PDELab::GridFunctionSpace<MDGV,FEM,NOCON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS;

  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM0,NOCON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS0;

  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM1,NOCON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS1;

  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM2,NOCON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS2;


  typedef GFS0::ConstraintsContainer<R>::Type C;
  C cg;

  GFS gfs(mdgv,fem);
  GFS0 gfs0(sdgv0,fem0);
  GFS1 gfs1(sdgv1,fem1);
  GFS2 gfs2(sdgv2,fem2);

  typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,GFS,GFS0,GFS1,GFS2> MultiGFS;

  MultiGFS multigfs(grid,gfs,gfs0,gfs1,gfs2);


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

  typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> EC;
  typedef Dune::PDELab::MultiDomain::SubDomainSupersetCondition<Grid> SC;

  SC c0;
  SC c1(2);
  EC c2(0,1);
  EC c3(1);
  SC c4(0);

  typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP,SC,GFS> SubProblem0;
  typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP,SC,GFS2> SubProblem1;
  typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP,EC,GFS0,GFS1> SubProblem2;
  typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP,EC,GFS1> SubProblem3;
  typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP,SC,GFS0> SubProblem4;

  SubProblem0 sp0(con,con,lop,c0);
  SubProblem1 sp1(con,con,lop,c1);
  SubProblem2 sp2(con,con,lop,c2);
  SubProblem3 sp3(con,con,lop,c3);
  SubProblem4 sp4(con,con,lop,c4);

  ProportionalFlowCoupling proportionalFlowCoupling(atof(argv[2]));

  typedef Dune::PDELab::MultiDomain::Coupling<SubProblem0,SubProblem1,ProportionalFlowCoupling> Coupling;
  Coupling coupling(sp0,sp1,proportionalFlowCoupling);

  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

  typedef Dune::PDELab::MultiDomain::MultiDomainGridOperatorSpace<MultiGFS,MultiGFS,MBE,SubProblem0,SubProblem1,Coupling,SubProblem2,SubProblem4,SubProblem3> MultiGOS;

  MultiGOS multigos(multigfs,multigfs,cg,cg,sp0,sp1,coupling,sp2,sp4,sp3);

  typedef MultiGOS::MatrixContainer<R>::Type M;
  M m(multigos);
  m = 0.0;

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
