#include "config.h"

#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/constraints/conforming.hh>
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
struct BoundaryType :
  public Dune::PDELab::DirichletConstraintsParameters
{

  enum BCType { dirichlet, neumann, none };

  template<typename I>
  BCType bc_type(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    Dune::FieldVector<typename I::ctype,I::dimension>
      xg = ig.geometry().global(x);

    if (!ig.boundary())
      {
        return none;
      }

    if (xg[1]<1E-6 || xg[1]>1.0-1E-6)
      {
        return neumann;
      }

    if (xg[0]>1.0-1E-6 && xg[1]>0.5+1E-6)
      {
        return neumann;
      }

    return dirichlet;
  }

  template<typename I>
  bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    return bc_type(ig,x) == dirichlet;
  }

  template<typename I>
  bool isNeumann(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    return bc_type(ig,x) == neumann;
  }


  template<typename R>
  void setTime(const R& r)
  {}
};



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
  const Dune::FieldVector<double,dim> h(1.0);
  const Dune::array<int,dim> s = { {1,1} };
  const std::bitset<dim> p(false);
  BaseGrid baseGrid(h,s,p,0);
  baseGrid.globalRefine(6);
  typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,4> > Grid;
  Grid grid(baseGrid,false);
  typedef Grid::SubDomainGrid SubDomainGrid;
  SubDomainGrid& sdg0 = grid.subDomain(0);
  SubDomainGrid& sdg1 = grid.subDomain(1);
  SubDomainGrid& sdg2 = grid.subDomain(2);
  typedef Grid::LeafGridView MDGV;
  typedef SubDomainGrid::LeafGridView SDGV;
  typedef MDGV::Grid::ctype DF;
  MDGV mdgv = grid.leafGridView();
  SDGV sdgv0 = sdg0.leafGridView();
  SDGV sdgv1 = sdg1.leafGridView();
  SDGV sdgv2 = sdg2.leafGridView();
  grid.startSubDomainMarking();
  for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
    {
      Dune::FieldVector<DF,dim> center = it->geometry().center();
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


  typedef Dune::PDELab::QkLocalFiniteElementMap<MDGV,DF,double,2> FEM;
  typedef Dune::PDELab::QkLocalFiniteElementMap<SDGV,DF,double,1> FEM0;
  typedef Dune::PDELab::QkLocalFiniteElementMap<SDGV,DF,double,1> FEM1;
  typedef Dune::PDELab::QkLocalFiniteElementMap<SDGV,DF,double,2> FEM2;

  typedef FEM0::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  FEM fem(mdgv);
  FEM0 fem0(sdgv0);
  FEM1 fem1(sdgv1);
  FEM2 fem2(sdgv2);
  typedef Dune::PDELab::NoConstraints NOCON;
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;

  typedef Dune::PDELab::GridFunctionSpace<MDGV,FEM,NOCON,VBE> GFS;
  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM0,NOCON,VBE> GFS0;
  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM1,NOCON,VBE> GFS1;
  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM2,NOCON,VBE> GFS2;

  typedef GFS0::ConstraintsContainer<R>::Type C;
  C cg;

  GFS gfs(mdgv,fem);
  GFS0 gfs0(sdgv0,fem0);
  GFS1 gfs1(sdgv1,fem1);
  GFS2 gfs2(sdgv2,fem2);

  typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<
    Grid,
    VBE,
    Dune::PDELab::LexicographicOrderingTag,
    GFS,
    GFS0,
    GFS1,
    GFS2
    > MultiGFS;

  MultiGFS multigfs(grid,gfs,gfs0,gfs1,gfs2);


  typedef BoundaryType BType;
  BType b;

  typedef F<MDGV,R> FType;
  FType f(mdgv);

  typedef G<MDGV,R> GType;
  GType g(mdgv);

  typedef J<MDGV,R> JType;
  JType j(mdgv);

  typedef Dune::PDELab::Poisson<FType,BType,JType> LOP;
  LOP lop(f,b,j,4);

  typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> EC;
  typedef Dune::PDELab::MultiDomain::SubDomainSupersetCondition<Grid> SC;

  SC c0;
  SC c1(2);
  EC c2(0,1);
  EC c3(1);
  SC c4(0);

  typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,SC,0> SubProblem0;
  typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,SC,3> SubProblem1;
  typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,EC,1,2> SubProblem2;
  typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,EC,2> SubProblem3;
  typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,SC,1> SubProblem4;

  SubProblem0 sp0(lop,c0);
  SubProblem1 sp1(lop,c1);
  SubProblem2 sp2(lop,c2);
  SubProblem3 sp3(lop,c3);
  SubProblem4 sp4(lop,c4);

  ProportionalFlowCoupling proportionalFlowCoupling(atof(argv[2]));

  typedef Dune::PDELab::MultiDomain::Coupling<SubProblem0,SubProblem1,ProportionalFlowCoupling> Coupling;
  Coupling coupling(sp0,sp1,proportionalFlowCoupling);

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(64);

  typedef Dune::PDELab::MultiDomain::GridOperator<
    MultiGFS,MultiGFS,
    MBE,
    R,R,R,
    C,C,
    SubProblem0,
    SubProblem1,
    Coupling,
    SubProblem2,
    SubProblem4,
    SubProblem3
    > GO;

  GO go(multigfs,multigfs,
        cg,cg,mbe,
        sp0,sp1,coupling,sp2,sp4,sp3);

  typedef GO::Traits::Jacobian M;
  M m(go);
  m = 0.0;

  std::cout << "******************** Pattern statistics ********************" << std::endl
            << m.patternStatistics() << std::endl;

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
