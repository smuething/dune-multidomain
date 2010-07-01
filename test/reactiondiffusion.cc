#include "config.h"

#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/multidomain/subproblemgridfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/multidomaingridoperatorspace.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include <dune/pdelab/localoperator/diffusion.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include<typeinfo>

#include "functionmacros.hh"
#include "../../dune-pdelab-howto/examples/example05_operator.hh"

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

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <refinement level>" << std::endl;
    return 1;
  }

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

  typedef MDGV::Grid::ctype DF;

  typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM;

  typedef FEM::Traits::LocalFiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  FEM fem0, fem1;
  typedef Dune::PDELab::NoConstraints NOCON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;

  typedef Dune::PDELab::GridFunctionSpace<MDGV,FEM,NOCON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS0;

  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM,NOCON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS1;

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

  typedef Example05LocalOperator LOP1;
  LOP1 lop1(d_0,d_1,lambda,sigma,kappa,2);

  typedef MDGV::IndexSet::SubDomainSet SDS;
  typedef Dune::PDELab::MultiDomain::EqualsSubDomains<SDS> EC;

  EC ec0(0);
  EC ec1(1);

  NOCON nocon;

  typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP0,EC,0> SubProblem0;
  typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP1,EC,0,1> SubProblem1;
  SubProblem0 sp0(nocon,nocon,lop0,ec0);
  SubProblem1 sp1(nocon,nocon,lop1,ec1);

  /*
  SubProblem::Traits::LocalTrialFunctionSpace
    splfs0(multigfs,sp0,sp0.trialGridFunctionSpaceConstraints()),
    splfs1(multigfs,sp1,sp1.trialGridFunctionSpaceConstraints());

  constraints(b,multigfs,cg,b,splfs0,b,splfs1);
  */
  // make coefficent Vector and initialize it from a function
  typedef MultiGFS::VectorContainer<R>::Type V;
  V x0(multigfs);
  x0 = 0.0;
  /*
  Dune::PDELab::MultiDomain::interpolate(multigfs,x0,g,splfs0,g,splfs1);
  Dune::PDELab::set_shifted_dofs(cg,0.0,x0);
  */
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

  typedef Dune::PDELab::MultiDomain::MultiDomainGridOperatorSpace<MultiGFS,MultiGFS,MBE,SubProblem0,SubProblem1> MultiGOS;

  MultiGOS multigos(multigfs,multigfs,cg,cg,sp0,sp1);


  typedef MultiGOS::MatrixContainer<R>::Type M;
  M m(multigos);
  m = 0.0;

  for(int i = 0; i < m.base().N(); ++i) {
    for(int j = 0; j < m.base().M(); ++j) {
      std::cout << (m.base().exists(i,j) ? "X " : ". ");
    }
    std::cout << std::endl;
  }

  multigos.jacobian(x0,m);

  V r(multigfs);

  r = 0.0;

  multigos.residual(x0,r);
  /*
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
  */
}
