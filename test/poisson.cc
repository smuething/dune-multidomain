#include "config.h"

#include <dune/grid/sgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemgridfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/multidomaingridoperatorspace.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

#include<typeinfo>

#define SIMPLE_ANALYTIC_FUNCTION(_NAME_,_CODE_) \
template<typename GV, typename RF> \
class _NAME_ \
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>, \
                                                  _NAME_<GV,RF> > \
{ \
public: \
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits; \
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,_NAME_<GV,RF> > BaseT; \
\
  _NAME_ (const GV& gv) : BaseT(gv) {} \
  inline void evaluateGlobal (const typename Traits::DomainType& x, \
                              typename Traits::RangeType& y) const \
  _CODE_ \
};

#define SIMPLE_BOUNDARY_FUNCTION(_NAME_,_CODE_) \
template<typename GV> \
class _NAME_ \
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab:: \
                                                  BoundaryGridFunctionTraits<GV,int,1, \
                                                                             Dune::FieldVector<int,1> >, \
                                                  _NAME_<GV> > \
{ \
  const GV& gv; \
\
public: \
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits; \
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,_NAME_<GV> > BaseT; \
\
  _NAME_ (const GV& gv_) : gv(gv_) {} \
\
  template<typename I>\
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, \
                        const typename Traits::DomainType& x, \
                        typename Traits::RangeType& y) const \
  _CODE_ \
  inline const GV& getGridView () \
  { \
    return gv; \
  } \
};





// source term
SIMPLE_ANALYTIC_FUNCTION(F,
{
  if (x[0]>0.25 && x[0]<0.375 && x[1]>0.25 && x[1]<0.375)
    y = 50.0;
  else
    y = 0.0;
  y=0;
})


// boundary condition type
SIMPLE_BOUNDARY_FUNCTION(B,
{({
  Dune::FieldVector<typename GV::Grid::ctype,GV::dimension>
    xg = ig.geometry().global(x);

  if (xg[1]<1E-6 || xg[1]>1.0-1E-6)
    {
      y = 0; // Neumann
      return;
    }
  if (xg[0]>1.0-1E-6 && xg[1]>0.5+1E-6)
    {
      y = 0; // Neumann
      return;
    }
  y = 1; // Dirichlet
});})


// dirichlet bc
SIMPLE_ANALYTIC_FUNCTION(G,
{
  typename Traits::DomainType center;
  for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
  center -= x;
  y = exp(-center.two_norm2());
})


// neumann bc
SIMPLE_ANALYTIC_FUNCTION(J,
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
})


template<typename GV>
class B2
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::BoundaryTypeGridFunctionTraits<GV>,
                                                  B2<GV> >
{

public:

  typedef Dune::PDELab::BoundaryTypeGridFunctionTraits<GV> Traits;

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& bct) const
  {
    Dune::FieldVector<typename GV::Grid::ctype,GV::dimension>
      xg = ig.geometry().global(x);

    if (xg[0]<1E-6 || xg[0]>1.0-1E-6 || xg[1]<1E-6 || xg[1]>1.0-1E-6)
      {
        bct = Traits::Dirichlet;
        return;
      }
    bct = Traits::None; // no boundary conditions on subproblem-subproblem interface
  }

};


int main(int argc, char** argv) {

  const int dim = 2;
  typedef Dune::SGrid<dim,dim> BaseGrid;
  const int s[2] = {4,4};
  const double h[2] = {1.0,1.0};
  BaseGrid baseGrid(s,h);
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

  typedef B2<MDGV> BT;

  BT bt;
  typedef BT BType;

  typedef MDGV::Grid::ctype DF;

  typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM;

  typedef FEM::Traits::LocalFiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  FEM fem;
  typedef Dune::PDELab::NoConstraints NOCON;
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;

  CON con;

  typedef Dune::PDELab::GridFunctionSpace<MDGV,FEM,NOCON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS;

  typedef GFS::ConstraintsContainer<R>::Type C;
  C cg;

  GFS gfs(mdgv,fem);

  typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,GFS> MultiGFS;

  MultiGFS multigfs(grid,gfs);

  //typedef B<MDGV> BType;
  //BType b(mdgv);

  typedef F<MDGV,R> FType;
  FType f(mdgv);

  typedef G<MDGV,R> GType;
  GType g(mdgv);

  typedef J<MDGV,R> JType;
  JType j(mdgv);

  typedef Dune::PDELab::Poisson<FType,BType,JType,2> LOP;
  LOP lop(f,bt,j);

  typedef MDGV::IndexSet::SubDomainSet SDS;
  typedef Dune::PDELab::MultiDomain::EqualsSubDomains<SDS> EC;

  EC ec0(0);
  EC ec1(1);

  typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,CON,MultiGFS,CON,LOP,EC,0> SubProblem;
  SubProblem sp0(con,con,lop,ec0);
  SubProblem sp1(con,con,lop,ec1);

  SubProblem::Traits::LocalTrialFunctionSpace
    splfs0(multigfs,sp0,sp0.trialGridFunctionSpaceConstraints()),
    splfs1(multigfs,sp1,sp1.trialGridFunctionSpaceConstraints());

  constraints(bt,multigfs,cg,bt,splfs0,bt,splfs1);

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

  /*
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    LOP,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  typedef GOS::MatrixContainer<R>::Type M2;
  M2 m2(gos);
  m2 = 0.0;

  gos.jacobian(x0,m2);
  */
  //Dune::printmatrix(std::cout,m.base(),"","");
  //Dune::printmatrix(std::cout,m2.base(),"","");

}
