#include "config.h"

#include <dune/grid/sgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemgridfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/multidomaingridoperatorspace.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/poisson.hh>

#include<typeinfo>

#define SIMPLE_ANALYTIC_FUNCTION(_NAME_,_CODE_) \
template<typename GV, typename RF> \
class _NAME_ \
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>, \
                                                  F<GV,RF> > \
{ \
public: \
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits; \
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F<GV,RF> > BaseT; \
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
                                                  B<GV> > \
{ \
  const GV& gv; \
\
public: \
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits; \
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B<GV> > BaseT; \
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

  typedef typename MDGV::Grid::ctype DF;

  typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM;

  typedef typename FEM::Traits::LocalFiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  FEM fem;
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;

  CON con;

  typedef Dune::PDELab::GridFunctionSpace<MDGV,FEM,CON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS;

  typedef GFS::ConstraintsContainer<R>::Type C;
  C cg;

  GFS gfs(mdgv,fem,con);

  typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,GFS> MultiGFS;

  MultiGFS multigfs(grid,gfs);

  typedef typename MultiGFS::VectorContainer<R>::Type V;
  V x0(multigfs);
  x0 = 0.0;

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

  typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,C,MultiGFS,C,LOP,EC,0> SubProblem;
  SubProblem sp0(cg,cg,lop,ec0);
  SubProblem sp1(cg,cg,lop,ec1);

  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

  typedef Dune::PDELab::MultiDomain::MultiDomainGridOperatorSpace<MultiGFS,MultiGFS,MBE,SubProblem,SubProblem> MultiGOS;

  MultiGOS multigos(multigfs,multigfs,sp0,sp1);

  typedef MultiGOS::MatrixContainer<R>::Type M;
  M m(multigos);
  m = 0.0;

  multigos.jacobian(x0,m);

}
