#include "config.h"

#include<typeinfo>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/localoperator/poisson.hh>

#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinggridfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinglocalfunctionspace.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/multidomain/coupling.hh>


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



class EnrichedCouplingOperator :
  public Dune::PDELab::MultiDomain::CouplingOperatorDefaultFlags,
  public Dune::PDELab::MultiDomain::NumericalJacobianEnrichedCoupling<EnrichedCouplingOperator>,
  public Dune::PDELab::MultiDomain::NumericalJacobianApplyCoupling<EnrichedCouplingOperator>,
  public Dune::PDELab::MultiDomain::FullEnrichedCouplingPattern,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{

public:

  static const bool doAlphaEnrichedCoupling = true;
  static const bool doPatternEnrichedCoupling = true;

  template<typename IG, typename LFSU1, typename LFSU2, typename X, typename LFSV1, typename LFSV2,
           typename CouplingLFSU, typename CouplingLFSV, typename R>
  void alpha_enriched_coupling
  ( const IG& ig,
    const LFSU1& lfsu_s, const X& x_s, const LFSV1& lfsv_s,
    const LFSU2& lfsu_n, const X& x_n, const LFSV2& lfsv_n,
    const CouplingLFSU& lfsu_c, const X& x_c, const CouplingLFSV& lfsv_c,
    R& r_s, R& r_n, R& r_c) const
  {
    for (auto it = r_s.begin(); it != r_s.end(); ++it)
      *it += 5;
    for (auto it = r_n.begin(); it != r_n.end(); ++it)
      *it += 5;
    for (auto it = r_c.begin(); it != r_c.end(); ++it)
      *it += 10;
  }

};

template<int dim>
struct PseudoGV
{
  static const int dimension = dim;
};

int main(int argc, char** argv) {

  try {

  Dune::MPIHelper::instance(argc,argv);

  const int dim = 2;

#ifdef UGGRID
  typedef Dune::UGGrid<dim> BaseGrid;
  BaseGrid baseGrid(500);
  std::vector<int> boundaryIndexToPhysicalGroup, elementIndexToPhysicalGroup;
  Dune::GmshReader<BaseGrid> gmshreader;
  gmshreader.read(baseGrid,"gmshtest.msh",boundaryIndexToPhysicalGroup,elementIndexToPhysicalGroup,true,false);
#else
  typedef Dune::YaspGrid<dim> BaseGrid;
  const Dune::FieldVector<double,dim> h(1.0);
  const Dune::array<int,dim> s = { {1,1} };
  const std::bitset<dim> p(false);
  BaseGrid baseGrid(h,s,p,0);
  baseGrid.globalRefine(2);
#endif

  typedef BaseGrid::LeafGridView GV;

  GV gv = baseGrid.leafGridView();
  const GV::IndexSet& is = gv.indexSet();

  typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,4> > Grid;
  Grid grid(baseGrid,false);
  typedef Grid::SubDomainGrid SubDomainGrid;
  SubDomainGrid& sdg0 = grid.subDomain(0);
  SubDomainGrid& sdg1 = grid.subDomain(1);
  SubDomainGrid& sdg2 = grid.subDomain(2);
  typedef Grid::ctype ctype;
  typedef Grid::LeafGridView MDGV;
  typedef SubDomainGrid::LeafGridView SDGV;
  MDGV mdgv = grid.leafGridView();
  SDGV sdgv0 = sdg0.leafGridView();
  SDGV sdgv1 = sdg1.leafGridView();
  SDGV sdgv2 = sdg2.leafGridView();
  grid.startSubDomainMarking();
  for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
    {
#ifdef UGGRID
      grid.addToSubDomain(elementIndexToPhysicalGroup[mdgv.indexSet().index(*it)],*it);
#else
      Dune::FieldVector<ctype,dim> center = it->geometry().center();
      if (center[0] > 0.5)
        grid.addToSubDomain(0,*it);
      if (center[0] < 0.5)
        grid.addToSubDomain(1,*it);
      if (center[1] > 0.3 && center[1] < 0.7)
        grid.addToSubDomain(2,*it);
#endif
    }
  grid.preUpdateSubDomains();
  grid.updateSubDomains();
  grid.postUpdateSubDomains();

  typedef MDGV::Grid::ctype DF;

  // we need something to feed into the QkFEM, only needs to return a correct dimension
  using CouplingGV = PseudoGV<dim-1>;

  typedef Dune::PDELab::QkLocalFiniteElementMap<CouplingGV,DF,double,1> COUPLINGFEM;

  typedef COUPLINGFEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  COUPLINGFEM couplingfem({});
  typedef Dune::PDELab::NoConstraints NOCON;
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;

  NOCON con;

  typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> EC;
  typedef Dune::PDELab::MultiDomain::SubDomainSupersetCondition<Grid> SC;

  EC c0({0});
  EC c1({1});

  typedef Dune::PDELab::MultiDomain::SubProblemSubProblemInterface<MDGV,EC,EC> Pred;

  Pred pred(mdgv,c0,c1);

  typedef Dune::PDELab::MultiDomain::CouplingGridFunctionSpace<MDGV,COUPLINGFEM,Pred,NOCON,VBE> CouplingGFS;

  typedef double RF;

  typedef CouplingGFS::ConstraintsContainer<RF>::Type C;
  C cg;

  CouplingGFS couplinggfs(mdgv,couplingfem,pred);
  std::cout << couplinggfs.size() << std::endl;

  typedef CouplingGFS::VectorContainer<double>::Type Vec;
  Vec global(couplinggfs,0.0);

  typedef Dune::PDELab::MultiDomain::CouplingLocalFunctionSpaceNode<CouplingGFS> CouplingLFS;
  CouplingLFS couplinglfs(couplinggfs);

  for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
    {
      for(MDGV::IntersectionIterator iit = mdgv.ibegin(*it); iit != mdgv.iend(*it); ++iit)
        {
          if (couplinglfs.bind(couplinglfs,*iit))
            {
              std::vector<double> local(couplinglfs.size(),0.0);
              couplinglfs.vread(global,local);
              for (int i = 0; i < couplinglfs.size(); ++i)
                {
                  local[couplinglfs.localIndex(i)] = i+1;
                }
              couplinglfs.debug();
              couplinglfs.vadd(local,global);
            }
        }
    }

  std::cout << std::endl;

  for(auto it = global.begin(); it != global.end(); ++it)
    std::cout << *it << " ";
  std::cout << std::endl;

  typedef Dune::PDELab::PkLocalFiniteElementMap<SDGV,DF,RF,2> FEM0;
  FEM0 fem0(sdgv0);
  typedef Dune::PDELab::PkLocalFiniteElementMap<SDGV,DF,RF,1> FEM1;
  FEM1 fem1(sdgv1);

  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM0,NOCON,VBE> GFS0;
  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM1,NOCON,VBE> GFS1;

  GFS0 gfs0(sdgv0,fem0);
  GFS1 gfs1(sdgv1,fem1);

  typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,GFS0,GFS1,CouplingGFS> MultiGFS;

  MultiGFS multigfs(grid,gfs0,gfs1,couplinggfs);

  typedef MultiGFS::LocalFunctionSpace LFS;
  typedef MultiGFS::CouplingLocalFunctionSpace CLFS;
  LFS lfs(multigfs);
  CLFS clfs(multigfs);

  for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
    {
      lfs.bind(*it);
      lfs.debug();
      for(MDGV::IntersectionIterator iit = mdgv.ibegin(*it); iit != mdgv.iend(*it); ++iit)
        {
          clfs.bind(*iit);
          clfs.debug();
        }
    }


  typedef B<MDGV> BType;
  BType b(mdgv);

  typedef F<MDGV,R> FType;
  FType f(mdgv);

  typedef G<MDGV,R> GType;
  GType g(mdgv);

  typedef J<MDGV,R> JType;
  JType j(mdgv);

  typedef Dune::PDELab::Poisson<FType,BType,JType,4> LOP;
  LOP lop(f,b,j);

  typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP,EC,GFS0> SubProblem0;
  typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP,EC,GFS1> SubProblem1;

  SubProblem0 sp0(con,con,lop,c0);
  SubProblem1 sp1(con,con,lop,c1);

  EnrichedCouplingOperator enrichedCouplingOperator;

  typedef Dune::PDELab::MultiDomain::EnrichedCoupling<SubProblem0,SubProblem1,EnrichedCouplingOperator,2> Coupling;
  Coupling coupling(sp0,sp1,enrichedCouplingOperator);

  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

  typedef Dune::PDELab::MultiDomain::MultiDomainGridOperatorSpace<MultiGFS,MultiGFS,MBE,SubProblem0,SubProblem1,Coupling> MultiGOS;

  MultiGOS multigos(multigfs,multigfs,cg,cg,sp0,sp1,coupling);

  typedef MultiGOS::MatrixContainer<R>::Type M;
  M m(multigos);
  m = 0.0;

  typedef MultiGFS::VectorContainer<R>::Type V;

  V u(multigfs,0.0);
  V r(u);

  multigos.residual(u,r);
  multigos.jacobian(u,m);

  std::cout << r << std::endl;

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
