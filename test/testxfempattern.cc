#include "config.h"

#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/finiteelementmap/q12dfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/common/function.hh>

#include<typeinfo>

struct PatternOperator
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::FullVolumePattern
{
  static const bool doPatternVolume = true;
};

int main(int argc, char** argv) {

  try {

  if (argc < 2)
    {
      std::cerr << "Usage: " << argv[0] << " <refinement steps>" << std::endl;
      return 1;
    }

  Dune::MPIHelper::instance(argc,argv);

  const int dim = 2;
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
  typedef Grid::ctype ctype;
  typedef Grid::LeafGridView MDGV;
  typedef SubDomainGrid::LeafGridView SDGV;
  MDGV mdgv = grid.leafView();
  SDGV sdgv0 = sdg0.leafView();

  grid.startSubDomainMarking();
  for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
    {
      Dune::FieldVector<ctype,dim> center = it->geometry().center();
      if (center[1] < 0.5)
        grid.addToSubDomain(0,*it);
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
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;

  NOCON con;

  typedef Dune::PDELab::GridFunctionSpace<MDGV,FEM,NOCON,
    Dune::PDELab::ISTLVectorBackend<1> > BGFS;

  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM,NOCON,
    Dune::PDELab::ISTLVectorBackend<1> > XGFS;

  typedef BGFS::ConstraintsContainer<R>::Type C;
  C cg;

  BGFS bgfs(mdgv,fem);
  XGFS xgfs(sdgv0,fem);

  typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,BGFS,XGFS> MultiGFS;

  MultiGFS multigfs(grid,bgfs,xgfs);

  typedef PatternOperator LOP;
  LOP lop;

  typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> BPC;
  typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> XPC;

  BPC bpc;
  XPC xpc(0);

  typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,MultiGFS,LOP,BPC,BGFS> BProblem;
  typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,MultiGFS,LOP,XPC,BGFS,XGFS> XProblem;

  BProblem bproblem(lop,bpc);
  XProblem xproblem(lop,xpc);

  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

  typedef Dune::PDELab::MultiDomain::GridOperator<
    MultiGFS,MultiGFS,
    MBE,
    R,R,R,
    C,C,
    BProblem,
    XProblem
    > GO;

  GO go(multigfs,multigfs,
        cg,cg,
        bproblem,
        xproblem);

  typedef GO::Traits::Jacobian M;
  M m(go);
  m = 0.0;

  for (int i = 0; i < m.N(); ++i)
    {
      for (int j = 0; j < m.M(); ++j)
        std::cout << (m.exists(i,j) ? "X " : ". ");
      std::cout << std::endl;
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
