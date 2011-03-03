#include "config.h"

#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
//#include <dune/pdelab/multidomain/constraints.hh>
//#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/multidomain/coupling.hh>

#include "functionmacros.hh"

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


struct DirichletBoundary :
  public Dune::PDELab::DirichletConstraintsParameters
{
  template<typename I>
  bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x)
  {
    Dune::FieldVector<typename I::ctype,I::dimension>
      xg = ig.geometry().global(x);

    if (!ig.boundary())
      {
        return false;
      }

    if (xg[1]<1E-6 || xg[1]>1.0-1E-6)
      {
        return false;
      }

    if (xg[0]>1.0-1E-6 && xg[1]>0.5+1E-6)
      {
        return false;
      }

    return true;
  }
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

template<typename MDLFS, typename SDS, typename... SubProblems>
void instantiateLocalFunctionSpaces(const MDLFS& mdlfs, const SDS& sds, SubProblems&... subProblems)
{
}

template<typename MDLFS, typename SDS, typename SubProblem, typename... SubProblems>
void instantiateLocalFunctionSpaces(const MDLFS& mdlfs, const SDS& sds, SubProblem& subProblem, SubProblems&... subProblems)
{
  if (subProblem.appliesTo(sds))
    {
      typename SubProblem::Traits::LocalTrialFunctionSpace lfs(mdlfs,subProblem,subProblem.trialGridFunctionSpaceConstraints());
      /*int status;
      std::unique_ptr<char> name(abi::__cxa_demangle(typeid(lfs).name(),0,0,&status));
      std::cout << name.get() << std::endl;*/
    }
  instantiateLocalFunctionSpaces(mdlfs,sds,subProblems...);
}

int main(int argc, char** argv) {

  try {

    Dune::MPIHelper::instance(argc,argv);

    const int dim = 2;
    typedef Dune::YaspGrid<dim> BaseGrid;
    const Dune::FieldVector<int,dim> s(1);
    const Dune::FieldVector<double,dim> h(1.0);
    const Dune::FieldVector<bool,dim> p(false);
    BaseGrid baseGrid(h,s,p,0);
    baseGrid.globalRefine(1);
    typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,8> > Grid;
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

    typedef Dune::PDELab::Q22DLocalFiniteElementMap<ctype,double> FEM0;
    typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM1;

    typedef FEM0::Traits::FiniteElementType::Traits::
      LocalBasisType::Traits::RangeFieldType R;

    FEM0 fem0;
    FEM1 fem1;
    typedef Dune::PDELab::NoConstraints NOCON;
    typedef Dune::PDELab::ISTLVectorBackend<1> VBE;

    typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM0,NOCON,
      Dune::PDELab::ISTLVectorBackend<1> > GFS0;

    typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM1,NOCON,
      Dune::PDELab::ISTLVectorBackend<1> > GFS1;

    typedef Dune::PDELab::GridFunctionSpace<MDGV,FEM0,NOCON,VBE> GFS2;

    typedef Dune::PDELab::PowerGridFunctionSpace<GFS0,2,Dune::PDELab::GridFunctionSpaceLexicographicMapper> PowerGFS;

    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,GFS0,PowerGFS> CompositeGFS;

    typedef GFS0::ConstraintsContainer<R>::Type C;
    C cg;

    GFS0 gfs0(sdgv0,fem0);
    GFS1 gfs1(sdgv1,fem1);
    GFS2 gfs2(mdgv,fem0);
    PowerGFS powergfs(gfs0);
    CompositeGFS compositegfs(gfs0,powergfs);

    typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,GFS0,GFS1,GFS2,PowerGFS,CompositeGFS> MultiGFS;

    MultiGFS multigfs(grid,gfs0,gfs1,gfs2,powergfs,compositegfs);

    Dune::PDELab::LocalFunctionSpace<MultiGFS> multilfs(multigfs);

    multilfs.bind(*mdgv.begin<0>());


    typedef DirichletBoundary BType;
    BType b;

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

    typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP,Condition,GFS0> SubProblem0;
    typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP,Condition,GFS1> SubProblem1;
    typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP,Condition,GFS0,PowerGFS> SubProblem2;
    typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP,Condition,PowerGFS> SubProblem3;
    typedef Dune::PDELab::MultiDomain::TypeBasedSubProblem<MultiGFS,NOCON,MultiGFS,NOCON,LOP,Condition,CompositeGFS> SubProblem4;

    NOCON nocon;
    SubProblem0 sp0(nocon,nocon,lop,c0);
    SubProblem1 sp1(nocon,nocon,lop,c1);
    SubProblem2 sp2(nocon,nocon,lop,c0);
    SubProblem3 sp3(nocon,nocon,lop,c0);
    SubProblem4 sp4(nocon,nocon,lop,c0);

    SubProblem0::Traits::LocalTrialFunctionSpace
      splfs0(multilfs,sp0,sp0.trialGridFunctionSpaceConstraints());
    SubProblem1::Traits::LocalTrialFunctionSpace
      splfs1(multilfs,sp1,sp1.trialGridFunctionSpaceConstraints());
    SubProblem2::Traits::LocalTrialFunctionSpace
      splfs2(multilfs,sp2,sp2.trialGridFunctionSpaceConstraints());
    SubProblem3::Traits::LocalTrialFunctionSpace
      splfs3(multilfs,sp3,sp3.trialGridFunctionSpaceConstraints());
    SubProblem4::Traits::LocalTrialFunctionSpace
      splfs4(multilfs,sp4,sp4.trialGridFunctionSpaceConstraints());

    for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
      {
        multilfs.bind(*it);
        auto sds = mdgv.indexSet().subDomains(*it);
        instantiateLocalFunctionSpaces(multilfs,sds,sp0,sp1,sp2,sp3,sp4);
        if (sp0.appliesTo(sds))
          {
            SubProblem0::Traits::LocalTrialFunctionSpace lfs(multilfs,sp0,sp0.trialGridFunctionSpaceConstraints());
            lfs.finiteElement();
            lfs.localVectorSize();
          }
        if (sp1.appliesTo(sds))
          {
            SubProblem1::Traits::LocalTrialFunctionSpace lfs(multilfs,sp1,sp1.trialGridFunctionSpaceConstraints());
            lfs.finiteElement();
            lfs.localVectorSize();
          }
        if (sp2.appliesTo(sds))
          {
            SubProblem2::Traits::LocalTrialFunctionSpace lfs(multilfs,sp2,sp2.trialGridFunctionSpaceConstraints());
            lfs.child<0>().finiteElement();
            lfs.child<1>().child(0).finiteElement();
            lfs.child<0>().localVectorSize();
            lfs.child<1>().child(1).localVectorSize();
          }
        if (sp3.appliesTo(sds))
          {
            SubProblem3::Traits::LocalTrialFunctionSpace lfs(multilfs,sp3,sp3.trialGridFunctionSpaceConstraints());
            lfs.child(1).finiteElement();
            lfs.child(0).localVectorSize();
          }
        if (sp4.appliesTo(sds))
          {
            SubProblem4::Traits::LocalTrialFunctionSpace lfs(multilfs,sp4,sp4.trialGridFunctionSpaceConstraints());
            lfs.child<0>().finiteElement();
            lfs.child<1>().child(0).finiteElement();
            lfs.child<0>().localVectorSize();
            lfs.child<1>().child(1).localVectorSize();
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
