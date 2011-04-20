#include "config.h"

#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/backend/ovlpistlsolverbackend.hh>
#include <dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include <dune/pdelab/localoperator/l2.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/instationary/onestep.hh>

#include<typeinfo>

#include "functionmacros.hh"
#include "proportionalflowcoupling.hh"

// source term
INSTATIONARY_ANALYTIC_FUNCTION(F,x,y)
{
  if (x[0]>0.4 && x[0]<0.5 && x[1]>0.4 && x[1]<0.55)
    y = 350.0;
  else if ((x[0] > 0.1 && x[0] < 0.2 && x[1] > 0.75 && x[1] < 0.85) || (x[0] > 0.8 && x[0] < 0.95 && x[1] > 0.05 && x[1] < 0.2))
    y = -100;
  else
    y = 0.0;
}
END_INSTATIONARY_ANALYTIC_FUNCTION


struct DirichletBoundary :
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
INSTATIONARY_ANALYTIC_FUNCTION(G,x,y)
{
  typename Traits::DomainType center;
  for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
  center -= x;
  y = -exp(-center.two_norm2());
}
END_INSTATIONARY_ANALYTIC_FUNCTION


// neumann bc
INSTATIONARY_ANALYTIC_FUNCTION(J,x,y)
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
  assert(false);
}
END_INSTATIONARY_ANALYTIC_FUNCTION


int main(int argc, char** argv) {

  try {

    Dune::MPIHelper& mpihelper = Dune::MPIHelper::instance(argc,argv);

    if (argc < 5) {
      std::cerr << "Usage: " << argv[0] << " <refinement level> <coupling intensity> <dt> <end time>" << std::endl;
      return 1;
    }

    Dune::Timer timer;
    Dune::Timer totalTimer;
    timer.start();
    const int dim = 2;
    typedef Dune::YaspGrid<dim> BaseGrid;
    Dune::FieldVector<int,dim> s(1 << atoi(argv[1]));
    const Dune::FieldVector<double,dim> h(1.0);
    const Dune::FieldVector<bool,dim> p(false);
    const int overlap = 1;
    BaseGrid baseGrid(mpihelper.getCommunicator(),h,s,p,overlap);
    //baseGrid.globalRefine();
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
    typedef Grid::Partition<Dune::Interior_Partition>::LeafGridView IMDGV;
    IMDGV imdgv = grid.leafView<Dune::Interior_Partition>();
    for (IMDGV::Codim<0>::Iterator it = imdgv.begin<0>(); it != imdgv.end<0>(); ++it)
      {
        if (it->geometry().center()[0] > 0.5)
          grid.addToSubDomain(1,*it);
        else
          grid.addToSubDomain(0,*it);
      }
    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();

    std::cout << "grid setup: " << timer.elapsed() << " sec" << std::endl;
    timer.reset();

    typedef MDGV::Grid::ctype DF;

    typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM0;
    //typedef Dune::PDELab::Q22DLocalFiniteElementMap<ctype,double> FEM1;
    typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM1;

    typedef FEM0::Traits::FiniteElementType::Traits::
      LocalBasisType::Traits::RangeFieldType R;

    FEM0 fem0;
    FEM1 fem1;
    typedef Dune::PDELab::NoConstraints NOCON;
    //typedef Dune::PDELab::ConformingDirichletConstraints CON;
    typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
    typedef Dune::PDELab::ISTLVectorBackend<1> VBE;

    CON con;

    typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM0,CON,
      Dune::PDELab::ISTLVectorBackend<1> > GFS0;

    typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM1,CON,
      Dune::PDELab::ISTLVectorBackend<1> > GFS1;

    typedef GFS0::ConstraintsContainer<R>::Type C;
    C cg;

    GFS0 gfs0(sdgv0,fem0,con);
    GFS1 gfs1(sdgv1,fem1,con);

    typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,GFS0,GFS1> MultiGFS;

    MultiGFS multigfs(grid,gfs0,gfs1);

    std::cout << "function space setup: " << timer.elapsed() << " sec" << std::endl;
    timer.reset();

    typedef DirichletBoundary BType;
    BType b;

    typedef F<MDGV,R,R> FType;
    FType f(mdgv);

    typedef G<MDGV,R,R> GType;
    GType g(mdgv);

    typedef J<MDGV,R,R> JType;
    JType j(mdgv);

    typedef Dune::PDELab::InstationaryPoisson<R,FType,BType,JType,4> LOP;
    LOP lop(f,b,j);

    typedef Dune::PDELab::L2 TLOP;
    TLOP tlop(4);

    typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> Condition;

    Condition c0(0);
    Condition c1(1);

    typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,Condition,0> LeftSubProblem_dt0;
    typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,TLOP,Condition,0> LeftSubProblem_dt1;

    LeftSubProblem_dt0 left_sp_dt0(lop,c0);
    LeftSubProblem_dt1 left_sp_dt1(tlop,c0);

    typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,Condition,1> RightSubProblem_dt0;
    typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,TLOP,Condition,1> RightSubProblem_dt1;

    RightSubProblem_dt0 right_sp_dt0(lop,c1);
    RightSubProblem_dt1 right_sp_dt1(tlop,c1);

    ContinuousValueContinuousFlowCoupling<R> proportionalFlowCoupling(4,atof(argv[2]));

    typedef Dune::PDELab::MultiDomain::Coupling<LeftSubProblem_dt0,RightSubProblem_dt0,ContinuousValueContinuousFlowCoupling<R> > Coupling;
    Coupling coupling(left_sp_dt0,right_sp_dt0,proportionalFlowCoupling);

    std::cout << "subproblem / coupling setup: " << timer.elapsed() << " sec" << std::endl;
    timer.reset();

    auto constraints = Dune::PDELab::MultiDomain::constraints<R>(multigfs,
                                                                 Dune::PDELab::MultiDomain::constrainSubProblem(left_sp_dt0,b),
                                                                 Dune::PDELab::MultiDomain::constrainSubProblem(right_sp_dt0,b));

    constraints.assemble(cg);

    std::cout << "constraints evaluation: " << timer.elapsed() << " sec" << std::endl;
    timer.reset();

    typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

    typedef Dune::PDELab::MultiDomain::GridOperator<
      MultiGFS,MultiGFS,
      MBE,R,R,R,C,C,
      LeftSubProblem_dt0,
      RightSubProblem_dt0,
      Coupling
      > GridOperator_dt0;

    typedef Dune::PDELab::MultiDomain::GridOperator<
      MultiGFS,MultiGFS,
      MBE,R,R,R,C,C,
      LeftSubProblem_dt1,
      RightSubProblem_dt1
      > GridOperator_dt1;

    typedef Dune::PDELab::OneStepGridOperator<GridOperator_dt0,GridOperator_dt1> GridOperator;

    // make coefficent Vector and initialize it from a function
    typedef GridOperator::Traits::Domain V;
    V uold(multigfs);
    uold = 0.0;
    Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,uold,g,left_sp_dt0,g,right_sp_dt0);

    std::cout << "interpolation: " << timer.elapsed() << " sec" << std::endl;
    std::cout << uold.N() << " dof total, " << cg.size() << " dof constrained" << std::endl;
    timer.reset();

    GridOperator_dt0 go_dt_0(multigfs,multigfs,
                             cg,cg,
                             left_sp_dt0,
                             right_sp_dt0,
                             coupling);

    GridOperator_dt1 go_dt_1(multigfs,multigfs,
                             cg,cg,
                             left_sp_dt1,
                             right_sp_dt1);

    GridOperator gridoperator(go_dt_0,go_dt_1);

    typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,0> SGFS0;
    typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,1> SGFS1;
    SGFS0 sgfs0(multigfs);
    SGFS1 sgfs1(multigfs);

    std::cout << "operator space setup: " << timer.elapsed() << " sec" << std::endl;
    timer.reset();


    // <<<5>>> Select a linear solver backend
    typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<MultiGFS,C> LS;
    LS ls(multigfs,cg);

    //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
    //LS ls(5000,false);

    // <<<6>>> Solver for linear problem per stage
    typedef Dune::PDELab::StationaryLinearProblemSolver<GridOperator,LS,V> PDESOLVER;
    PDESOLVER pdesolver(gridoperator,ls,1e-10);

    // <<<7>>> time-stepper
    Dune::PDELab::Alexander2Parameter<R> method;               // defines coefficients
    Dune::PDELab::OneStepMethod<R,GridOperator,PDESOLVER,V,V> osm(method,gridoperator,pdesolver);
    osm.setVerbosityLevel(2);                                     // time stepping scheme

    // <<<8>>> graphics for initial guess
    Dune::PDELab::FilenameHelper fn_left("testinstationarypoisson_left");
    Dune::PDELab::FilenameHelper fn_right("testinstationarypoisson_right");

    {
      typedef Dune::PDELab::DiscreteGridFunction<SGFS0,V> DGF;
      DGF dgf(sgfs0,uold);
      Dune::VTKWriter<SDGV> vtkwriter(sdgv0,Dune::VTKOptions::conforming);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
      vtkwriter.write(fn_left.getName(),Dune::VTKOptions::ascii);
      fn_left.increment();
    }

    {
      typedef Dune::PDELab::DiscreteGridFunction<SGFS1,V> DGF;
      DGF dgf(sgfs1,uold);
      Dune::SubsamplingVTKWriter<SDGV> vtkwriter(sdgv1,2);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
      vtkwriter.write(fn_right.getName(),Dune::VTKOptions::ascii);
      fn_right.increment();
    }

    // <<<9>>> time loop
    R time = 0;
    const R dt = atof(argv[3]);
    const R tend = atof(argv[4]);
    V unew(uold);                                              // solution to be computed
    while (time<tend-1e-8) {
      // do time step
      /*
        bctype.setTime(time+dt);                                       // compute constraints
        cc.clear();                                               // for this time step
        Dune::PDELab::constraints(bctype,gfs,cc);
      */

      auto f_ = Dune::PDELab::MultiDomain::interpolateOnSubProblems(g,left_sp_dt0,g,right_sp_dt0);

      osm.apply(time,dt,uold,f_,unew);                           // do one time step

      // graphics
      {
        typedef Dune::PDELab::DiscreteGridFunction<SGFS0,V> DGF;
        DGF dgf(sgfs0,unew);
        Dune::VTKWriter<SDGV> vtkwriter(sdgv0,Dune::VTKOptions::conforming);
        vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
        vtkwriter.write(fn_left.getName(),Dune::VTKOptions::ascii);
        fn_left.increment();
      }

      {
        typedef Dune::PDELab::DiscreteGridFunction<SGFS1,V> DGF;
        DGF dgf(sgfs1,unew);
        Dune::SubsamplingVTKWriter<SDGV> vtkwriter(sdgv1,2);
        vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
        vtkwriter.write(fn_right.getName(),Dune::VTKOptions::ascii);
        fn_right.increment();
      }

      uold = unew;                                              // advance time step
      time += dt;
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
