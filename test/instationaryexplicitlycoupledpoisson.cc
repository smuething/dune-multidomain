#include "config.h"

#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/instationarymultidomaingridoperatorspace.hh>
#include <dune/pdelab/multidomain/instationarysubproblem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/multidomain/instationarycoupling.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include "functionmacros.hh"
#include "proportionalflowcoupling.hh"
#include "simpletimeoperator.hh"
#include "instationarypoissonoperator.hh"

// source term
INSTATIONARY_ANALYTIC_FUNCTION(F,x,y)
{
  if (x[0]>0.4 && x[0]<0.5 && x[1]>0.4 && x[1]<0.55)
    y = sin(time()*0.1)*350.0;
  else if ((x[0] > 0.1 && x[0] < 0.2 && x[1] > 0.75 && x[1] < 0.85) || (x[0] > 0.8 && x[0] < 0.95 && x[1] > 0.05 && x[1] < 0.2))
    y = cos(time()*0.15)*(-100);
  else
    y = 0.0;
}
END_INSTATIONARY_ANALYTIC_FUNCTION


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
INSTATIONARY_ANALYTIC_FUNCTION(G,x,y)
{
  typename Traits::DomainType center;
  for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
  center -= x;
  y = exp(-center.two_norm2());
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
      y = -5.0*sin(time()*M_PI/50);
      return;
    }
}
END_INSTATIONARY_ANALYTIC_FUNCTION


int main(int argc, char** argv) {

  try {

    Dune::MPIHelper::instance(argc,argv);

    if (argc < 5) {
      std::cerr << "Usage: " << argv[0] << " <refinement level> <coupling intensity> <dt> <tend>" << std::endl;
      return 1;
    }

    Dune::Timer timer;
    Dune::Timer totalTimer;

    const int dim = 2;

    /*
     * Create underlying grid
     */

    typedef Dune::YaspGrid<dim> BaseGrid;
    const Dune::FieldVector<double,dim> h(1.0);
    const Dune::array<int,dim> s = { {1,1} };
    const std::bitset<dim> p(false);
    BaseGrid baseGrid(h,s,p,0);
    baseGrid.globalRefine(atoi(argv[1]));

    /*
     * Create MultiDomainGrid and obtain references to SubDomainGrids
     * Get leaf views for all grids
     */

    typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,4> > Grid;
    typedef Grid::SubDomainGrid SubDomainGrid;
    typedef Grid::ctype ctype;
    typedef Grid::LeafGridView MDGV;
    typedef SubDomainGrid::LeafGridView SDGV;

    Grid grid(baseGrid,false);
    SubDomainGrid& sdg0 = grid.subDomain(0);
    SubDomainGrid& sdg1 = grid.subDomain(1);

    MDGV mdgv = grid.leafGridView();
    SDGV sdgv0 = sdg0.leafGridView();
    SDGV sdgv1 = sdg1.leafGridView();

    /*
     * Initialize subdomains. The right half of the domain is assigned to
     * subdomain 0, the left half to subdomain 1.
     */

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

    std::cout << "grid setup: " << timer.elapsed() << " sec" << std::endl;
    timer.reset();

    /*
     * Finite Element Map Setup
     * FEM0 will be for the right subdomain, FEM1 for the left
     */

    typedef ctype DF;
    typedef double RF;
    typedef double TReal;

    typedef Dune::PDELab::Q22DLocalFiniteElementMap<DF,RF> FEM0;
    typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,RF,dim> FEM1;

    typedef FEM0::Traits::LocalFiniteElementType::Traits::
      LocalBasisType::Traits::RangeFieldType R;

    FEM0 fem0;
    FEM1 fem1;

    /*
     * Constraints Setup
     * The empty constraints will be used for the underlying grid function spaces, as
     * we place the actual constraints on the subproblems.
     */

    typedef Dune::PDELab::NoConstraints NOCON;
    typedef Dune::PDELab::ConformingDirichletConstraints CON;

    CON con;

    /*
     * Grid Function Space Setup
     * Both grid function space are only defined on the part of the grid assigned
     * to the SubDomainGrid they are given at construction time
     */

    typedef Dune::PDELab::ISTLVectorBackend<> VBE;

    typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM0,NOCON,VBE> GFS0;
    typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM1,NOCON,VBE> GFS1;

    GFS0 gfs0(sdgv0,fem0);
    GFS1 gfs1(sdgv1,fem1);

    /*
     * The MultiDomainGridFunctionSpace needs to be parameterized on the underlying
     * MultiDomainGrid and the (basic PDELab) GridFunctionSpaces it should contain.
     * It does in general behave similar to a CompositeGridFunctionSpace.
     */

    typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,GFS0,GFS1> MultiGFS;

    MultiGFS multigfs(grid,gfs0,gfs1);

    std::cout << "function space setup: " << timer.elapsed() << " sec" << std::endl;
    timer.reset();

    /*
     * Parameter functions for the Poisson problem
     */

    typedef B<MDGV> BType;
    BType b(mdgv);

    typedef F<MDGV,R,double> FType;
    FType f(mdgv);

    typedef G<MDGV,R,double> GType;
    GType g(mdgv);

    typedef J<MDGV,R,double> JType;
    JType j(mdgv);

    /*
     * The local Poisson operators for the two subproblems.
     * In principle, these two are identical, but since we
     * are using different discretizations on the two subdomains,
     * they differ in their integration order.
     *
     * The temporal operator takes the integration order in its
     * constructor.
     */

    typedef InstationaryPoisson<FType,BType,JType,4> LOP0;
    LOP0 lop0(f,b,j);
    typedef InstationaryPoisson<FType,BType,JType,2> LOP1;
    LOP1 lop1(f,b,j);

    typedef SimpleTimeOperator TOP;
    TOP top0(4);
    TOP top1(2);

    /*
     * SubProblem Setup
     * The subproblems delegate the decision of whether or not they
     * apply to a given cell of the grid to a 'condition class'.
     * Currently there are two such classes: SubDomainEqualityCondition
     * and SubDomainSupersetCondition. They are simply passed a number
     * of subdomain identifiers at construction time.
     */

    typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> Condition;

    Condition c0(0);
    Condition c1(1);

    /*
     * The actual SubProblems are typed on the (ansatz and test) MultiDomainGridFunctionSpaces
     * that they are based on along with corresponding constraints assemblers, their local
     * operator(s), the condition class that defines the subdomain of the grid the problems apply
     * to and finally the components of the MultiDomainGridFunctionSpaces where the corresponding
     * functions are defined.
     * There are two ways of specifying those components:
     * - The TypeBased... variants can be passed the type of the component spaces. Obviously, this
     *   will only work if all component spaces have a distinct type.
     * - The basic variants (without TypeBased...) take the child indices of the components within
     *   the MultiDomainGridFunctionSpace. This always works, but may be prone to hard-to-find errors.
     * The ordering of the component function spaces is important, as it will be identical to the ordering
     * of the local function spaces in the local operator.
     */

    typedef Dune::PDELab::MultiDomain::TypeBasedInstationarySubProblem<TReal,MultiGFS,CON,MultiGFS,CON,LOP0,TOP,Condition,GFS0> SubProblem0;
    typedef Dune::PDELab::MultiDomain::TypeBasedInstationarySubProblem<TReal,MultiGFS,CON,MultiGFS,CON,LOP1,TOP,Condition,GFS1> SubProblem1;
    SubProblem0 sp0(con,con,lop0,top0,c0);
    SubProblem1 sp1(con,con,lop1,top1,c1);

    /*
     * The local operator coupling the two subproblems. This operator only defines a method
     * alphaCoupling() that will be called for all intersections on the interface between the
     * two subproblems.
     */

    typedef ContinuousValueContinuousFlowCoupling<TReal> COP;
    COP cop(4,atof(argv[2]));

    /*
     * The class defining the coupling of the two subproblems. It gets passed the two subproblems
     * and the local operator for calculating the coupling term. Note that the system will enforce
     * the subproblem ordering specified here: The local operator will always be called with the
     * first subproblem defined on the inside of the intersection and the second subproblem on its
     * outside.
     */

    typedef Dune::PDELab::MultiDomain::InstationaryCoupling<TReal,SubProblem1,SubProblem0,COP> Coupling;
    Coupling coupling(sp1,sp0,cop);

    std::cout << "subproblem / coupling setup: " << timer.elapsed() << " sec" << std::endl;
    timer.reset();


    /*
     * Constraints Evaluation
     * We have a Galerkin scheme, so it is sufficient to only evaluate
     * the constraints on the trial space and reuse them on the test space.
     * The extended constraints() function takes a boundary condition type
     * function for the MultiDomainGridFunctionSpace, the corresponding
     * MultiDomainGridFunctionSpace, the container in which to save the resulting
     * constraints, followed by pairs of boundary condition type functions and
     * their associated subproblems.
     */

    typedef GFS0::ConstraintsContainer<R>::Type C;
    C cg;

    Dune::PDELab::MultiDomain::trialSpaceConstraints(b,multigfs,cg,b,sp0,b,sp1);

    std::cout << "constraints evaluation: " << timer.elapsed() << " sec" << std::endl;
    timer.reset();

    /*
     * Interpolation of initial guess
     */

    typedef MultiGFS::VectorContainer<R>::Type V;
    V xold(multigfs);
    xold = 0.0;
    Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,xold,g,sp0,g,sp1);

    std::cout << "interpolation: " << timer.elapsed() << " sec" << std::endl;
    std::cout << xold.size() << " dof total, " << cg.size() << " dof constrained" << std::endl;
    timer.reset();

    /*
     * MultiDomainGridOperatorSpace Setup
     * Just like the basic InstationaryGridOperatorSpace, the
     * InstationaryMultiDomainGridOperatorSpace takes the type representing time values,
     * the type for storing residuals, the trial and test function spaces and
     * the matrix backend. After those basic parameters, you simply list all of the
     * subproblems and couplings to include in the calculations.
     */

    typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

    typedef Dune::PDELab::MultiDomain::InstationaryMultiDomainGridOperatorSpace<TReal,V,MultiGFS,MultiGFS,MBE,SubProblem0,SubProblem1,Coupling> MultiGOS;

    MultiGOS multigos(multigfs,multigfs,cg,cg,sp0,sp1,coupling);

    std::cout << "operator space setup: " << timer.elapsed() << " sec" << std::endl;
    timer.reset();

    /*
     * Time Stepping Method Setup
     * This is identical to standard PDELab.
     */

    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
    LS ls(5000,false);

    typedef Dune::PDELab::StationaryLinearProblemSolver<MultiGOS,LS,V> PDESOLVER;
    PDESOLVER pdesolver(multigos,ls,1e-10);

    Dune::PDELab::Alexander2Parameter<double> method;
    Dune::PDELab::OneStepMethod<double,MultiGOS,PDESOLVER,V,V> osm(method,multigos,pdesolver);
    osm.setVerbosityLevel(2);

    /*
     * Grid function sub spaces for VTK output
     * As for the subproblems, there is an enhanced version of the sub space
     * which allows specifying the child by type instead of by index (the same
     * restrictions apply).
     */

    typedef Dune::PDELab::MultiDomain::TypeBasedGridFunctionSubSpace<MultiGFS,GFS0> SGFS0;
    typedef Dune::PDELab::MultiDomain::TypeBasedGridFunctionSubSpace<MultiGFS,GFS1> SGFS1;
    SGFS0 sgfs0(multigfs);
    SGFS1 sgfs1(multigfs);

    /*
     * Output initial guess
     */

    Dune::PDELab::FilenameHelper fn0("instationaryexplicitlycoupledpoisson-right");
    {
      typedef Dune::PDELab::DiscreteGridFunction<SGFS0,V> U0DGF;
      U0DGF u0dgf(sgfs0,xold);
      Dune::SubsamplingVTKWriter<SDGV> vtkwriter(sdgv0,2);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"c0"));
      vtkwriter.write(fn0.getName(),Dune::VTKOptions::binaryappended);
      fn0.increment();
    }
    Dune::PDELab::FilenameHelper fn1("instationaryexplicitlycoupledpoisson-left");
    {
      typedef Dune::PDELab::DiscreteGridFunction<SGFS1,V> U1DGF;
      U1DGF u1dgf(sgfs1,xold);
      Dune::VTKWriter<SDGV> vtkwriter(sdgv1,Dune::VTKOptions::conforming);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"c1"));
      vtkwriter.write(fn1.getName(),Dune::VTKOptions::binaryappended);
      fn1.increment();
    }

    /*
     * Time stepping - also identical to standard PDELab version.
     */

    V xnew(xold);
    double time = 0;
    const double dt = atof(argv[3]);
    const double tend = atof(argv[4]);

    while (time < tend)
      {
        osm.apply(time,dt,xold,xnew);
        xold = xnew;
        time += dt;

        {
          typedef Dune::PDELab::DiscreteGridFunction<SGFS0,V> U0DGF;
          U0DGF u0dgf(sgfs0,xold);
          Dune::SubsamplingVTKWriter<SDGV> vtkwriter(sdgv0,2);
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"c0"));
          vtkwriter.write(fn0.getName(),Dune::VTKOptions::binaryappended);
          fn0.increment();
        }

        {
          typedef Dune::PDELab::DiscreteGridFunction<SGFS1,V> U1DGF;
          U1DGF u1dgf(sgfs1,xold);
          Dune::VTKWriter<SDGV> vtkwriter(sdgv1,Dune::VTKOptions::conforming);
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"c1"));
          vtkwriter.write(fn1.getName(),Dune::VTKOptions::binaryappended);
          fn1.increment();
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
