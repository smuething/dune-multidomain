#include "config.h"

#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/istlhelpers.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinggridfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinglocalfunctionspace.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/multidomain/couplingutilities.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>

#include <dune/pdelab/multidomain/vtk.hh>


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



class MortarPoissonCoupling :
  public Dune::PDELab::MultiDomain::CouplingOperatorDefaultFlags,
  public Dune::PDELab::MultiDomain::NumericalJacobianEnrichedCouplingToFirstSubProblem<MortarPoissonCoupling>,
  public Dune::PDELab::MultiDomain::NumericalJacobianEnrichedCouplingToSecondSubProblem<MortarPoissonCoupling>,
  public Dune::PDELab::MultiDomain::FullEnrichedCouplingFirstSubProblemPattern,
  public Dune::PDELab::MultiDomain::FullEnrichedCouplingSecondSubProblemPattern,
  public Dune::PDELab::MultiDomain::FullEnrichedCouplingPattern,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{

public:

  static const bool doAlphaEnrichedCouplingToSubProblems = true;
  static const bool doPatternEnrichedCouplingToSubProblems = true;
  static const bool doPatternEnrichedCoupling = true;

  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename LFSU_C, typename LFSV_C, typename R>
  void alpha_enriched_coupling_first
  ( const IG& ig,
    const LFSU& lfsu, const X& x, const LFSV& lfsv,
    const LFSU_C& lfsu_c, const X& x_c, const LFSV_C& lfsv_c,
    R& r, R& r_c) const
  {
    this->template alpha_generic<1>(ig,lfsu,x,lfsv,lfsu_c,x_c,lfsv_c,r,r_c);
  }

  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename LFSU_C, typename LFSV_C, typename R>
  void alpha_enriched_coupling_second
  ( const IG& ig,
    const LFSU& lfsu, const X& x, const LFSV& lfsv,
    const LFSU_C& lfsu_c, const X& x_c, const LFSV_C& lfsv_c,
    R& r, R& r_c) const
  {
    this->template alpha_generic<-1>(ig,lfsu,x,lfsv,lfsu_c,x_c,lfsv_c,r,r_c);
  }


  template<int sign,
           typename IG, typename LFSU, typename X, typename LFSV,
           typename LFSU_C, typename LFSV_C, typename R>
  void alpha_generic
  ( const IG& ig,
    const LFSU& lfsu, const X& x, const LFSV& lfsv,
    const LFSU_C& lfsu_c, const X& x_c, const LFSV_C& lfsv_c,
    R& r, R& r_c) const
  {
    const int dimWorld = IG::dimensionworld;
    const int dimIF = IG::dimension - 1;

    typedef Dune::FiniteElementInterfaceSwitch<
      typename LFSU::Traits::FiniteElementType
      > LFSU_FESwitch;
    typedef Dune::BasisInterfaceSwitch<
      typename LFSU_FESwitch::Basis
      > LFSU_BasisSwitch;

    typedef typename LFSU_BasisSwitch::DomainField DF;
    typedef typename LFSU_BasisSwitch::Range Range;
    typedef typename LFSU::Traits::SizeType size_type;

    const size_type lfsu_size(lfsu.size());

    typedef Dune::FiniteElementInterfaceSwitch<
      typename LFSU_C::Traits::FiniteElementType
      > LFSU_C_FESwitch;

    const size_type lfsu_c_size(lfsu_c.size());

    typedef typename IG::Geometry::GlobalCoordinate GC;

    const Dune::GeometryType gt = ig.geometry().type();

    const double gamma = gamma_0 / ig.geometry().volume(); // TODO: replace with cell diameter

    const size_type qorder = 2 * std::max(lfsu.finiteElement().localBasis().order(),
                                          lfsu_c.finiteElement().localBasis().order()
                                          );

    const Dune::QuadratureRule<DF,dimIF>& rule = Dune::QuadratureRules<DF,dimIF>::rule(gt,qorder);
    typedef typename Dune::QuadratureRule<DF,dimIF>::const_iterator RuleIterator;

    const RuleIterator qend = rule.end();

    for (RuleIterator qit = rule.begin(); qit != qend; ++qit)
      {
        const GC element_pos = (sign > 0 ? ig.geometryInInside() : ig.geometryInOutside()).global(qit->position());
        GC normal = ig.centerUnitOuterNormal();

        std::vector<Range> v(lfsu_size);
        LFSU_FESwitch::basis(lfsu.finiteElement()).evaluateFunction(element_pos,v);

        std::vector<Range> mu(lfsu_c_size);
        LFSU_C_FESwitch::basis(lfsu_c.finiteElement()).evaluateFunction(qit->position(),mu);

        std::vector<Dune::FieldMatrix<DF,1,dimWorld> > gradv(lfsu_size);
        LFSU_BasisSwitch::gradient(LFSU_FESwitch::basis(lfsu.finiteElement()),
                                   (sign > 0 ? ig.insideElement().geometry() : ig.outsideElement().geometry()),
                                   element_pos,
                                   gradv
                                   );

        Range u(0.0);
        GC gradu(0.0);
        for (size_type i = 0; i < lfsu_size; ++i)
          {
            u += x(lfsu,i) * v[i];
            gradu.axpy(x(lfsu,i),gradv[i][0]);
          }

        Range lambda(0.0);
        for (size_type i = 0; i < lfsu_c_size; ++i)
          lambda += x_c(lfsu_c,i) * mu[i];

        const double factor = qit->weight() * ig.geometry().integrationElement(qit->position());

        for(size_t i = 0; i < lfsu_size; ++i)
          {
            r.accumulate(lfsv,i,(-(normal*gradu)*v[i]-(u-lambda)*(normal*gradv[i][0])+2*gamma*(u-lambda)*v[i])*factor);
          }

        for(size_t i = 0; i < lfsu_c_size; ++i)
          {
            r_c.accumulate(lfsv_c,i,((normal*gradu)-2*gamma*(u-lambda))*mu[i]*factor);
          }
      }
  }

  MortarPoissonCoupling(double gamma_0_)
    : gamma_0(gamma_0_)
  {}

  const double gamma_0;

};


template<typename SubDomain>
struct subdomain_predicate
{

  template<typename LFS>
  bool operator()(const LFS& lfs) const
  {
    return lfs.gridFunctionSpace().gridView().grid().domain() == _subDomain;
  }

  subdomain_predicate(const SubDomain& subDomain)
    : _subDomain(subDomain)
  {}

private:
  SubDomain _subDomain;

};

template<int dim>
struct PseudoGV
{
  static const int dimension = dim;
};

int main(int argc, char** argv) {

  try {

    Dune::MPIHelper::instance(argc,argv);

    if (argc != 3)
      std::cerr << "Usage: " << argv[0] << " <refinement> <scaling factor for coupling> " << std::endl
                << "Using default values refinement=5 scaling=1.0" << std::endl;
    const int refinement = argc == 3 ? atoi(argv[1]) : 5;
    const double scaling = argc == 3 ? atof(argv[2]) : 1.0;

    const int dim = 2;

    typedef Dune::YaspGrid<dim> BaseGrid;
    const Dune::FieldVector<double,dim> h(1.0);
    const Dune::array<int,dim> s = { {1,1} };
    const std::bitset<dim> p(false);
    BaseGrid baseGrid(h,s,p,0);
    baseGrid.globalRefine(refinement);

#ifdef UGGRID
    typedef Dune::UGGrid<dim> BaseGrid;
    BaseGrid baseGrid(500);
    std::vector<int> boundaryIndexToPhysicalGroup, elementIndexToPhysicalGroup;
    Dune::GmshReader<BaseGrid> gmshreader;
    gmshreader.read(baseGrid,"gmshtest.msh",boundaryIndexToPhysicalGroup,elementIndexToPhysicalGroup,true,false);
#endif

    typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,4> > Grid;
    Grid grid(baseGrid,false);
    typedef Grid::SubDomainGrid SubDomainGrid;
    SubDomainGrid& sdg0 = grid.subDomain(0);
    SubDomainGrid& sdg1 = grid.subDomain(1);

    typedef Grid::ctype ctype;
    typedef Grid::LeafGridView MDGV;
    typedef SubDomainGrid::LeafGridView SDGV;
    MDGV mdgv = grid.leafGridView();
    SDGV sdgv0 = sdg0.leafGridView();
    SDGV sdgv1 = sdg1.leafGridView();
    grid.startSubDomainMarking();
    for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
      {
#ifdef UGGRID
        grid.addToSubDomain(elementIndexToPhysicalGroup[mdgv.indexSet().index(*it)],*it);
#else
        Dune::FieldVector<ctype,dim> center = it->geometry().center();
        if (center[0] > 0.5)
          grid.addToSubDomain(1,*it);
        if (center[0] < 0.5)
          grid.addToSubDomain(0,*it);
#endif
      }
    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();

    typedef MDGV::Grid::ctype DF;

    // we need something to feed into the QkFEM, only needs to return a correct dimension
    using CouplingGV = PseudoGV<dim-1>;

    typedef Dune::PDELab::QkLocalFiniteElementMap<MDGV,DF,double,1> FEM;
    typedef Dune::PDELab::QkLocalFiniteElementMap<CouplingGV,DF,double,1> COUPLINGFEM;

    typedef FEM::Traits::FiniteElementType::Traits::
      LocalBasisType::Traits::RangeFieldType R;

    FEM fem(mdgv);
    COUPLINGFEM couplingfem({});

    typedef Dune::PDELab::NoConstraints NOCON;
    typedef Dune::PDELab::ConformingDirichletConstraints CON;

    typedef Dune::PDELab::ISTLVectorBackend<> VBE;

    CON con;

    typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM,CON,VBE> GFS0;
    typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM,CON,VBE> GFS1;

    GFS0 gfs0(sdgv0,fem,con);
    GFS1 gfs1(sdgv1,fem,con);
    gfs0.name("u");
    gfs1.name("u");

    typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> EC;
    EC c0({0});
    EC c1({1});

    typedef Dune::PDELab::MultiDomain::SubProblemSubProblemInterface<MDGV,EC,EC> Pred;
    Pred pred(mdgv,c0,c1);

    typedef Dune::PDELab::MultiDomain::CouplingGridFunctionSpace<MDGV,COUPLINGFEM,Pred,NOCON,VBE> CouplingGFS;
    CouplingGFS couplinggfs(mdgv,couplingfem,pred);

    typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<
      Grid,
      VBE,
      Dune::PDELab::LexicographicOrderingTag,
      GFS0,
      GFS1,
      CouplingGFS
      > MultiGFS;
    MultiGFS multigfs(grid,VBE(),gfs0,gfs1,couplinggfs);

    std::cout << multigfs.ordering().size() << std::endl;

    typedef DirichletBoundary BType;
    BType b;

    typedef F<MDGV,R> FType;
    FType f(mdgv);

    typedef G<MDGV,R> GType;
    GType g(mdgv);

    typedef J<MDGV,R> JType;
    JType j(mdgv);

    typedef Dune::PDELab::Poisson<FType,BType,JType> LOP;
    LOP lop(f,b,j,4);

    typedef MortarPoissonCoupling CouplingLOP;
    CouplingLOP coupling_lop(scaling);

    typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,EC,0> SubProblem0;
    typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,EC,1> SubProblem1;
    SubProblem0 subproblem0(lop,c0);
    SubProblem1 subproblem1(lop,c1);

    typedef Dune::PDELab::MultiDomain::EnrichedCoupling<SubProblem0,SubProblem1,CouplingLOP,2> Coupling;
    Coupling coupling(subproblem0,subproblem1,coupling_lop);

    typedef MultiGFS::ConstraintsContainer<R>::Type C;
    C cg;

    typedef Dune::PDELab::ISTLMatrixBackend MBE;

    typedef Dune::PDELab::MultiDomain::GridOperator<
      MultiGFS,MultiGFS,
      MBE,R,R,R,C,C,
      SubProblem0,
      SubProblem1,
      Coupling> GridOperator;

    typedef GridOperator::Traits::Domain V;
    V u(multigfs,0.0);

    Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,u,g,subproblem0,g,subproblem1);

    Dune::PDELab::MultiDomain::constraints<R>(
      multigfs,
      Dune::PDELab::MultiDomain::constrainSubProblem(subproblem0,b),
      Dune::PDELab::MultiDomain::constrainSubProblem(subproblem1,b)
    ).assemble(cg);

    GridOperator gridoperator(multigfs,multigfs,
                              cg,cg,
                              MBE(),
                              subproblem0,
                              subproblem1,
                              coupling);

    GridOperator::Traits::Range r(multigfs,0.0);

    gridoperator.residual(u,r);

    GridOperator::Traits::Jacobian J(gridoperator);

    gridoperator.jacobian(u,J);

    Dune::writeMatrixToMatlab(J.base(),"jacobian.mat");

    // <<<5>>> Select a linear solver backend
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
    LS ls(5000,false);

    // <<<6>>> Solver for linear problem per stage
    typedef Dune::PDELab::StationaryLinearProblemSolver<GridOperator,LS,V> PDESOLVER;
    PDESOLVER pdesolver(gridoperator,ls,1e-10);

    pdesolver.apply(u);

    {
      Dune::VTKWriter<SDGV> vtkwriter(sdgv0,Dune::VTK::conforming);
      Dune::PDELab::MultiDomain::addSolutionToVTKWriter(
        vtkwriter,
        multigfs,
        u,
        subdomain_predicate<Grid::SubDomainIndex>(0)
      );
      vtkwriter.write("testmortarpoisson-left",Dune::VTK::ascii);
    }

    {
      Dune::VTKWriter<SDGV> vtkwriter(sdgv1,Dune::VTK::conforming);
      Dune::PDELab::MultiDomain::addSolutionToVTKWriter(
        vtkwriter,
        multigfs,
        u,
        subdomain_predicate<Grid::SubDomainIndex>(1)
      );
      vtkwriter.write("testmortarpoisson-right",Dune::VTK::ascii);
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
