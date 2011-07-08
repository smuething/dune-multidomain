#include "config.h"

#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinggridfunctionspace.hh>
#include <dune/pdelab/multidomain/couplinglocalfunctionspace.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/multidomain/couplingutilities.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>


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
  public Dune::PDELab::MultiDomain::NumericalJacobianEnrichedCoupling<MortarPoissonCoupling>,
  public Dune::PDELab::MultiDomain::NumericalJacobianApplyCoupling<MortarPoissonCoupling>,
  public Dune::PDELab::MultiDomain::FullEnrichedCouplingFirstPattern,
  public Dune::PDELab::MultiDomain::FullEnrichedCouplingSecondPattern,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{

public:

  static const bool doAlphaEnrichedCoupling = true;
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
    const int dim = IG::dimension;
    const int dimWorld = IG::dimensionworld;
    const int dimIF = IG::dimension - 1;

    typedef Dune::FiniteElementInterfaceSwitch<
      typename LFSU::Traits::FiniteElementType
      > LFSU_FESwitch;
    typedef Dune::BasisInterfaceSwitch<
      typename LFSU_FESwitch::Basis
      > LFSU_BasisSwitch;

    typedef typename LFSU_BasisSwitch::DomainField DF;
    typedef typename LFSU_BasisSwitch::DomainField RF;
    typedef typename LFSU_BasisSwitch::Range R;
    typedef typename LFSU::Traits::SizeType size_type;

    const size_type lfsu_size(lfsu.size());

    typedef Dune::FiniteElementInterfaceSwitch<
      typename LFSU_C::Traits::FiniteElementType
      > LFSU_C_FESwitch;
    typedef Dune::BasisInterfaceSwitch<
      typename LFSU_C_FESwitch::Basis
      > LFSU_C_BasisSwitch;

    const size_type lfsu_c_size(lfsu_c.size());

    typedef typename IG::Geometry::LocalCoordinate LC;
    typedef typename IG::Geometry::GlobalCoordinate GC;

    const Dune::GeometryType gt = ig.geometry().type();

    const double gamma = gamma_0 / ig.geometry().volume(); // TODO: replace with cell diameter

    const size_type qorder = 2 * std::max(lfsu.finiteElement().localBasis().order(),
                                          lfsu_c.finiteElement().localBasis().order()
                                          );

    const Dune::QuadratureRule<RF,dimIF>& rule = Dune::QuadratureRules<RF,dimIF>::rule(gt,qorder);
    typedef typename Dune::QuadratureRule<RF,dimIF>::const_iterator RuleIterator;

    const RuleIterator qend = rule.end();

    for (RuleIterator qit = rule.begin(); qit != qend; ++qit)
      {
        const GC element_pos = (sign > 0 ? ig.geometryInInside() : ig.geometryInOutside()).global(qit->position());
        GC normal = ig.centerUnitOuterNormal();

        std::vector<R> v(lfsu_size);
        LFSU_FESwitch::basis(lfsu.finiteElement()).evaluateFunction(element_pos,v);

        std::vector<R> mu(lfsu_c_size);
        LFSU_C_FESwitch::basis(lfsu_c.finiteElement()).evaluateFunction(qit->position(),mu);

        std::vector<Dune::FieldMatrix<RF,1,dimWorld> > gradv(lfsu_size);
        LFSU_BasisSwitch::gradient(LFSU_FESwitch::basis(lfsu.finiteElement()),
                                   (sign > 0 ? ig.insideElement().geometry() : ig.outsideElement().geometry()),
                                   element_pos,
                                   gradv
                                   );

        R u(0.0);
        GC gradu(0.0);
        for (size_type i = 0; i < lfsu_size; ++i)
          {
            u += x(lfsu,i) * v[i];
            gradu.axpy(x(lfsu,i),gradv[i][0]);
          }

        R lambda(0.0);
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


int main(int argc, char** argv) {

  try {

  Dune::MPIHelper::instance(argc,argv);

  const int dim = 2;

  typedef Dune::YaspGrid<dim> BaseGrid;
  const Dune::FieldVector<int,dim> s(1);
  const Dune::FieldVector<double,dim> h(1.0);
  const Dune::FieldVector<bool,dim> p(false);
  BaseGrid baseGrid(h,s,p,0);
  baseGrid.globalRefine(2);

#ifdef UGGRID
  typedef Dune::UGGrid<dim> BaseGrid;
  BaseGrid baseGrid(500);
  std::vector<int> boundaryIndexToPhysicalGroup, elementIndexToPhysicalGroup;
  Dune::GmshReader<BaseGrid> gmshreader;
  gmshreader.read(baseGrid,"gmshtest.msh",boundaryIndexToPhysicalGroup,elementIndexToPhysicalGroup,true,false);
#endif

  typedef BaseGrid::LeafGridView GV;

  GV gv = baseGrid.leafView();
  const GV::IndexSet& is = gv.indexSet();

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

  typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM;
  typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim-1> COUPLINGFEM;

  typedef FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  FEM fem;
  COUPLINGFEM couplingfem;

  typedef Dune::PDELab::NoConstraints NOCON;
  typedef Dune::PDELab::ConformingDirichletConstraints CON;

  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;

  NOCON nocon;
  CON con;

  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM,CON,VBE> GFS0;
  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM,CON,VBE> GFS1;

  GFS0 gfs0(sdgv0,fem,con);
  GFS1 gfs1(sdgv1,fem,con);

  typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<Grid> EC;
  EC c0(0);
  EC c1(1);

  typedef Dune::PDELab::MultiDomain::SubProblemSubProblemInterface<MDGV,EC,EC> Pred;
  Pred pred(mdgv,c0,c1);

  typedef Dune::PDELab::MultiDomain::CouplingGridFunctionSpace<MDGV,COUPLINGFEM,Pred,NOCON,VBE> CouplingGFS;
  CouplingGFS couplinggfs(mdgv,couplingfem,pred);

  typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,GFS0,GFS1,CouplingGFS> MultiGFS;
  MultiGFS multigfs(grid,gfs0,gfs1,couplinggfs);

  typedef DirichletBoundary BType;
  BType b;

  typedef F<MDGV,R> FType;
  FType f(mdgv);

  typedef G<MDGV,R> GType;
  GType g(mdgv);

  typedef J<MDGV,R> JType;
  JType j(mdgv);

  typedef Dune::PDELab::Poisson<FType,BType,JType,4> LOP;
  LOP lop(f,b,j);

  typedef MortarPoissonCoupling CouplingLOP;
  CouplingLOP coupling_lop(1.0);

  typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,EC,0> SubProblem0;
  typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,EC,1> SubProblem1;
  SubProblem0 subproblem0(lop,c0);
  SubProblem1 subproblem1(lop,c1);

  typedef Dune::PDELab::MultiDomain::EnrichedCoupling<SubProblem0,SubProblem1,CouplingLOP,2> Coupling;
  Coupling coupling(subproblem0,subproblem1,coupling_lop);

  typedef MultiGFS::ConstraintsContainer<R>::Type C;
  C cg;

  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

  typedef Dune::PDELab::MultiDomain::GridOperator<
    MultiGFS,MultiGFS,
    MBE,R,R,R,C,C,
    SubProblem0,
    SubProblem1,
    Coupling> GridOperator;

  typedef GridOperator::Traits::Domain V;
  V u(multigfs,0.0);

  Dune::PDELab::MultiDomain::interpolateOnTrialSpace(multigfs,u,g,subproblem0,g,subproblem1);

  Dune::PDELab::MultiDomain::constraints<R>(multigfs,
                                            Dune::PDELab::MultiDomain::constrainSubProblem(subproblem0,b),
                                            Dune::PDELab::MultiDomain::constrainSubProblem(subproblem1,b)
                                            ).assemble(cg);

  GridOperator gridoperator(multigfs,multigfs,
                            cg,cg,
                            subproblem0,
                            subproblem1,
                            coupling);

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,false);

  // <<<6>>> Solver for linear problem per stage
  typedef Dune::PDELab::StationaryLinearProblemSolver<GridOperator,LS,V> PDESOLVER;
  PDESOLVER pdesolver(gridoperator,ls,1e-10);

  pdesolver.apply(u);

    typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,0> SGFS0;
    typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,1> SGFS1;
    SGFS0 sgfs0(multigfs);
    SGFS1 sgfs1(multigfs);

    {
      typedef Dune::PDELab::DiscreteGridFunction<SGFS0,V> DGF0;
      DGF0 dgf0(sgfs0,u);
      Dune::VTKWriter<SDGV> vtkwriter(sdgv0,Dune::VTKOptions::conforming);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF0>(dgf0,"u"));
      vtkwriter.write("testmortarpoisson-left",Dune::VTKOptions::ascii);
    }

    {
      typedef Dune::PDELab::DiscreteGridFunction<SGFS1,V> DGF1;
      DGF1 dgf1(sgfs1,u);
      Dune::VTKWriter<SDGV> vtkwriter(sdgv1,Dune::VTKOptions::conforming);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF1>(dgf1,"u"));
      vtkwriter.write("testmortarpoisson-right",Dune::VTKOptions::ascii);
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
