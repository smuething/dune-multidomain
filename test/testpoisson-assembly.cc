#include "config.h"

#include <dune/common/parametertreeparser.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include <dune/pdelab/localoperator/poisson.hh>
#include <dune/pdelab/localoperator/l2.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/multidomain/vtk.hh>

#include<typeinfo>

#include "functionmacros.hh"
#include "proportionalflowcoupling.hh"
#include "neumanndirichletcoupling.hh"

#include <dune/pdelab/common/benchmarkhelper.hh>


template<typename GV, typename RF>
class CouplingParameters
{

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1 : 0;
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

};


// source term
SIMPLE_ANALYTIC_FUNCTION(F,x,y)
{
  if (x[0]>0.25 && x[0]<0.375 && x[1]>0.25 && x[1]<0.375)
    y = 50.0;
  else {
    auto xm = x;
    xm = 0.5;
    xm -= x;
    if (xm.two_norm() < 0.2)
      y = 60 * std::exp(-10*x.two_norm2());
    else
      y = 0.0;
  }
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


  template<typename R>
  void setTime(const R& r)
  {}
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





template<typename GV, typename RF>
class ConvectionDiffusionModelProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1 : 0;
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& xl) const
  {
    auto x = e.geometry().global(xl);
    if (x[0]>0.25 && x[0]<0.375 && x[1]>0.25 && x[1]<0.375)
      return 50.0;
    else {
      auto xm = x;
      xm = 0.5;
      xm -= x;
      if (xm.two_norm() < 0.2)
        return 60 * std::exp(10*x.two_norm2());
    }
    return 0.0;
  }

  //! boundary condition type function
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    auto xg = is.geometry().global(x);

    if (!is.boundary())
      {
        return BCType::None;
      }

    if (xg[1]<1E-6 || xg[1]>1.0-1E-6)
      {
        return BCType::Neumann;
      }

    if (xg[0]>1.0-1E-6 && xg[1]>0.5+1E-6)
      {
        return BCType::Neumann;
      }

    return BCType::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType center(0.5);
    auto xg = e.geometry().global(x);
    center -= xg;
    return std::exp(-center.two_norm2());
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xl) const
  {
    auto x = is.geometry().global(xl);
    if (x[1]<1E-6 || x[1]>1.0-1E-6)
      {
        return 0;
      }
    if (x[0]>1.0-1E-6 && x[1]>0.5+1E-6)
      {
        return -5.0;
      }
    return 0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }
};



struct NoConstLeafOrderingParams
  : public Dune::PDELab::NoConstOrderingSize<true>
  , public Dune::PDELab::AllPartitionSelector
{};

using LeafTag = Dune::PDELab::LeafOrderingTag<NoConstLeafOrderingParams>;


template<template<class,class,class,int> class Preconditioner,
         template<class> class Solver>
class ISTL_SEQ_Subblock_Backend
  : public Dune::PDELab::SequentialNorm
  , public Dune::PDELab::LinearResultStorage
{
public:
  /*! \brief make a linear solver object

    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit ISTL_SEQ_Subblock_Backend(unsigned block, unsigned maxiter_=5000, int verbose_=1)
    : _block(block)
    , maxiter(maxiter_)
    , verbose(verbose_)
  {}

  unsigned block() const
  {
    return _block;
  }

  void setBlock(unsigned block)
  {
    _block = block;
  }

  /*! \brief solve the given linear system

    \param[in] A the given matrix
    \param[out] z the solution vector to be computed
    \param[in] r right hand side
    \param[in] reduction to be achieved
  */
  template<class M, class V, class W>
  void apply(M& A, V& z, W& r, typename W::ElementType reduction)
  {
    static_assert(std::is_same<V,W>::value,"V and W must be identical");

    typedef typename Dune::PDELab::istl::raw_type<M>::type ISTLMatrix;
    typedef typename Dune::PDELab::istl::raw_type<V>::type ISTLVector;

    typedef typename ISTLMatrix::block_type MatrixBlock;
    typedef typename ISTLVector::block_type VectorBlock;

    /*
    std::cout << "incoming residuals:" << std::endl
              << "block 0: " << Dune::PDELab::istl::raw(r)[0].two_norm() << std::endl
              << "block 1: " << Dune::PDELab::istl::raw(r)[1].two_norm() << std::endl;
    */

    Dune::MatrixAdapter<
      MatrixBlock,
      VectorBlock,
      VectorBlock
      > opa(Dune::PDELab::istl::raw(A)[_block][_block]);
    Preconditioner<
      MatrixBlock,
      VectorBlock,
      VectorBlock,
      1
      > prec(Dune::PDELab::istl::raw(A)[_block][_block], 3, 1.0);
    Solver<
      VectorBlock
      > solver(opa, prec, reduction, maxiter, verbose);
    Dune::InverseOperatorResult stat;
    solver.apply(Dune::PDELab::istl::raw(z)[_block], Dune::PDELab::istl::raw(r)[_block], stat);
    res.converged  = stat.converged;
    res.iterations = stat.iterations;
    res.elapsed    = stat.elapsed;
    res.reduction  = stat.reduction;
    res.conv_rate  = stat.conv_rate;
  }

private:
  unsigned _block;
  unsigned maxiter;
  int verbose;
};

template<
  typename GV,
  typename FEM,
  typename LOP,
  typename B,
  typename G
  >
void assemble(
  std::string name,
  std::size_t runs,
  std::string separator,
  const GV& gv,
  const FEM& fem,
  LOP& lop,
  const B& b,
  const G& g,
  std::size_t pattern_size)
{

  Dune::PDELab::BenchmarkHelper<> bh(name,runs,separator);

  for (std::size_t run = 0; run < runs; ++run)
    {

      bh.start_run(std::cout);

      typedef typename FEM::Traits::FiniteElementType::Traits::
        LocalBasisType::Traits::RangeFieldType R;

      typedef Dune::PDELab::ConformingDirichletConstraints CON;
      typedef Dune::PDELab::ISTLVectorBackend<> VBE;

      CON con;

      typedef Dune::PDELab::GridFunctionSpace<
        GV,
        FEM,
        CON,
        VBE,
        LeafTag
        > GFS;

      GFS gfs(gv,fem,con);
      gfs.name("u");

      typedef typename GFS::template ConstraintsContainer<R>::Type C;
      C cg;

      bh.start("constraints",std::cout);
      Dune::PDELab::constraints(b,gfs,cg);
      bh.end("constraints",std::cout);

      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
      MBE mbe(pattern_size);

      typedef Dune::PDELab::GridOperator<
        GFS,GFS,LOP,
        MBE,R,R,R,C,C
        > GridOperator;
      GridOperator grid_operator(gfs,cg,gfs,cg,lop,mbe);

      using V = typename GridOperator::Traits::Range;
      using M = typename GridOperator::Traits::Jacobian;

      bh.start("vector creation",std::cout);
      V u(gfs);
      bh.end("vector creation",std::cout);

      bh.start("matrix creation",std::cout);
      M m(grid_operator);
      bh.end("matrix creation",std::cout);

      V r(u);

      bh.start("residual",std::cout);
      grid_operator.residual(u,r);
      bh.end("residual",std::cout);

      bh.start("jacobian",std::cout);
      grid_operator.jacobian(u,m);
      bh.end("jacobian",std::cout);

      bh.end_run(std::cout);

      if (run == 0)
        {
          std::cout << std::endl << gfs.ordering().size() << " dof total, " << cg.size() << " dof constrained" << std::endl;
          std::cout << "max local size: " << gfs.ordering().maxLocalSize() << std::endl;
          std::cout << std::endl << m.patternStatistics() << std::endl;
        }

    }
  bh.print(std::cout);
}




template<
  typename Grid,
  typename GV,
  typename FEM,
  typename LOP,
  typename B,
  typename G
  >
void assemble_wrapped(
  std::string name,
  std::size_t runs,
  std::string separator,
  Grid& grid,
  const GV& gv,
  const FEM& fem,
  LOP& lop,
  const B& b,
  const G& g,
  std::size_t pattern_size)
{

  Dune::PDELab::BenchmarkHelper<> bh(name,runs,separator);

  for (std::size_t run = 0; run < runs; ++run)
    {

      bh.start_run(std::cout);

      typedef typename FEM::Traits::FiniteElementType::Traits::
        LocalBasisType::Traits::RangeFieldType R;

      typedef Dune::PDELab::ConformingDirichletConstraints CON;
      typedef Dune::PDELab::ISTLVectorBackend<> VBE;

      CON con;

      typedef Dune::PDELab::GridFunctionSpace<typename GV::Grid::SubDomainGrid::LeafGridView,FEM,CON,VBE> GFS;

      GFS gfs(gv.grid().subDomain(0).leafGridView(),fem,con);
      gfs.name("u");

      typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<
        Grid,
        VBE,
        Dune::PDELab::LexicographicOrderingTag,
        GFS
        > MultiGFS;

      MultiGFS multigfs(grid,gfs);

      typedef typename MultiGFS::template ConstraintsContainer<R>::Type C;
      C cg;

      auto condition = [](const typename GV::IndexSet::SubDomainSet&) { return true; };

      typedef Dune::PDELab::MultiDomain::SubProblem<
        MultiGFS,
        MultiGFS,
        LOP,
        decltype(condition),
        0
        > SubProblem;

      SubProblem sp(lop,condition);

      bh.start("constraints",std::cout);
      auto constraints = Dune::PDELab::MultiDomain::constraints<R>(
        multigfs,
        Dune::PDELab::MultiDomain::constrainSubProblem(
          sp,
          b)
        );

      constraints.assemble(cg);
      bh.end("constraints",std::cout);

      using MBE = Dune::PDELab::istl::BCRSMatrixBackend<>;
      MBE mbe(pattern_size);

      typedef Dune::PDELab::MultiDomain::GridOperator<
        MultiGFS,MultiGFS,
        MBE,R,R,R,C,C,
        SubProblem
        > GridOperator;
      GridOperator grid_operator(
        multigfs,multigfs,
        cg,cg,
        mbe,
        sp);

      using V = typename GridOperator::Traits::Range;
      using M = typename GridOperator::Traits::Jacobian;

      bh.start("vector creation",std::cout);
      V u(multigfs);
      bh.end("vector creation",std::cout);

      bh.start("matrix creation",std::cout);
      M m(grid_operator);
      bh.end("matrix creation",std::cout);

      V r(u);

      bh.start("residual",std::cout);
      grid_operator.residual(u,r);
      bh.end("residual",std::cout);

      bh.start("jacobian",std::cout);
      grid_operator.jacobian(u,m);
      bh.end("jacobian",std::cout);

      bh.end_run(std::cout);

      if (run == 0)
        {
          std::cout << std::endl << multigfs.ordering().size() << " dof total, " << cg.size() << " dof constrained" << std::endl;
          std::cout << "max local size: " << multigfs.ordering().maxLocalSize() << std::endl;
          std::cout << std::endl << m.patternStatistics() << std::endl;
        }

    }
  bh.print(std::cout);
}



template<
  typename Grid,
  typename MDGV,
  typename SDGV,
  typename FEM0, typename LOP0,
  typename FEM1, typename LOP1,
  typename B, typename G,
  typename MDVBE, typename MBE>
void assemble_multidomain(
  std::string name, std::size_t runs, std::string separator,
  Grid& grid, MDGV mdgv,
  SDGV sdgv0, const FEM0& fem0, LOP0& lop0,
  SDGV sdgv1, const FEM1& fem1, LOP1& lop1,
  const Dune::ParameterTree& parameters,
  const B& b, const G& g,
  MDVBE mdvbe, MBE mbe
  )
{

  Dune::PDELab::BenchmarkHelper<> bh(name,runs,separator);

  for (std::size_t run = 0; run < runs; ++run)
    {

      bh.start_run(std::cout);

      typedef typename FEM0::Traits::FiniteElementType::Traits::
        LocalBasisType::Traits::RangeFieldType R;

      typedef Dune::PDELab::ConformingDirichletConstraints CON;
      typedef Dune::PDELab::ISTLVectorBackend<> VBE;

      CON con;

      typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM0,CON,VBE> GFS0;

      typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM1,CON,VBE> GFS1;

      GFS0 gfs0(sdgv0,fem0,con);
      GFS1 gfs1(sdgv1,fem1,con);
      gfs0.name("u");
      gfs1.name("u");

      typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<
        typename MDGV::Grid,
        MDVBE,
        Dune::PDELab::LexicographicOrderingTag,
        GFS0,
        GFS1
        > MultiGFS;

      MultiGFS multigfs(grid,gfs0,gfs1);

      typedef typename MultiGFS::template ConstraintsContainer<R>::Type C;
      C cg;

      typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<typename MDGV::Grid> Condition;

      Condition c0({0});
      Condition c1({1});

      typedef Dune::PDELab::MultiDomain::SubProblem<
        MultiGFS,
        MultiGFS,
        LOP0,
        Condition,
        0
        > LeftSubProblem;

      LeftSubProblem left_sp(lop0,c0);

      typedef Dune::PDELab::MultiDomain::SubProblem<
        MultiGFS,
        MultiGFS,
        LOP1,
        Condition,
        1
        > RightSubProblem;

      RightSubProblem right_sp(lop1,c1);

      ContinuousValueContinuousFlowCoupling<R> proportionalFlowCoupling(4,1.0);

      typedef Dune::PDELab::MultiDomain::Coupling<
        LeftSubProblem,
        RightSubProblem,
        ContinuousValueContinuousFlowCoupling<R>
        > Coupling;
      Coupling coupling(left_sp,right_sp,proportionalFlowCoupling);

      bh.start("constraints",std::cout);
      auto constraints = Dune::PDELab::MultiDomain::constraints<R>(
        multigfs,
        Dune::PDELab::MultiDomain::constrainSubProblem(
          left_sp,
          b),
        Dune::PDELab::MultiDomain::constrainSubProblem(
          right_sp,
          b)
        );
      constraints.assemble(cg);
      bh.end("constraints",std::cout);

      typedef Dune::PDELab::MultiDomain::GridOperator<
        MultiGFS,MultiGFS,
        MBE,R,R,R,C,C,
        LeftSubProblem,
        RightSubProblem,
        Coupling
        > GridOperator;
      GridOperator grid_operator(
        multigfs,multigfs,
        cg,cg,
        mbe,
        left_sp,
        right_sp,
        coupling);

      using V = typename GridOperator::Traits::Range;
      using M = typename GridOperator::Traits::Jacobian;

      bh.start("vector creation",std::cout);
      V u(multigfs);
      bh.end("vector creation",std::cout);

      bh.start("matrix creation",std::cout);
      M m(grid_operator);
      bh.end("matrix creation",std::cout);

      V r(u);

      bh.start("residual",std::cout);
      grid_operator.residual(u,r);
      bh.end("residual",std::cout);

      bh.start("jacobian",std::cout);
      grid_operator.jacobian(u,m);
      bh.end("jacobian",std::cout);

      bh.end_run(std::cout);

      if (run == 0)
        {
          std::cout << std::endl << multigfs.ordering().size() << " dof total, " << cg.size() << " dof constrained" << std::endl;
          std::cout << "max local size: " << multigfs.ordering().maxLocalSize() << std::endl;
        }

    }
  bh.print(std::cout);
}


template<
  typename Grid,
  typename MDGV,
  typename SDGV,
  typename FEM0, typename LOP0,
  typename FEM1, typename LOP1,
  typename B, typename G>
void solve(
  Grid& grid, MDGV mdgv,
  SDGV sdgv0, const FEM0& fem0, const LOP0& lop0_,
  SDGV sdgv1, const FEM1& fem1, const LOP1& lop1_,
  const Dune::ParameterTree& parameters,
  const B& b, const G& g
  )
{
  auto lop0 = lop0_;
  auto lop1 = lop1_;


  std::size_t columns_left = parameters.get<int>("discretization.matrixcolumns.left");
  std::size_t columns_right = parameters.get<int>("discretization.matrixcolumns.right");
  std::size_t columns_multi = parameters.get<int>("discretization.matrixcolumns");
  std::size_t runs = parameters.get<int>("benchmark.runs");
  std::string separator = parameters.get("benchmark.separator","");

  if (parameters.get("benchmark.left-standard",true))
    assemble("FEM0 standard",runs,separator,mdgv,fem0,lop0,b,g,columns_left);

  if (parameters.get("benchmark.left-wrapped",true))
    assemble_wrapped("FEM0 wrapped",runs,separator,grid,mdgv,fem0,lop0,b,g,columns_left);

  if (parameters.get("benchmark.right-standard",true))
    assemble("FEM1 standard",runs,separator,mdgv,fem1,lop1,b,g,columns_right);

  if (parameters.get("benchmark.right-wrapped",true))
    assemble_wrapped("FEM1 wrapped",runs,separator,grid,mdgv,fem1,lop1,b,g,columns_right);

  if (parameters.get("benchmark.multidomain-flat",true))
    assemble_multidomain(
      "Multi domain flat",
      runs,
      separator,
      grid,mdgv,
      sdgv0,fem0,lop0,
      sdgv1,fem1,lop1,
      parameters,
      b,g,
      Dune::PDELab::ISTLVectorBackend<>(),
      Dune::PDELab::istl::BCRSMatrixBackend<>(columns_multi)
      );

  using BlockedVBE =
    Dune::PDELab::ISTLVectorBackend<
      Dune::PDELab::ISTLParameters::dynamic_blocking
    >;

  using EntriesPerRow = std::array<std::array<std::size_t,2>,2>;
  using BlockedMBE =  Dune::PDELab::istl::BCRSMatrixBackend<EntriesPerRow>;
  BlockedMBE blocked_mbe(
    {{
        {columns_left,  0},
        {0, columns_right}
      }});

  if (parameters.get("benchmark.multidomain-blocked",true))
    assemble_multidomain(
      "Multi domain blocked",
      runs,
      separator,
      grid,mdgv,
      sdgv0,fem0,lop0,
      sdgv1,fem1,lop1,
      parameters,
      b,g,
      BlockedVBE(),
      blocked_mbe
    );

}


template<typename FEM, typename Params>
Dune::PDELab::ConvectionDiffusionDG<
  Params,
  FEM
  >
dglop(const FEM& fem, Params& params)
{
  return Dune::PDELab::ConvectionDiffusionDG<
    Params,
    FEM
    >(
      params,
      Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,
      Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,
      2.0
      );
}


int main(int argc, char** argv) {

  try {

    Dune::MPIHelper::instance(argc,argv);

    if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <ini file>" << std::endl;
      return 1;
    }

    Dune::ParameterTree parameters;
    Dune::ParameterTreeParser::readINITree(argv[1],parameters);

    Dune::Timer timer;
    Dune::Timer totalTimer;
    timer.start();
    const int dim = 2;
    typedef Dune::YaspGrid<dim> BaseGrid;
    const Dune::FieldVector<double,dim> h(1.0);
    const Dune::array<int,dim> s = { {1,1} };
    const std::bitset<dim> p(false);
    BaseGrid baseGrid(h,s,p,0);
    baseGrid.globalRefine(parameters.get<int>("mesh.refine"));
    typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,4> > Grid;
    Grid grid(baseGrid,false);
    typedef Grid::SubDomainGrid SubDomainGrid;
    SubDomainGrid& sdg0 = grid.subDomain(0);
    SubDomainGrid& sdg1 = grid.subDomain(1);
    typedef Grid::LeafGridView MDGV;
    typedef SubDomainGrid::LeafGridView SDGV;
    MDGV mdgv = grid.leafGridView();
    SDGV sdgv0 = sdg0.leafGridView();
    SDGV sdgv1 = sdg1.leafGridView();
    sdg0.hostEntityPointer(*sdgv0.begin<0>());
    grid.startSubDomainMarking();
    for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
      {
        if (it->geometry().center()[0] > 0.5)
          grid.addToSubDomain(0,*it);
        else
          grid.addToSubDomain(0,*it);
      }
    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();

    std::cout << "grid setup: " << timer.elapsed() << " sec" << std::endl;
    timer.reset();

    typedef MDGV::Grid::ctype DF;

    typedef Dune::PDELab::QkLocalFiniteElementMap<SDGV,DF,double,1> FEM1;
    typedef Dune::PDELab::QkDGLocalFiniteElementMap<DF,double,1,dim> DGFEM1;
    typedef Dune::PDELab::QkDGLocalFiniteElementMap<DF,double,2,dim> DGFEM2;
    typedef Dune::PDELab::QkLocalFiniteElementMap<SDGV,DF,double,2> FEM2;

    typedef FEM1::Traits::FiniteElementType::Traits::
      LocalBasisType::Traits::RangeFieldType R;

    typedef DirichletBoundary BType;
    BType b;

    typedef F<MDGV,R> FType;
    FType f(mdgv);

    typedef G<MDGV,R> GType;
    GType g(mdgv);

    typedef J<MDGV,R> JType;
    JType j(mdgv);

    typedef Dune::PDELab::Poisson<FType,BType,JType> CGLOP;
    CGLOP cglop(f,b,j,6);

    using Params = ConvectionDiffusionModelProblem<
      MDGV,
      double
      >;
    Params params;

    solve(
      grid,mdgv,
      sdgv0,FEM1(sdgv0),cglop,
      //sdgv0,FEM2(sdgv0),cglop,
      //sdgv0,DGFEM1(),dglop(DGFEM1(),params),
      //sdgv0,DGFEM2(),dglop(DGFEM2(),params),
      sdgv1,FEM1(sdgv1),cglop,
      //sdgv1,FEM2(sdgv1),cglop,
      //sdgv1,DGFEM1(),dglop(DGFEM1(),params),
      //sdgv1,DGFEM2(),dglop(DGFEM2(),params),

      parameters,
      b,g
      );

    return 0;
  }
  catch (int &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }/*
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
    }*/

}
