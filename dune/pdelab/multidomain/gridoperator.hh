// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_GRIDOPERATOR_HH
#define DUNE_PDELAB_MULTIDOMAIN_GRIDOPERATOR_HH

#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridoperator/common/borderdofexchanger.hh>

#include <dune/pdelab/multidomain/localassembler.hh>
#include <dune/pdelab/multidomain/globalassembler.hh>
#include <dune/pdelab/multidomain/interpolate.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {

namespace {

  template<typename GFS, typename XG, typename Tuple, typename... InterpolationPairs>
  typename enable_if<(sizeof...(InterpolationPairs) == tuple_size<Tuple>::value)>::type
  interpolate_from_tuple(const GFS& gfs, XG& xg, const Tuple& tuple, const InterpolationPairs&... pairs)
  {
    Dune::PDELab::MultiDomain::interpolateOnTrialSpace(gfs,xg,pairs...);
  }

  template<typename GFS, typename XG, typename Tuple, typename... InterpolationPairs>
  typename enable_if<(sizeof...(InterpolationPairs) < tuple_size<Tuple>::value)>::type
  interpolate_from_tuple(const GFS& gfs, XG& xg, const Tuple& tuple, const InterpolationPairs&... pairs)
  {
    interpolate_from_tuple(gfs,xg,tuple,pairs...,get<sizeof...(InterpolationPairs)>(tuple),get<sizeof...(InterpolationPairs)+1>(tuple));
  }

} // anonymous namespace


template<typename GFSU,
         typename GFSV,
         typename MB,
         typename DF,
         typename RF,
         typename JF,
         typename CU,
         typename CV,
         typename... AssemblyParticipants>
class GridOperator
{

public:

  typedef GridOperatorTraits<
    GFSU,
    GFSV,
    MB,
    DF,
    RF,
    JF,
    CU,
    CV,
    GlobalAssembler<GFSU,GFSV>,
    LocalAssembler<GridOperator,AssemblyParticipants...>
    > Traits;

  static const bool nonoverlapping_mode = false;

  typedef typename Dune::conditional<
    nonoverlapping_mode,
    NonOverlappingBorderDOFExchanger<GridOperator>,
    OverlappingBorderDOFExchanger<GridOperator>
    >::type BorderDOFExchanger;


  template<typename P>
  void fill_pattern(P& globalpattern) const
  {
    _assembler.assemble(_localAssembler.localPatternAssemblerEngine(globalpattern));
  }

  void residual(const typename Traits::Domain& x, typename Traits::Range& r) const
  {
    _assembler.assemble(_localAssembler.localResidualAssemblerEngine(r,x));
  }

  void jacobian(const typename Traits::Domain& x, typename Traits::Jacobian& a) const
  {
    _assembler.assemble(_localAssembler.localJacobianAssemblerEngine(a,x));
  }

  typename Traits::Assembler& assembler()
  {
    return _assembler;
  }

  const typename Traits::Assembler& assembler() const
  {
    return _assembler;
  }

  typename Traits::LocalAssembler& localAssembler()
  {
    return _localAssembler;
  }

  const typename Traits::LocalAssembler& localAssembler() const
  {
    return _localAssembler;
  }

  const GFSU& trialGridFunctionSpace() const
  {
    return assembler().trialGridFunctionSpace();
  }

  const GFSV& testGridFunctionSpace() const
  {
    return assembler().testGridFunctionSpace();
  }

  typename GFSU::Traits::SizeType globalSizeU () const
  {
    return assembler().trialGridFunctionSpace().size();
  }

  typename GFSV::Traits::SizeType globalSizeV () const
  {
    return assembler().testGridFunctionSpace().size();
  }

  template<typename F>
  void interpolate(const typename Traits::Domain& xold,
                   const F& f,
                   typename Traits::Domain& xnew)
  {
    interpolate_from_tuple(_assembler.trialGridFunctionSpace(),xnew,f.base());

    // Copy non-constrained dofs from old time step
    Dune::PDELab::copy_nonconstrained_dofs(_localAssembler.trialConstraints(),xold,xnew);

    // Make solution consistent
    CopyDataHandle<typename Traits::TrialGridFunctionSpace,typename Traits::Range> cdh(trialGridFunctionSpace(),xnew);
    if (trialGridFunctionSpace().gridView().comm().size() > 1)
      trialGridFunctionSpace().gridView().communicate(cdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

  GridOperator(const GFSU& gfsu,
               const GFSV& gfsv,
               const CU& cu,
               const CV& cv,
               const MB& mb,
               AssemblyParticipants&... participants)
    : _assembler(gfsu,gfsv)
    , _localAssembler(cu,cv,participants...)
    , dof_exchanger(make_shared<BorderDOFExchanger>(*this))
  {}

  GridOperator(const GFSU& gfsu,
               const GFSV& gfsv,
               const CU& cu,
               const CV& cv,
               AssemblyParticipants&... participants)
    : _assembler(gfsu,gfsv)
    , _localAssembler(cu,cv,participants...)
    , dof_exchanger(make_shared<BorderDOFExchanger>(*this))
  {}

  template<typename Tuple, std::size_t k>
  static typename enable_if<(k == tuple_size<Tuple>::value - 1)>::type
  shareData(const Tuple& gridOperators)
  {
  }

  template<typename Tuple, std::size_t k>
  static typename enable_if<(k < tuple_size<Tuple>::value - 1)>::type
  shareData(const Tuple& gridOperators)
  {
    get<k+1>(gridOperators).localAssembler().shareData(get<k>(gridOperators).localAssembler());
    shareData<Tuple,k+1>(gridOperators);
  }

  template<typename... GridOperators>
  static void setupGridOperators(const tuple<GridOperators&...>& gridOperators)
  {
    shareData<tuple<GridOperators&...>,0>(gridOperators);
  }

  void make_consistent(typename Traits::Jacobian& a) const {
    dof_exchanger->accumulateBorderEntries(*this,a);
  }

 private:

  mutable typename Traits::Assembler _assembler;
  mutable typename Traits::LocalAssembler _localAssembler;
  shared_ptr<BorderDOFExchanger> dof_exchanger;

};


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_GRIDOPERATOR_HH
