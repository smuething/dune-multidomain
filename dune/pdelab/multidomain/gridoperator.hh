// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_GRIDOPERATOR_HH
#define DUNE_PDELAB_MULTIDOMAIN_GRIDOPERATOR_HH

#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>

#include <dune/pdelab/multidomain/multidomaingridoperatorspaceutilities.hh>
#include <dune/pdelab/multidomain/localassembler.hh>
#include <dune/pdelab/multidomain/globalassembler.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


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
    GlobalAssembler<GridOperator>,
    LocalAssembler<GridOperator,AssemblyParticipants...>
    > Traits;


  template<typename P>
  void fill_pattern(P& globalpattern) const
  {
    _assembler.assemble(_localAssembler.patternAssemblerEngine(globalpattern));
  }

  void residual(const typename Traits::Domain& x, typename Traits::Range& r) const
  {
    _assembler.assemble(_localAssembler.residualAssemblerEngine(x,r));
  }

  void jacobian(const typename Traits::Domain& x, typename Traits::Jacobian& a) const
  {
    _assembler.assemble(_localAssembler.jacobianAssemblerEngine(x,a));
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
    return assembler().gfsu();
  }

  const GFSV& testGridFunctionSpace() const
  {
    return assembler().gfsv();
  }

  typename GFSU::Traits::SizeType globalSizeU () const
  {
    return assembler().gfsu().size();
  }

  typename GFSV::Traits::SizeType globalSizeV () const
  {
    return assembler().gfsu().size();
  }

  template<typename F>
  void interpolate(const typename Traits::Domain& xold,
                   const F& f,
                   typename Traits::Domain& xnew)
  {
  }

  GridOperator(const GFSU& gfsu,
               const GFSV& gfsv,
               const CU& cu,
               const CV& cv,
               AssemblyParticipants&... participants)
    : _assembler(gfsu,gfsv)
    , _localAssembler(cu,cv,participants...)
  {}

 private:

  mutable typename Traits::Assembler _assembler;
  mutable typename Traits::LocalAssembler _localAssembler;


};


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_GRIDOPERATOR_HH
