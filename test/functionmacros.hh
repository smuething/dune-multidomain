#ifndef DUNE_MULTIDOMAIN_FUNCTIONMACROS_HH
#define DUNE_MULTIDOMAIN_FUNCTIONMACROS_HH

#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

#define SIMPLE_ANALYTIC_FUNCTION(_NAME_,_POINT_,_VALUE_)  \
template<typename GV, typename RF> \
class _NAME_ \
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>, \
                                                  _NAME_<GV,RF> > \
{ \
public: \
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits; \
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,_NAME_<GV,RF> > BaseT; \
\
  _NAME_ (const GV& gv) : BaseT(gv) {} \
  inline void evaluateGlobal (const typename Traits::DomainType& _POINT_, \
                              typename Traits::RangeType& _VALUE_) const

#define END_SIMPLE_ANALYTIC_FUNCTION };

#define SIMPLE_BOUNDARYTYPE_FUNCTION(_NAME_,_INTERSECTION_,_POINT_,_VALUE_)    \
template<typename GV> \
class _NAME_ \
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab:: \
                                                  BoundaryTypeGridFunctionTraits<GV>, \
                                                  _NAME_<GV> > \
{ \
  const GV& gv; \
\
public: \
  typedef Dune::PDELab::BoundaryTypeGridFunctionTraits<GV> Traits; \
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,_NAME_<GV> > BaseT; \
\
  _NAME_ (const GV& gv_) : gv(gv_) {} \
\
  template<typename I>\
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& _INTERSECTION_, \
                        const typename Traits::DomainType& _POINT_, \
                        typename Traits::RangeType& _VALUE_) const

#define END_SIMPLE_BOUNDARYTYPE_FUNCTION \
  inline const GV& getGridView () \
  { \
    return gv; \
  } \
};

#define INSTATIONARY_ANALYTIC_FUNCTION(_NAME_,_POINT_,_VALUE_)  \
  template<typename GV, typename RF, typename TReal>                           \
class _NAME_ \
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>, \
                                                  _NAME_<GV,RF,TReal> >      \
{ \
  TReal _time; \
\
public: \
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits; \
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,_NAME_<GV,RF,TReal> > BaseT; \
\
  void setTime(TReal time) \
  { \
    _time = time; \
  } \
\
  TReal time() const \
  { \
    return _time; \
  } \
\
  _NAME_ (const GV& gv) : BaseT(gv) {} \
  inline void evaluateGlobal (const typename Traits::DomainType& _POINT_, \
                              typename Traits::RangeType& _VALUE_) const

#define END_INSTATIONARY_ANALYTIC_FUNCTION };

#endif // DUNE_MULTIDOMAIN_FUNCTIONMACROS_HH
