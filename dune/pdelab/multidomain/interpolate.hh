#ifndef DUNE_MULTIDOMAIN_INTERPOLATE_HH
#define DUNE_MULTIDOMAIN_INTERPOLATE_HH

#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/multidomain/utility.hh>
#include <dune/common/tupleutility.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {

// end of recursion
template<typename IB, typename XG, typename Element>
void interpolate_subproblems(const IB& ib, XG& xg, const Element& element)
{}

template<typename IB, typename XG, typename Element, typename Function, typename SubProblemLFS, typename... SubProblems>
void interpolate_subproblems(const IB& ib,
                             XG& xg,
                             const Element& element,
                             const Function& f,
                             SubProblemLFS& subProblemLFS,
                             SubProblems&&... subProblems)
{
  if (subProblemLFS.appliesTo(element))
    {
      subProblemLFS.bind();
      Dune::PDELab::TypeTree::applyToTreePair(f,subProblemLFS,InterpolateVisitor<IB,typename Element::Entity,XG>(ib,element.entity(),xg));
    }
  interpolate_subproblems(ib,xg,element,subProblems...);
}

// interpolation from a given grid function
template<typename GFS, typename LFS, typename XG, typename... SubProblems>
void do_interpolate(const GFS& gfs, LFS& lfs, XG& xg, SubProblems&&... subProblems)
{

  // get some types
  typedef typename GFS::Traits::GridType::LeafGridView GV;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

  // get index set
  const typename GV::IndexSet& is=gfs.gridview().indexSet();

  // loop once over the grid
  for (ElementIterator it = gfs.gridview().template begin<0>();
       it!=gfs.gridview().template end<0>(); ++it)
    {
      // bind local function space to element
      lfs.bind(*it);
      // call interpolate
      interpolate_subproblems(InterpolateBackendStandard(),xg,
                              ElementWrapper<GV>(*it,is.subDomains(*it)),
                              subProblems...);
    }
}


template<typename GFS, typename XG, typename... SubProblems>
void interpolateOnTrialSpace(const GFS& gfs, XG& xg, const SubProblems&... subProblems)
{
  LocalFunctionSpace<GFS> lfs(gfs);
  do_interpolate(gfs,lfs,xg,extract_trial_lfs(lfs,subProblems)...);
}

template<typename GFS, typename XG, typename... SubProblems>
void interpolateOnTestSpace(const GFS& gfs, XG& xg, const SubProblems&... subProblems)
{
  LocalFunctionSpace<GFS> lfs(gfs);
  do_interpolate(gfs,lfs,xg,extract_test_lfs(lfs,subProblems)...);
}



template <typename R>
struct interpolation_descriptors_set_time {

  template <class T>
  void visit(T& elem) const
  {
    elem.setTime(time);
  }

  interpolation_descriptors_set_time(const R& t)
    : time(t)
  {}

  const R time;
};


template<typename... T>
struct interpolation_descriptors
  : public tuple<T&...>
{

  typedef tuple<T&...> BaseT;

  interpolation_descriptors(T&... t)
    : BaseT(t...)
  {}

  template<typename R>
  void setTime(const R& time)
  {
    Dune::ForEachValue<BaseT> forEach(*this);
    interpolation_descriptors_set_time<R> f(time);
    forEach.apply(f);
  }

  const BaseT& base() const
  {
    return static_cast<const BaseT&>(*this);
  }

};

template<typename... T>
interpolation_descriptors<T...>
interpolateOnSubProblems(T&... t)
{
  return interpolation_descriptors<T...>(t...);
}

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_INTERPOLATE_HH
