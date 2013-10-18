#ifndef DUNE_MULTIDOMAIN_INTERPOLATE_HH
#define DUNE_MULTIDOMAIN_INTERPOLATE_HH

#include <dune/common/tupleutility.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/multidomain/entitywrappers.hh>
#include <dune/pdelab/multidomain/utility.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {

// end of recursion
template<typename IB, typename XView, typename Element>
void interpolate_subproblems(const IB& ib, XView& x_view, const Element& element)
{}

template<typename IB, typename XView, typename Element, typename Function, typename SubProblemLFS, typename... SubProblems>
void interpolate_subproblems(const IB& ib,
                             XView& x_view,
                             const Element& element,
                             const Function& f,
                             SubProblemLFS& subProblemLFS,
                             SubProblems&&... subProblems)
{
  if (subProblemLFS.appliesTo(element))
    {
      subProblemLFS.bind();
      TypeTree::applyToTreePair(f,subProblemLFS,InterpolateVisitor<IB,typename Element::Entity,XView>(ib,element.entity(),x_view));
    }
  interpolate_subproblems(ib,x_view,element,subProblems...);
}

// interpolation from a given grid function
template<typename GFS, typename LFS, typename XG, typename... SubProblems>
void do_interpolate(const GFS& gfs, LFS& lfs, XG& xg, SubProblems&&... subProblems)
{

  // get some types
  typedef typename GFS::Traits::GridType::LeafGridView GV;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;

  // get index set
  const typename GV::IndexSet& is=gfs.gridView().indexSet();

  // caching
  typedef LFSIndexCache<LFS> LFSCache;
  LFSCache lfs_cache(lfs);
  typedef typename XG::template LocalView<LFSCache> XView;

  XView x_view(xg);

  // loop once over the grid
  for (ElementIterator it = gfs.gridView().template begin<0>(),
         endit = gfs.gridView().template end<0>();
       it != endit;
       ++it)
    {
      // bind local function space to element
      lfs.bind(*it);
      lfs_cache.update();
      x_view.bind(lfs_cache);

      // call interpolate
      interpolate_subproblems(InterpolateBackendStandard(),x_view,
                              ElementWrapper<GV>(*it,is.subDomains(*it)),
                              subProblems...);

      x_view.unbind();
    }

  x_view.detach();
}


template<typename GFS, typename XG, typename... SubProblems>
void interpolateOnTrialSpace(const GFS& gfs, XG& xg, SubProblems&&... subProblems)
{
  LocalFunctionSpace<GFS> lfs(gfs);
  do_interpolate(gfs,lfs,xg,extract_trial_lfs(lfs,subProblems)...);
}

template<typename GFS, typename XG, typename... SubProblems>
void interpolateOnTestSpace(const GFS& gfs, XG& xg, SubProblems&&... subProblems)
{
  LocalFunctionSpace<GFS> lfs(gfs);
  do_interpolate(gfs,lfs,xg,extract_test_lfs(lfs,subProblems)...);
}


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
    _setTime(time,TypeTree::index_range<sizeof...(T)>());
  }

  const BaseT& base() const
  {
    return static_cast<const BaseT&>(*this);
  }

private:

  template<typename R, std::size_t... i>
  void _setTime(const R& time, TypeTree::index_pack<i...>)
  {
    TypeTree::discard((std::get<i>(*this).setTime(time),0)...);
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
