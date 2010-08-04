#ifndef DUNE_MULTIDOMAIN_INTERPOLATE_HH
#define DUNE_MULTIDOMAIN_INTERPOLATE_HH

#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/multidomain/utility.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {

template<typename... SubProblems>
struct interpolation_pairs;

template<typename Function, typename SubProblemLFS, typename... SubProblems>
struct interpolation_pairs<Function,SubProblemLFS,SubProblems...>
{
  typedef interpolation_pairs<SubProblems...> next_pair;

  dune_static_assert(((is_subproblem<typename SubProblemLFS::Traits::SubProblem>::value == true)), "subproblem local function space parameter invalid");

  template<typename LFS>
  static void setupLFS(const LFS& lfs,
                       const Function& f,
                       const SubProblemLFS& subProblemLFS,
                       const SubProblems&... subProblems)
  {
    subProblemLFS.setup(lfs);
    next_pair::setupLFS(lfs,subProblems...);
  }

  template<typename LFS, typename IB, typename XG, typename Element, typename SubDomainSet>
  static void interpolate(const LFS& lfs,
                          const IB& ib,
                          XG& xg,
                          const Element& element,
                          const SubDomainSet& subDomainSet,
                          const Function& f,
                          const SubProblemLFS& subProblemLFS,
                          const SubProblems&... subProblems)
  {
    if (subProblemLFS.appliesTo(subDomainSet))
    {
      subProblemLFS.bind();
      InterpolateVisitNodeMetaProgram<IB,Function,Function::isLeaf,SubProblemLFS,SubProblemLFS::isLeaf>
        ::interpolate(ib,f,subProblemLFS,xg,element);
    }
    next_pair::interpolate(lfs,ib,xg,element,subDomainSet,subProblems...);
  }

};

template<typename UnpairedParameter>
struct interpolation_pairs<UnpairedParameter>
{
  dune_static_assert(AlwaysFalse<UnpairedParameter>::value, "incomplete Function/SubProblemLFS pair in interpolate()");
};

template<>
struct interpolation_pairs<>
{

  template<typename LFS>
  static void setupLFS(const LFS& lfs)
  {
  }

  template<typename LFS, typename IB, typename XG, typename Element, typename SubDomainSet>
  static void interpolate(const LFS& lfs,
                          const IB& ib,
                          XG& xg,
                          const Element& element,
                          const SubDomainSet& subDomainSet)
  {
  }

};


// interpolation from a given grid function
template<typename GFS, typename XG, typename... SubProblems>
void interpolate(const GFS& gfs, XG& xg, const SubProblems&... subProblems)
{
  typedef interpolation_pairs<SubProblems...> _SubProblems;

  // this is the leaf version now

  // get some types
  typedef typename GFS::Traits::GridType::LeafGridView GV;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

  // make local function space
  typedef typename GFS::LocalFunctionSpace LFS;
  LFS lfs(gfs);

  // get index set
  const typename GV::IndexSet& is=gfs.gridview().indexSet();

  _SubProblems::setupLFS(lfs,subProblems...);

  // loop once over the grid
  for (ElementIterator it = gfs.gridview().template begin<0>();
       it!=gfs.gridview().template end<0>(); ++it)
    {
      // bind local function space to element
      lfs.bind(*it);
      // call interpolate
      _SubProblems::interpolate(lfs,InterpolateBackendStandard(),xg,*it,is.subDomains(*it),subProblems...);
    }
}


template<typename GFS, typename XG, typename... SubProblems>
void interpolateOnTrialSpace(const GFS& gfs, XG& xg, const SubProblems&... subProblems)
{
  interpolate(gfs,xg,extract_trial_lfs(subProblems)...);
}

template<typename GFS, typename XG, typename... SubProblems>
void interpolateOnTestSpace(const GFS& gfs, XG& xg, const SubProblems&... subProblems)
{
  interpolate(gfs,xg,extract_test_lfs(subProblems)...);
}


// interpolation from a given grid function, using the global interface of the local finite element
template<typename GFS, typename XG, typename... SubProblems>
void interpolateGlobal(const GFS& gfs, XG& xg, const SubProblems&... subProblems)
{
  typedef interpolation_pairs<SubProblems...> _SubProblems;

  // this is the leaf version now

  // get some types
  typedef typename GFS::Traits::GridType::LeafGridView GV;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef InterpolateBackendGlobal<typename GV::Traits::template Codim<0>::Geometry> IB;

  // make local function space
  typedef typename GFS::LocalFunctionSpace LFS;
  LFS lfs(gfs);

  // get index set
  const typename GV::IndexSet& is=gfs.gridview().indexSet();

  _SubProblems::setupLFS(lfs,subProblems...);

  // loop once over the grid
  for (ElementIterator it = gfs.gridview().template begin<0>();
       it!=gfs.gridview().template end<0>(); ++it)
    {
      // bind local function space to element
      lfs.bind(*it);
      // call interpolate
      _SubProblems::interpolate(lfs,IB(it->geometry()),xg,*it,is.subDomains(*it),subProblems...);
    }
}


template<typename GFS, typename XG, typename... SubProblems>
void interpolateGlobalOnTrialSpace(const GFS& gfs, XG& xg, const SubProblems&... subProblems)
{
  interpolateGlobal(gfs,xg,extract_trial_lfs(subProblems)...);
}

template<typename GFS, typename XG, typename... SubProblems>
void interpolateGlobalOnTestSpace(const GFS& gfs, XG& xg, const SubProblems&... subProblems)
{
  interpolateGlobal(gfs,xg,extract_test_lfs(subProblems)...);
}

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_INTERPOLATE_HH
