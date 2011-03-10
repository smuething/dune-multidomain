#ifndef DUNE_PDELAB_MULTIDOMAIN_CONSTRAINTS_HH
#define DUNE_PDELAB_MULTIDOMAIN_CONSTRAINTS_HH

#include <dune/pdelab/constraints/constraints.hh>
#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/multidomain/operatorflagtests.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {



template<typename CG, typename... ConstraintsSpecifications>
class ConstraintsAssemblerEngine
  : public Dune::PDELab::TypeTree::VariadicCompositeNode<ConstraintsSpecifications...>
  , public Dune::PDELab::LocalAssemblerEngineBase
{

  typedef Dune::PDELab::TypeTree::VariadicCompositeNode<ConstraintsSpecifications...> NodeT;

public:

  bool requireIntersections() const
  {
    return requireVSkeleton() || requireVBoundary();
  }

  bool requireIntersectionsTwoSided() const
  {
    return requireVBoundary();
  }

  bool requireVVolume() const
  {
    return any_child<ConstraintsAssemblerEngine,do_constraints_volume>::value;
  }

  bool requireVSkeleton() const
  {
    return any_child<ConstraintsAssemblerEngine,do_constraints_skeleton>::value ||
      any_child<ConstraintsAssemblerEngine,do_constraints_boundary>::value ||
      any_child<ConstraintsAssemblerEngine,do_constraints_processor>::value;
  }

  bool requireVBoundary() const
  {
    return any_child<ConstraintsAssemblerEngine,do_constraints_boundary>::value;
  }

  template<typename EG, typename LFSV>
  void onBindLFSV(const EG& eg, const LFSV& lfsv)
  {
    TypeTree::applyToTree(subProblemConstraints,BindSubProblemLFS_S());
  }

  template<typename IG, typename LFSV_N>
  void onBindLFSVOutside(const IG& ig, const LFSV_N& lfsv_n)
  {
    TypeTree::applyToTree(subProblemConstraints,BindSubProblemLFS_N());
  }

  template<typename EG, typename LFSV>
  void assembleVVolume(const EG& eg, const LFSV& lfsv)
  {
    typedef visitor<functors::volume_constraints> Visitor;
    applyConstraints(Visitor::add_data(wrap_eg(eg),wrap_lfsv(lfsv),wrap_cg(cg)));
  }

  template<typename IG, typename LFSV_S, typename LFSV_N>
  void assembleVSkeleton(const IG& ig,
                         const LFSV_S& lfsv_s,
                         const LFSV_N& lfsv_n)
  {
    typedef visitor<functors::skeleton_or_processor_or_boundary_constraints> Visitor;
    applyConstraints(Visitor::add_data(wrap_ig(ig),wrap_lfsv_s(lfsv_s),wrap_lfsv_n(lfsv_n),wrap_cg(cg)));
  }

  template<typename IG, typename LFSV>
  void assembleVBoundary(const IG& ig, const LFSV& lfsv)
  {
    typedef visitor<functors::boundary_constraints> Visitor;
    applyConstraints(Visitor::add_data(wrap_eg(ig),wrap_lfsv(lfsv),wrap_cg(cg)));
  }

  ConstraintsAssemblerEngine(CG& cg_, ConstraintsSpecifications&... constraintsSpecifications)
    : NodeT(constraintsSpecifications...)
    , cg(cg_)
  {}

private:

  CG& cg;

};




//! construct constraints from given boundary condition function
/**
 * \tparam F   Type implementing a boundary condition function
 * \tparam GFS Type implementing the model GridFunctionSpace
 * \tparam CG  Type implementing the model
 *             GridFunctionSpace::ConstraintsContainer::Type
 *
 * \param f       The boundary condition function
 * \param gfs     The gridfunctionspace
 * \param cg      The constraints container
 * \param verbose Print information about the constaints at the end
 */
template<typename F, typename GFS, typename CG, typename... SubProblemBoundaries>
void constraints(const F& f, const GFS& gfs, CG& cg, const SubProblemBoundaries&... subProblemBoundaries)
{
  typedef constraints_pairs<SubProblemBoundaries...> SubProblemConstraints;

  // clear global constraints
  cg.clear();

  // get some types
  typedef typename GFS::Traits::GridType::LeafGridView GV;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;

  // make local function space
  typedef typename GFS::LocalFunctionSpace LFS;
  LFS lfs_e(gfs);
  LFS lfs_f(gfs);
  SubProblemConstraints::setupLFS(lfs_e,subProblemBoundaries...);

  // get index set
  const typename GV::IndexSet& is=gfs.gridview().indexSet();

  // helper to compute offset dependent on geometry type
  const int chunk=1<<28;
  int offset = 0;
  std::map<Dune::GeometryType,int> gtoffset;

  // loop once over the grid
  for (ElementIterator it = gfs.gridview().template begin<0>();
       it!=gfs.gridview().template end<0>(); ++it)
    {
      // assign offset for geometry type;
      if (gtoffset.find(it->type())==gtoffset.end())
        {
          gtoffset[it->type()] = offset;
          offset += chunk;
        }

      const typename GV::IndexSet::IndexType id = is.index(*it)+gtoffset[it->type()];

      // bind local function space to element
      lfs_e.bind(*it);

      ConstraintsVisitNodeMetaProgram2<LFS,LFS::isLeaf>
        ::volume(lfs_e,cg,ElementGeometry<Element>(*it));

      SubProblemConstraints::volume(lfs_e,cg,ElementGeometry<Element>(*it),is.subDomains(*it),subProblemBoundaries...);

      // iterate over intersections and call metaprogram
      unsigned int intersection_index = 0;
      IntersectionIterator endit = gfs.gridview().iend(*it);
      for (IntersectionIterator iit = gfs.gridview().ibegin(*it); iit!=endit; ++iit, ++intersection_index)
        {
          if (iit->boundary())
            {
              //ConstraintsVisitNodeMetaProgram<F,F::isLeaf,LFS,LFS::isLeaf>
              //  ::boundary(f,lfs_e,cg,IntersectionGeometry<Intersection>(*iit,intersection_index));
              SubProblemConstraints::boundary(lfs_e,cg,IntersectionGeometry<Intersection>(*iit,intersection_index),is.subDomains(*(iit->inside())),subProblemBoundaries...);
            }

          /* TODO: parallel
          // ParallelStuff: BEGIN support for processor boundaries.
          if ((!iit->boundary()) && (!iit->neighbor()))
            ConstraintsVisitNodeMetaProgram2<LFS,LFS::isLeaf>
              ::processor(lfs_e,cg,IntersectionGeometry<Intersection>(*iit,intersection_index));
          // END support for processor boundaries.
          */

          if (iit->neighbor()){

            Dune::GeometryType gtn = iit->outside()->type();
            const typename GV::IndexSet::IndexType idn = is.index(*(iit->outside()))+gtoffset[gtn];

            if(id>idn){
              // bind local function space to element in neighbor
              lfs_f.bind( *(iit->outside()) );

              //ConstraintsVisitNodeMetaProgram2<LFS,LFS::isLeaf>
              //  ::skeleton(lfs_e,lfs_f,cg,IntersectionGeometry<Intersection>(*iit,intersection_index));

              SubProblemConstraints::skeletonOrBoundary(lfs_e,lfs_f,cg,IntersectionGeometry<Intersection>(*iit,intersection_index),is.subDomains(*(iit->inside())),is.subDomains(*(iit->outside())),subProblemBoundaries...);
            }
          }
        }
    }

  /*
  // print result
  if(verbose){
    std::cout << "constraints:" << std::endl;
    typedef typename CG::iterator global_col_iterator;
    typedef typename CG::value_type::second_type global_row_type;
    typedef typename global_row_type::iterator global_row_iterator;

    std::cout << cg.size() << " constrained degrees of freedom" << std::endl;

    for (global_col_iterator cit=cg.begin(); cit!=cg.end(); ++cit)
      {
        std::cout << cit->first << ": ";
        for (global_row_iterator rit=(cit->second).begin(); rit!=(cit->second).end(); ++rit)
          std::cout << "(" << rit->first << "," << rit->second << ") ";
        std::cout << std::endl;
      }
      } */

} // constraints


/**
 * Construct constraints on the trial grid function space.
 */
template<typename F, typename GFS, typename CG, typename... SubProblemBoundaries>
void trialSpaceConstraints(const F& f, const GFS& gfs, CG& cg, const SubProblemBoundaries&... subProblemBoundaries)
{
  constraints(f,gfs,cg,extract_trial_lfs(subProblemBoundaries)...);
}


/**
 * Construct constraints on the test grid function space.
 */
template<typename F, typename GFS, typename CG, typename... SubProblemBoundaries>
void testSpaceConstraints(const F& f, const GFS& gfs, CG& cg, const SubProblemBoundaries&... subProblemBoundaries)
{
  constraints(f,gfs,cg,extract_test_lfs(subProblemBoundaries)...);
}


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_CONSTRAINTS_HH
