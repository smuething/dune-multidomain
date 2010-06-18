#ifndef DUNE_MULTIDOMAIN_CONSTRAINTS_HH
#define DUNE_MULTIDOMAIN_CONSTRAINTS_HH

#include <dune/pdelab/gridfunctionspace/constraints.hh>
#include <type_traits>

namespace Dune {
namespace PDELab {
namespace MultiDomain {

template<typename... SubProblemBoundaries>
struct constraints_pairs;

template<typename BoundaryConditionTypeFunction, typename SubProblemLFS, typename... SubProblemBoundaries>
struct constraints_pairs<BoundaryConditionTypeFunction,SubProblemLFS,SubProblemBoundaries...>
{
  typedef constraints_pairs<SubProblemBoundaries...> next_pair;

  dune_static_assert(((tag<typename SubProblemLFS::Traits::SubProblem>::value == true)), "subproblem local function space parameter invalid");

  template<typename LFS, typename CG, typename Geometry, typename SubDomainSet>
  static void volume(const LFS& lfs,
                     CG& cg,
                     const Geometry& geometry,
                     const SubDomainSet& subDomainSet,
                     const BoundaryConditionTypeFunction& boundaryType,
                     const SubProblemLFS& subProblemLFS,
                     const SubProblemBoundaries&... subProblemBoundaries)
  {
    if (subProblemLFS.appliesTo(subDomainSet))
    {
      ConstraintsVisitNodeMetaProgram2<SubProblemLFS,SubProblemLFS::isLeaf>
        ::volume(subProblemLFS,cg,geometry);
    }
    next_pair::volume(lfs,cg,geometry,subDomainSet,subProblemBoundaries...);
  }

  template<typename LFS, typename CG, typename IntersectionGeometry, typename SubDomainSet>
  static void do_boundary(const LFS& lfs,
                          CG& cg,
                          const IntersectionGeometry& intersectionGeometry,
                          const SubDomainSet& subDomainSet,
                          const BoundaryConditionTypeFunction& boundaryType,
                          const SubProblemLFS& subProblemLFS)
  {
    ConstraintsVisitNodeMetaProgram<BoundaryConditionTypeFunction,
                                    BoundaryConditionTypeFunction::isLeaf,
                                    SubProblemLFS,
                                    SubProblemLFS::isLeaf>
      ::boundary(boundaryType,subProblemLFS,cg,intersectionGeometry);
  }

  template<typename LFS, typename CG, typename IntersectionGeometry, typename SubDomainSet>
  static void boundary(const LFS& lfs,
                       CG& cg,
                       const IntersectionGeometry& intersectionGeometry,
                       const SubDomainSet& subDomainSet,
                       const BoundaryConditionTypeFunction& boundaryType,
                       const SubProblemLFS& subProblemLFS,
                       const SubProblemBoundaries&... subProblemBoundaries)
  {
    if (subProblemLFS.appliesTo(subDomainSet))
    {
      do_boundary(lfs,cg,intersectionGeometry,subDomainSet,boundaryType,subProblemLFS);
    }
    next_pair::boundary(lfs,cg,intersectionGeometry,subDomainSet,subProblemBoundaries...);
  }

  template<typename LFS, typename CG, typename IntersectionGeometry, typename SubDomainSet>
  static void skeletonOrBoundary(const LFS& lfs,
                                 const LFS& lfs_n,
                                 CG& cg,
                                 const IntersectionGeometry& intersectionGeometry,
                                 const SubDomainSet& subDomainSet,
                                 const SubDomainSet& neighborSubDomainSet,
                                 const BoundaryConditionTypeFunction& boundaryType,
                                 const SubProblemLFS& subProblemLFS,
                                 const SubProblemBoundaries&... subProblemBoundaries)
  {
    if (subProblemLFS.appliesTo(subDomainSet))
    {
      if (subProblemLFS.appliesTo(neighborSubDomainSet))
      {
        SubProblemLFS subProblemLFS_neighbor(lfs_n,subProblemLFS.subProblem(),subProblemLFS.constraints());
        ConstraintsVisitNodeMetaProgram2<SubProblemLFS,SubProblemLFS::isLeaf>
          ::skeleton(subProblemLFS,subProblemLFS_neighbor,cg,intersectionGeometry);
      } else {
        do_boundary(lfs,cg,intersectionGeometry,subDomainSet,boundaryType,subProblemLFS);
      }
    }
    next_pair::skeletonOrBoundary(lfs,lfs_n,cg,intersectionGeometry,subDomainSet,neighborSubDomainSet,subProblemBoundaries...);
  }

};

template<typename UnpairedParameter>
struct constraints_pairs<UnpairedParameter>
{
  dune_static_assert(AlwaysFalse<UnpairedParameter>::value, "incomplete BoundaryTypeFunction/SubProblemLFS pair in constraints()");
};

template<>
struct constraints_pairs<>
{

  template<typename LFS, typename CG, typename Geometry, typename SubDomainSet>
  static void volume(const LFS& lfs,
                     CG& cg,
                     const Geometry& geometry,
                     const SubDomainSet& subDomainSet)
  {
  }

  template<typename LFS, typename CG, typename IntersectionGeometry, typename SubDomainSet>
  static void boundary(const LFS& lfs,
                       CG& cg,
                       const IntersectionGeometry& intersectionGeometry,
                       const SubDomainSet& subDomainSet)
  {
  }

  template<typename LFS, typename CG, typename IntersectionGeometry, typename SubDomainSet>
  static void skeletonOrBoundary(const LFS& lfs,
                                 const LFS& lfs_n,
                                 CG& cg,
                                 const IntersectionGeometry& intersectionGeometry,
                                 const SubDomainSet& subDomainSet,
                                 const SubDomainSet& neighborsubDomainSet)
  {
  }

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
  typedef typename GFS::Traits::GridViewType GV;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;

  // make local function space
  typedef typename GFS::LocalFunctionSpace LFS;
  LFS lfs_e(gfs);
  LFS lfs_f(gfs);

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

      SubProblemConstraints::volume(lfs_e,cg,*it,subProblemBoundaries...);

      // iterate over intersections and call metaprogram
      unsigned int intersection_index = 0;
      IntersectionIterator endit = gfs.gridview().iend(*it);
      for (IntersectionIterator iit = gfs.gridview().ibegin(*it); iit!=endit; ++iit, ++intersection_index)
        {
          if (iit->boundary())
            {
              ConstraintsVisitNodeMetaProgram<F,F::isLeaf,LFS,LFS::isLeaf>
                ::boundary(f,lfs_e,cg,IntersectionGeometry<Intersection>(*iit,intersection_index));
              SubProblemConstraints::boundary(lfs_e,cg,*iit,intersection_index,subProblemBoundaries...);
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

              ConstraintsVisitNodeMetaProgram2<LFS,LFS::isLeaf>
                ::skeleton(lfs_e,lfs_f,cg,IntersectionGeometry<Intersection>(*iit,intersection_index));

              SubProblemConstraints::skeletonOrBoundary(lfs_e,lfs_f,cg,*iit,intersection_index,subProblemBoundaries...);
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

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_CONSTRAINTS_HH
