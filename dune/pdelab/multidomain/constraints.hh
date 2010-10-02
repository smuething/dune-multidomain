#ifndef DUNE_MULTIDOMAIN_CONSTRAINTS_HH
#define DUNE_MULTIDOMAIN_CONSTRAINTS_HH

#include <map>
#include <dune/pdelab/gridfunctionspace/constraints.hh>
#include <dune/pdelab/multidomain/utility.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {


//FIXME: This should be handled better....

template<typename C, typename F, bool FisLeaf, typename LFS, bool LFSisLeaf>
struct SubProblemConstraintsVisitNodeMetaProgram;

template<typename C, typename F, typename LFS, int n, int i>
struct SubProblemConstraintsVisitChildMetaProgram // visit i'th child of inner node
{
  template<typename CG, typename I>
  static void boundary (const C& c, const F& f, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
  {
    // vist children of both nodes in pairs
    typedef typename F::template Child<i>::Type FC;
    typedef typename LFS::template Child<i>::Type LFSC;

    const FC& fc=f.template getChild<i>();
    const LFSC& lfsc=lfs.template getChild<i>();

    SubProblemConstraintsVisitNodeMetaProgram<C,FC,FC::isLeaf,LFSC,LFSC::isLeaf>::boundary(c,fc,lfsc,cg,ig);
    SubProblemConstraintsVisitChildMetaProgram<C,F,LFS,n,i+1>::boundary(c,f,lfs,cg,ig);
  }
};

template<typename C, typename F, typename LFS, int n>
struct SubProblemConstraintsVisitChildMetaProgram<C,F,LFS,n,n> // end of child recursion
{
  // end of child recursion
  template<typename CG, typename I>
  static void boundary (const C& c, const F& f, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
  {
    return;
  }
};

template<typename C, typename F, bool FisLeaf, typename LFS, bool LFSisLeaf>
struct SubProblemConstraintsVisitNodeMetaProgram // visit inner node
{
  template<typename CG, typename I>
  static void boundary (const C& c, const F& f, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
  {
    // both are inner nodes, visit all children
    // check that both have same number of children
    dune_static_assert((static_cast<int>(F::CHILDREN)==static_cast<int>(LFS::CHILDREN)),
                       "both nodes must have same number of children");

    // start child recursion
    SubProblemConstraintsVisitChildMetaProgram<C,F,LFS,F::CHILDREN,0>::boundary(c,f,lfs,cg,ig);
  }
};

template<typename C, typename F, typename LFS>
struct SubProblemConstraintsVisitNodeMetaProgram<C,F,true,LFS,false> // try to interpolate components from vector valued function
{
  template<typename CG, typename I>
  static void boundary (const C& c, const F& f, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
  {
    dune_static_assert((static_cast<int>(LFS::isPower)==1),
                       "specialization only for power");
    dune_static_assert((static_cast<int>(LFS::template Child<0>::Type::isLeaf)==1),
                       "children must be leaves");
    dune_static_assert((static_cast<int>(F::Traits::dimRange)==static_cast<int>(LFS::CHILDREN)),
                       "number of components must coincide with number of children");

    // extract constraints type
    //typedef typename LFS::template Child<0>::Type::Traits::ConstraintsType C;

    for (int k=0; k<LFS::CHILDREN; k++)
      {
        // allocate empty local constraints map
        CG cl;

        // call boundary condition evaluation of child k with component k
        typedef BoundaryGridFunctionSelectComponentAdapter<F> FCOMP;
        FCOMP fcomp(f,k);

        ConstraintsCallBoundary<C,C::doBoundary>::boundary(c, // lfs.getChild(k).constraints(),
                                                           fcomp,ig,lfs.getChild(k),cl);

        // write coefficients into local vector
        lfs.getChild(k).mwrite(cl,cg);
      }
  }
};

template<typename C, typename F, typename LFS>
struct SubProblemConstraintsVisitNodeMetaProgram<C,F,true,LFS,true> // leaf node in both trees
{
  template<typename CG, typename I>
  static void boundary (const C& c, const F& f, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
  {
    // now we are at a single component local function space
    // which is part of a multi component local function space

    // allocate local constraints map
    CG cl;

    // extract constraints type
    // typedef typename LFS::Traits::ConstraintsType C;

    // iterate over boundary, need intersection iterator
    ConstraintsCallBoundary<C,C::doBoundary>::boundary(c /*lfs.constraints()*/,f,ig,lfs,cl);

    // write coefficients into local vector
    lfs.mwrite(cl,cg);
  }
};


// second metaprogram that iterates over local function space only

template<typename C, typename LFS, bool LFSisLeaf>
struct SubProblemConstraintsVisitNodeMetaProgram2;

template<typename C, typename LFS, int n, int i>
struct SubProblemConstraintsVisitChildMetaProgram2 // visit i'th child of inner node
{
  template<typename CG, typename I>
  static void processor (const C&c, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
  {
    typedef typename LFS::template Child<i>::Type LFSC;
    const LFSC& lfsc=lfs.template getChild<i>();

    SubProblemConstraintsVisitNodeMetaProgram2<C,LFSC,LFSC::isLeaf>::processor(c,lfsc,cg,ig);
    SubProblemConstraintsVisitChildMetaProgram2<C,LFS,n,i+1>::processor(c,lfs,cg,ig);
  }
  template<typename CG, typename I>
  static void skeleton (const C& c, const LFS& lfs_e, const LFS& lfs_f, CG& cg,
                        const IntersectionGeometry<I>& ig)
  {
    typedef typename LFS::template Child<i>::Type LFSC;
    const LFSC& lfsc_e=lfs_e.template getChild<i>();
    const LFSC& lfsc_f=lfs_f.template getChild<i>();

    SubProblemConstraintsVisitNodeMetaProgram2<C,LFSC,LFSC::isLeaf>::skeleton(c,lfsc_e,lfsc_f,cg,ig);
    SubProblemConstraintsVisitChildMetaProgram2<C,LFS,n,i+1>::skeleton(c,lfs_e,lfs_f,cg,ig);
  }
  template<typename CG, typename E>
  static void volume (const C& c, const LFS& lfs, CG& cg, const ElementGeometry<E>& eg)
  {
    typedef typename LFS::template Child<i>::Type LFSC;
    const LFSC& lfsc=lfs.template getChild<i>();

    SubProblemConstraintsVisitNodeMetaProgram2<C,LFSC,LFSC::isLeaf>::volume(c,lfsc,cg,eg);
    SubProblemConstraintsVisitChildMetaProgram2<C,LFS,n,i+1>::volume(c,lfs,cg,eg);
  }
};

template<typename C, typename LFS, int n>
struct SubProblemConstraintsVisitChildMetaProgram2<C,LFS,n,n> // end of child recursion
{
  template<typename CG, typename I>
  static void processor (const C& c, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
  {
    return;
  }
  template<typename CG, typename I>
  static void skeleton (const C& c, const LFS& lfs_e, const LFS& lfs_f, CG& cg,
                        const IntersectionGeometry<I>& ig)
  {
    return;
  }
  template<typename CG, typename E>
  static void volume (const C& c, const LFS& lfs, CG& cg, const ElementGeometry<E>& eg)
  {
    return;
  }
};


template<typename C, typename LFS, bool LFSisLeaf>
struct SubProblemConstraintsVisitNodeMetaProgram2 // visit inner node
{
  template<typename CG, typename I>
  static void processor (const C& c, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
  {
    // start child recursion
    SubProblemConstraintsVisitChildMetaProgram2<C,LFS,LFS::CHILDREN,0>::processor(c,lfs,cg,ig);
  }
  template<typename CG, typename I>
  static void skeleton (const C& c, const LFS& lfs_e, const LFS& lfs_f, CG& cg,
                        const IntersectionGeometry<I>& ig)
  {
    // start child recursion
    SubProblemConstraintsVisitChildMetaProgram2<C,LFS,LFS::CHILDREN,0>::skeleton(c,lfs_e,lfs_f,cg,ig);
  }
  template<typename CG, typename E>
  static void volume (const C& c, const LFS& lfs, CG& cg, const ElementGeometry<E>& eg)
  {
    // start child recursion
    SubProblemConstraintsVisitChildMetaProgram2<C,LFS,LFS::CHILDREN,0>::volume(c,lfs,cg,eg);
  }
};


template<typename C, typename LFS>
struct SubProblemConstraintsVisitNodeMetaProgram2<C,LFS,true> // leaf node
{
  template<typename CG, typename I>
  static void processor (const C& c, const LFS& lfs, CG& cg, const IntersectionGeometry<I>& ig)
  {
    // now we are at a single component local function space
    // which is part of a multi component local function space

    // allocate local constraints map
    CG cl;

    // extract constraints type
    // typedef typename LFS::Traits::ConstraintsType C;

    // iterate over boundary, need intersection iterator
    ConstraintsCallProcessor<C,C::doProcessor>::processor(c /*lfs.constraints()*/,ig,lfs,cl);

    // write coefficients into local vector
    lfs.mwrite(cl,cg);
  }
  template<typename CG, typename I>
  static void skeleton (const C& c, const LFS& lfs_e, const LFS& lfs_f, CG& cg,
                        const IntersectionGeometry<I>& ig)
  {
    // now we are at a single component local function space
    // which is part of a multi component local function space

    // allocate local constraints map for both elements adjacent
    // to this intersection
    CG cl_e;
    CG cl_f;

    // extract constraints type
    // typedef typename LFS::Traits::ConstraintsType C;

    // as LFS::constraints() just returns the constraints of the
    // GridFunctionSpace, lfs_e.constraints() is equivalent to
    // lfs_f.constraints()
    //const C & c = lfs_e.constraints();

    // iterate over boundary, need intersection iterator
    ConstraintsCallSkeleton<C,C::doSkeleton>::skeleton(c,ig,lfs_e,lfs_f,cl_e,cl_f);

    // write coefficients into local vector
    lfs_e.mwrite(cl_e,cg);
    lfs_f.mwrite(cl_f,cg);
  }
  template<typename CG, typename E>
  static void volume (const C& c, const LFS& lfs, CG& cg, const ElementGeometry<E>& eg)
  {
    // now we are at a single component local function space
    // which is part of a multi component local function space

    // allocate local constraints map
    CG cl;

    // extract constraints type
    //typedef typename LFS::Traits::ConstraintsType C;
    //const C & c = lfs.constraints();

    // iterate over boundary, need intersection iterator
    ConstraintsCallVolume<C,C::doVolume>::volume(c,eg,lfs,cl);

    // write coefficients into local vector
    lfs.mwrite(cl,cg);
  }
};


template<typename... SubProblemBoundaries>
struct constraints_pairs;

template<typename BoundaryConditionTypeFunction, typename SubProblemLFS, typename... SubProblemBoundaries>
struct constraints_pairs<BoundaryConditionTypeFunction,SubProblemLFS,SubProblemBoundaries...>
{
  typedef constraints_pairs<SubProblemBoundaries...> next_pair;

  dune_static_assert(((is_subproblem<typename SubProblemLFS::Traits::SubProblem>::value == true)), "subproblem local function space parameter invalid");

  template<typename LFS>
  static void setupLFS(const LFS& lfs,
                       const BoundaryConditionTypeFunction& boundaryType,
                       const SubProblemLFS& subProblemLFS,
                       const SubProblemBoundaries&... subProblemBoundaries)
  {
    subProblemLFS.setup(lfs);
    next_pair::setupLFS(lfs,subProblemBoundaries...);
  }

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
      subProblemLFS.bind();
      typedef typename SubProblemLFS::Traits::Constraints Constraints;
      const Constraints& constraints = subProblemLFS.constraints();
      SubProblemConstraintsVisitNodeMetaProgram2<Constraints,SubProblemLFS,SubProblemLFS::isLeaf>
        ::volume(constraints,subProblemLFS,cg,geometry);
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
    typedef typename SubProblemLFS::Traits::Constraints Constraints;
    const Constraints& constraints = subProblemLFS.constraints();
    SubProblemConstraintsVisitNodeMetaProgram<Constraints,
                                              BoundaryConditionTypeFunction,
                                              BoundaryConditionTypeFunction::isLeaf,
                                              SubProblemLFS,
                                              SubProblemLFS::isLeaf>
      ::boundary(constraints,boundaryType,subProblemLFS,cg,intersectionGeometry);
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
      subProblemLFS.bind();
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
        subProblemLFS.bind();
        SubProblemLFS subProblemLFS_neighbor(lfs_n,subProblemLFS.subProblem(),subProblemLFS.constraints());
        subProblemLFS_neighbor.bind();
        typedef typename SubProblemLFS::Traits::Constraints Constraints;
        const Constraints& constraints = subProblemLFS.constraints();
        SubProblemConstraintsVisitNodeMetaProgram2<Constraints,SubProblemLFS,SubProblemLFS::isLeaf>
          ::skeleton(constraints,subProblemLFS,subProblemLFS_neighbor,cg,intersectionGeometry);
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
  dune_static_assert((AlwaysFalse<UnpairedParameter>::value), "incomplete BoundaryTypeFunction/SubProblemLFS pair in constraints()");
};

template<>
struct constraints_pairs<>
{

  template<typename LFS>
  static void setupLFS(const LFS& lfs)
  {
  }

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
                                 const SubDomainSet& neighborSubDomainSet)
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

#endif // DUNE_MULTIDOMAIN_CONSTRAINTS_HH
