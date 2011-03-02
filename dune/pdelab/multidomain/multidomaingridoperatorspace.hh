// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_MULTIDOMAIN_MULTIDOMAINGRIDOPERATORSPACE_HH
#define DUNE_MULTIDOMAIN_MULTIDOMAINGRIDOPERATORSPACE_HH

#include<map>
#include<tuple>

#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>

#include <dune/pdelab/common/geometrywrapper.hh>
//#include"../gridfunctionspace/gridfunctionspace.hh"
#include <dune/pdelab/gridfunctionspace/constraints.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>

#include <dune/pdelab/multidomain/multidomaingridoperatorspaceutilities.hh>
#include <dune/pdelab/multidomain/operatorapplier.hh>

#include <dune/pdelab/multidomain/multidomaingridoperatorspaceinvocationhelpers.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


  void residual (const X& x, R& r) const
  {
    // visit each face only once
    const int chunk=1<<28;
    int offset = 0;
    const typename GV::IndexSet& is=gfsu.gridview().indexSet();
    std::map<Dune::GeometryType,int> gtoffset;

    // make local function spaces
    typedef typename GFSU::LocalFunctionSpace LFSU;
    LFSU lfsu(gfsu);
    typedef typename GFSV::LocalFunctionSpace LFSV;
    LFSV lfsv(gfsv);

    operator_applier_all_data apply_operator(*this);
    apply_operator.setlfsu(lfsu);
    apply_operator.setlfsv(lfsv);

    // traverse grid view
    for (ElementIterator it = gfsu.gridview().template begin<0>();
         it!=gfsu.gridview().template end<0>(); ++it)
      {
        // assign offset for geometry type;
        if (gtoffset.find(it->type())==gtoffset.end())
          {
            gtoffset[it->type()] = offset;
            offset += chunk;
          }

        // compute unique id
        int id = is.index(*it)+gtoffset[it->type()];

        // skip ghost and overlap
        if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
          continue;

        // bind local function spaces to element
        lfsu.bind(*it);
        lfsv.bind(*it);
        apply_operator.setElement(*it);
        apply_operator.setElementSubDomains(is.subDomains(*it));

        // allocate local data container
        typedef std::vector<typename X::ElementType> XL;
        XL xl(lfsu.size());
        typedef std::vector<typename R::ElementType> RL;
        RL rl(lfsv.size(),0.0);

        // read coefficents
        lfsu.vread(x,xl);

        // volume evaluation
        apply_operator.template conditional<SubProblems,do_alpha_volume<> >(InvokeAlphaVolume<XL,RL>(xl,rl));
        apply_operator.template conditional<SubProblems,do_lambda_volume<> >(InvokeLambdaVolume<RL>(rl));

        // skip if no intersection iterator is needed
        if (any_child<MultiDomainGridOperatorSpace,SubProblems,do_alpha_skeleton<> >::value ||
            any_child<MultiDomainGridOperatorSpace,SubProblems,do_alpha_boundary<> >::value ||
            any_child<MultiDomainGridOperatorSpace,SubProblems,do_lambda_boundary<> >::value ||
            any_child<MultiDomainGridOperatorSpace,Couplings,do_alpha_coupling<> >::value ||
            any_child<MultiDomainGridOperatorSpace,Couplings,do_alpha_enriched_coupling<> >::value)
          {
            // local function spaces in neighbor
            LFSU lfsun(gfsu);
            LFSV lfsvn(gfsv);
            apply_operator.setlfsun(lfsun);
            apply_operator.setlfsvn(lfsvn);

            typedef typename GFSU::CouplingLocalFunctionSpace CouplingLFSU;
            CouplingLFSU couplinglfsu(gfsu);
            typedef typename GFSV::CouplingLocalFunctionSpace CouplingLFSV;
            CouplingLFSV couplinglfsv(gfsv);
            apply_operator.setcouplinglfsu(couplinglfsu);
            apply_operator.setcouplinglfsv(couplinglfsv);

            // traverse intersections
            unsigned int intersection_index = 0;
            IntersectionIterator endit = gfsu.gridview().iend(*it);
            for (IntersectionIterator iit = gfsu.gridview().ibegin(*it);
                 iit!=endit; ++iit, ++intersection_index)
              {
                apply_operator.setIntersection(*iit,intersection_index);
                // skeleton term
                if (iit->neighbor())
                  {
                    // assign offset for geometry type;
                    Dune::GeometryType gtn = iit->outside()->type();
                    if (gtoffset.find(gtn)==gtoffset.end())
                      {
                        gtoffset[gtn] = offset;
                        offset += chunk;
                      }

                    // compute unique id for neighbor
                    int idn = is.index(*(iit->outside()))+gtoffset[gtn];
                    lfsun.bind(*(iit->outside()));
                    lfsvn.bind(*(iit->outside()));
                    apply_operator.setNeighborSubDomains(gfsu.gridview().indexSet().subDomains(*(iit->outside())));

                    // allocate local data container
                    XL xn(lfsun.size());
                    RL rn(lfsvn.size(),0.0);

                    XL xcoupling;
                    RL rcoupling;

                    // only bind the coupling local function spaces if anyone is going to use them
                    if (any_child<MultiDomainGridOperatorSpace,Couplings,do_alpha_enriched_coupling<> >::value)
                      {
                        couplinglfsu.bind(*iit);
                        couplinglfsv.bind(*iit);
                        couplinglfsu.vread(x,xcoupling);
                        rcoupling.resize(couplinglfsv.size(),0.0);
                      }

                    // read coefficents
                    lfsun.vread(x,xn);

                    // unique vist of intersection
                    apply_operator.template conditional<SubProblems,do_alpha_skeleton_or_boundary<> >
                      (InvokeAlphaSkeletonOrBoundary<XL,RL>(xl,xn,rl,rn,
                                                            id > idn ||
                                                            (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                            )
                       );
                    apply_operator.template conditional<Couplings,do_alpha_coupling<> >
                      (InvokeAlphaCoupling<XL,RL>(xl,xn,rl,rn));
                    apply_operator.template conditional<Couplings,do_alpha_enriched_coupling<> >
                      (InvokeAlphaEnrichedCoupling<XL,RL>(xl,xn,xcoupling,rl,rn,rcoupling));
                    if (apply_operator.alphaSkeletonInvoked())
                      {
                        lfsvn.vadd(rn,r);
                        apply_operator.clearAlphaSkeletonInvoked();
                      }
                    if (apply_operator.alphaEnrichedCouplingInvoked())
                      {
                        couplinglfsv.vadd(rcoupling,r);
                        apply_operator.clearAlphaEnrichedCouplingInvoked();
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    apply_operator.template conditional<SubProblems,do_alpha_boundary<> >(InvokeAlphaBoundary<XL,RL>(xl,rl));
                    apply_operator.template conditional<SubProblems,do_lambda_boundary<> >(InvokeLambdaBoundary<RL>(rl));
                  }
              }
          }

        apply_operator.template conditional<SubProblems,do_alpha_volume_post_skeleton<> >(InvokeAlphaVolumePostSkeleton<XL,RL>(xl,rl));
        apply_operator.template conditional<SubProblems,do_lambda_volume_post_skeleton<> >(InvokeLambdaVolumePostSkeleton<RL>(rl));

        // accumulate result (note: r needs to be cleared outside)
        lfsv.vadd(rl,r);
      }

    // set residual to zero on constrained dofs
    Dune::PDELab::constrain_residual(*pconstraintsv,r);












//================================================
// The operator
//================================================

//! The generic assembler ...
/**
 * \tparam GFSU GridFunctionSpace for ansatz functions
 * \tparam GFSV GridFunctionSpace for test functions
 * \tparam LP   local pattern assembler (provided by user)
 * \tparam LA   local operator assembler (provided by user)
 */
template<typename GFSU, typename GFSV,
         typename B,
         typename... SubProblemsAndCouplings>
class MultiDomainGridOperatorSpace : public VariadicCompositeNode<CopyStoragePolicy,SubProblemsAndCouplings...>
{

  template<typename,typename>
  friend struct BuildVolumePattern;

  template<typename,typename>
  friend struct BuildSkeletonPattern;

  template<typename,typename>
  friend struct BuildCouplingPattern;

  template<typename,typename>
  friend struct BuildEnrichedCouplingPattern;

  typedef VariadicCompositeNode<CopyStoragePolicy,SubProblemsAndCouplings...> BaseT;

  typedef typename extract_problems<SubProblemsAndCouplings...>::type SubProblemList;
  typedef typename extract_problems<SubProblemsAndCouplings...>::map_type SubProblemMap;
  typedef typename extract_couplings<SubProblemsAndCouplings...>::type CouplingList;

  static const std::size_t subProblemCount = std::tuple_size<SubProblemList>::value;
  static const std::size_t couplingCount = std::tuple_size<CouplingList>::value;

  typedef typename GFSU::template ConstraintsContainer<double>::Type CU;
  typedef typename GFSV::template ConstraintsContainer<double>::Type CV;

  template<std::size_t k>
  struct SubProblem {
    typedef typename std::tuple_element<k,SubProblemList>::type::type Type;
  };

  template<std::size_t k>
  const typename SubProblem<k>::Type& subProblem() const {
    return this->template getChild<std::tuple_element<k,SubProblemList>::type::globalPos>();
  }

  template<std::size_t k>
  struct Coupling {
    typedef typename std::tuple_element<k,CouplingList>::type::type Type;
  };

  template<std::size_t k>
  const typename Coupling<k>::Type& coupling() const {
    return this->template getChild<std::tuple_element<k,CouplingList>::type::globalPos>();
  }

  struct SubProblems
  {

    static const std::size_t N = subProblemCount;

    template<std::size_t k>
    struct GetType
    {
      dune_static_assert(k < N, "subProblem index too large");
      typedef typename SubProblem<k>::Type Type;
    };

    template<std::size_t k>
    static const typename SubProblem<k>::Type& get(const MultiDomainGridOperatorSpace& mgos)
    {
      dune_static_assert(k < N, "subProblem index too large");
      return mgos.template subProblem<k>();
    }

    template<std::size_t k>
    static typename SubProblem<k>::Type& get(MultiDomainGridOperatorSpace& mgos)
    {
      dune_static_assert(k < N, "subProblem index too large");
      return mgos.template subProblem<k>();
    }

  };

  struct Couplings
  {

    static const std::size_t N = couplingCount;

    template<std::size_t k>
    struct GetType
    {
      dune_static_assert(k < N, "coupling index too large");
      typedef typename Coupling<k>::Type Type;
    };

    template<std::size_t k>
    static const typename Coupling<k>::Type& get(const MultiDomainGridOperatorSpace& mgos)
    {
      dune_static_assert(k < N, "coupling index too large");
      return mgos.template coupling<k>();
    }

    template<std::size_t k>
    static typename Coupling<k>::Type& get(MultiDomainGridOperatorSpace& mgos)
    {
      dune_static_assert(k < N, "coupling index too large");
      return mgos.template coupling<k>();
    }
  };

  struct AllChildren
  {

    static const std::size_t N = BaseT::CHILDREN;

    template<std::size_t k>
    struct GetType
    {
      dune_static_assert(k < N, "child index too large");
      typedef typename BaseT::template Child<k>::Type Type;
    };

    template<std::size_t k>
    static const typename BaseT::template Child<k>::Type& get(const MultiDomainGridOperatorSpace& mgos)
    {
      dune_static_assert(k < N, "child index too large");
      return mgos.template getChild<k>();
    }

    template<std::size_t k>
    static typename BaseT::template Child<k>::Type& get(MultiDomainGridOperatorSpace& mgos)
    {
      dune_static_assert(k < N, "child index too large");
      return mgos.template getChild<k>();
    }

  };


  // extract useful types
  typedef typename GFSU::Traits::GridType Grid;
  typedef typename Grid::LeafGridView GV;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;

public:
  typedef MultiDomainGridOperatorSpaceTraits<GFSU,GFSV,B> Traits;

  template<typename E>
  struct MatrixContainer
  {
    //! \brief define Type as the Type of a Matrix of E's
    typedef typename B::template Matrix<MultiDomainGridOperatorSpace,E> Type;
  private:
    MatrixContainer () {}
  };

  //! construct GridOperatorSpace
  MultiDomainGridOperatorSpace (const GFSU& gfsu_, const GFSV& gfsv_, const CU& cu, const CV& cv, SubProblemsAndCouplings&... subProblems_)
    :  BaseT(subProblems_...), gfsu(gfsu_), gfsv(gfsv_)
  {
    pconstraintsu = &cu;
    pconstraintsv = &cv;
  }

  //! get dimension of space u
  typename GFSU::Traits::SizeType globalSizeU () const
  {
    return gfsu.globalSize();
  }

  //! get dimension of space v
  typename GFSV::Traits::SizeType globalSizeV () const
  {
    return gfsv.globalSize();
  }

  //! get the trial grid function space
  const GFSU& trialGridFunctionSpace() const
  {
    return gfsu;
  }

  //! get the test grid function space
  const GFSV& testGridFunctionSpace() const
  {
    return gfsv;
  }

  typedef operator_applier<
    const MultiDomainGridOperatorSpace,
    data::ElementData,
    data::NeighborData,
    data::CouplingFunctionSpaces,
    data::IntersectionReference,
    data::SkeletonInvocationTracker,
    data::EnrichedCouplingInvocationTracker
    > operator_applier_all_data;


  /**\brief Construct global sparsity pattern from local description

     This function can be called by the Matrix to get the sparsity pattern.
     Assumes that the pattern is initially empty.
  */
  template<typename P>
  void fill_pattern (P& globalpattern) const
  {
    // make local function spaces
    typedef typename GFSU::LocalFunctionSpace LFSU;
    LFSU lfsu(gfsu);
    typedef typename GFSV::LocalFunctionSpace LFSV;
    LFSV lfsv(gfsv);

    operator_applier<
      const MultiDomainGridOperatorSpace,
      data::ElementData,
      data::NeighborData,
      data::CouplingFunctionSpaces
      > apply_operator(*this);
    apply_operator.setlfsu(lfsu);
    apply_operator.setlfsv(lfsv);

    for (ElementIterator it = gfsu.gridview().template begin<0>();
         it!=gfsu.gridview().template end<0>(); ++it)
      {
        // bind local function spaces to element
        lfsu.bind(*it);
        lfsv.bind(*it);
        apply_operator.setElement(*it);
        apply_operator.setElementSubDomains(gfsu.gridview().indexSet().subDomains(*it));

        apply_operator.template conditional<SubProblems,do_pattern_volume<> >(BuildVolumePattern<P>(globalpattern));

        // skeleton and boundary pattern
        if (!(any_child<MultiDomainGridOperatorSpace,SubProblems,do_pattern_skeleton<> >::value ||
              any_child<MultiDomainGridOperatorSpace,Couplings,do_pattern_coupling<> >::value ||
              any_child<MultiDomainGridOperatorSpace,Couplings,do_pattern_enriched_coupling<> >::value))
          continue;

        // local function spaces in neighbor
        LFSU lfsun(gfsu);
        LFSV lfsvn(gfsv);
        apply_operator.setlfsun(lfsun);
        apply_operator.setlfsvn(lfsvn);

        typedef typename GFSU::CouplingLocalFunctionSpace CouplingLFSU;
        CouplingLFSU couplinglfsu(gfsu);
        typedef typename GFSV::CouplingLocalFunctionSpace CouplingLFSV;
        CouplingLFSV couplinglfsv(gfsv);
        apply_operator.setcouplinglfsu(couplinglfsu);
        apply_operator.setcouplinglfsv(couplinglfsv);

        IntersectionIterator endit = gfsu.gridview().iend(*it);
        for (IntersectionIterator iit = gfsu.gridview().ibegin(*it); iit!=endit; ++iit)
          {
            // skip if there is no neighbor
            if (!iit->neighbor()) continue;

            // bind local function spaces to neighbor element
            lfsun.bind(*(iit->outside()));
            lfsvn.bind(*(iit->outside()));
            apply_operator.setNeighborSubDomains(gfsu.gridview().indexSet().subDomains(*(iit->outside())));

            // only bind the coupling local function spaces if anyone is going to use them
            if (any_child<MultiDomainGridOperatorSpace,Couplings,do_pattern_enriched_coupling<> >::value)
              {
                couplinglfsu.bind(*iit);
                couplinglfsv.bind(*iit);
              }

            // get pattern
            apply_operator.template conditional<SubProblems,do_pattern_skeleton<> >(BuildSkeletonPattern<P>(globalpattern));
            apply_operator.template conditional<Couplings,do_pattern_coupling<> >(BuildCouplingPattern<P>(globalpattern));
            apply_operator.template conditional<Couplings,do_pattern_enriched_coupling<> >(BuildEnrichedCouplingPattern<P>(globalpattern));

          }
      }
  }


  //! generic evaluation of residual
  /**
   * \param r residual (needs to be cleared before this method is called)
   */
  template<typename X, typename R>
  void residual (const X& x, R& r) const
  {
    // visit each face only once
    const int chunk=1<<28;
    int offset = 0;
    const typename GV::IndexSet& is=gfsu.gridview().indexSet();
    std::map<Dune::GeometryType,int> gtoffset;

    // make local function spaces
    typedef typename GFSU::LocalFunctionSpace LFSU;
    LFSU lfsu(gfsu);
    typedef typename GFSV::LocalFunctionSpace LFSV;
    LFSV lfsv(gfsv);

    operator_applier_all_data apply_operator(*this);
    apply_operator.setlfsu(lfsu);
    apply_operator.setlfsv(lfsv);

    // traverse grid view
    for (ElementIterator it = gfsu.gridview().template begin<0>();
         it!=gfsu.gridview().template end<0>(); ++it)
      {
        // assign offset for geometry type;
        if (gtoffset.find(it->type())==gtoffset.end())
          {
            gtoffset[it->type()] = offset;
            offset += chunk;
          }

        // compute unique id
        int id = is.index(*it)+gtoffset[it->type()];

        // skip ghost and overlap
        if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
          continue;

        // bind local function spaces to element
        lfsu.bind(*it);
        lfsv.bind(*it);
        apply_operator.setElement(*it);
        apply_operator.setElementSubDomains(is.subDomains(*it));

        // allocate local data container
        typedef std::vector<typename X::ElementType> XL;
        XL xl(lfsu.size());
        typedef std::vector<typename R::ElementType> RL;
        RL rl(lfsv.size(),0.0);

        // read coefficents
        lfsu.vread(x,xl);

        // volume evaluation
        apply_operator.template conditional<SubProblems,do_alpha_volume<> >(InvokeAlphaVolume<XL,RL>(xl,rl));
        apply_operator.template conditional<SubProblems,do_lambda_volume<> >(InvokeLambdaVolume<RL>(rl));

        // skip if no intersection iterator is needed
        if (any_child<MultiDomainGridOperatorSpace,SubProblems,do_alpha_skeleton<> >::value ||
            any_child<MultiDomainGridOperatorSpace,SubProblems,do_alpha_boundary<> >::value ||
            any_child<MultiDomainGridOperatorSpace,SubProblems,do_lambda_boundary<> >::value ||
            any_child<MultiDomainGridOperatorSpace,Couplings,do_alpha_coupling<> >::value ||
            any_child<MultiDomainGridOperatorSpace,Couplings,do_alpha_enriched_coupling<> >::value)
          {
            // local function spaces in neighbor
            LFSU lfsun(gfsu);
            LFSV lfsvn(gfsv);
            apply_operator.setlfsun(lfsun);
            apply_operator.setlfsvn(lfsvn);

            typedef typename GFSU::CouplingLocalFunctionSpace CouplingLFSU;
            CouplingLFSU couplinglfsu(gfsu);
            typedef typename GFSV::CouplingLocalFunctionSpace CouplingLFSV;
            CouplingLFSV couplinglfsv(gfsv);
            apply_operator.setcouplinglfsu(couplinglfsu);
            apply_operator.setcouplinglfsv(couplinglfsv);

            // traverse intersections
            unsigned int intersection_index = 0;
            IntersectionIterator endit = gfsu.gridview().iend(*it);
            for (IntersectionIterator iit = gfsu.gridview().ibegin(*it);
                 iit!=endit; ++iit, ++intersection_index)
              {
                apply_operator.setIntersection(*iit,intersection_index);
                // skeleton term
                if (iit->neighbor())
                  {
                    // assign offset for geometry type;
                    Dune::GeometryType gtn = iit->outside()->type();
                    if (gtoffset.find(gtn)==gtoffset.end())
                      {
                        gtoffset[gtn] = offset;
                        offset += chunk;
                      }

                    // compute unique id for neighbor
                    int idn = is.index(*(iit->outside()))+gtoffset[gtn];
                    lfsun.bind(*(iit->outside()));
                    lfsvn.bind(*(iit->outside()));
                    apply_operator.setNeighborSubDomains(gfsu.gridview().indexSet().subDomains(*(iit->outside())));

                    // allocate local data container
                    XL xn(lfsun.size());
                    RL rn(lfsvn.size(),0.0);

                    XL xcoupling;
                    RL rcoupling;

                    // only bind the coupling local function spaces if anyone is going to use them
                    if (any_child<MultiDomainGridOperatorSpace,Couplings,do_alpha_enriched_coupling<> >::value)
                      {
                        couplinglfsu.bind(*iit);
                        couplinglfsv.bind(*iit);
                        couplinglfsu.vread(x,xcoupling);
                        rcoupling.resize(couplinglfsv.size(),0.0);
                      }

                    // read coefficents
                    lfsun.vread(x,xn);

                    // unique vist of intersection
                    apply_operator.template conditional<SubProblems,do_alpha_skeleton_or_boundary<> >
                      (InvokeAlphaSkeletonOrBoundary<XL,RL>(xl,xn,rl,rn,
                                                            id > idn ||
                                                            (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                            )
                       );
                    apply_operator.template conditional<Couplings,do_alpha_coupling<> >
                      (InvokeAlphaCoupling<XL,RL>(xl,xn,rl,rn));
                    apply_operator.template conditional<Couplings,do_alpha_enriched_coupling<> >
                      (InvokeAlphaEnrichedCoupling<XL,RL>(xl,xn,xcoupling,rl,rn,rcoupling));
                    if (apply_operator.alphaSkeletonInvoked())
                      {
                        lfsvn.vadd(rn,r);
                        apply_operator.clearAlphaSkeletonInvoked();
                      }
                    if (apply_operator.alphaEnrichedCouplingInvoked())
                      {
                        couplinglfsv.vadd(rcoupling,r);
                        apply_operator.clearAlphaEnrichedCouplingInvoked();
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    apply_operator.template conditional<SubProblems,do_alpha_boundary<> >(InvokeAlphaBoundary<XL,RL>(xl,rl));
                    apply_operator.template conditional<SubProblems,do_lambda_boundary<> >(InvokeLambdaBoundary<RL>(rl));
                  }
              }
          }

        apply_operator.template conditional<SubProblems,do_alpha_volume_post_skeleton<> >(InvokeAlphaVolumePostSkeleton<XL,RL>(xl,rl));
        apply_operator.template conditional<SubProblems,do_lambda_volume_post_skeleton<> >(InvokeLambdaVolumePostSkeleton<RL>(rl));

        // accumulate result (note: r needs to be cleared outside)
        lfsv.vadd(rl,r);
      }

    // set residual to zero on constrained dofs
    Dune::PDELab::constrain_residual(*pconstraintsv,r);
  }

  //! generic application of Jacobian
  template<typename X, typename Y>
  void jacobian_apply (X& x, Y& y) const
  {
    // visit each face only once
    const int chunk=1<<28;
    int offset = 0;
    const typename GV::IndexSet& is=gfsu.gridview().indexSet();
    std::map<Dune::GeometryType,int> gtoffset;

    // make local function spaces
    typedef typename GFSU::LocalFunctionSpace LFSU;
    LFSU lfsu(gfsu);
    typedef typename GFSV::LocalFunctionSpace LFSV;
    LFSV lfsv(gfsv);

    operator_applier_all_data apply_operator(*this);
    apply_operator.setlfsu(lfsu);
    apply_operator.setlfsv(lfsv);

    // traverse grid view
    for (ElementIterator it = gfsu.gridview().template begin<0>();
         it!=gfsu.gridview().template end<0>(); ++it)
      {
        // assign offset for geometry type;
        if (gtoffset.find(it->type())==gtoffset.end())
          {
            gtoffset[it->type()] = offset;
            offset += chunk;
          }

        // compute unique id
        int id = is.index(*it)+gtoffset[it->type()];

        // skip ghost and overlap
        if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
          continue;

        // bind local function spaces to element
        lfsu.bind(*it);
        lfsv.bind(*it);
        apply_operator.setElement(*it);
        apply_operator.setElementSubDomains(is.subDomains(*it));

        // allocate local data container
        typedef std::vector<typename X::ElementType> XL;
        XL xl(lfsu.size());
        typedef std::vector<typename Y::ElementType> YL;
        YL yl(lfsv.size(),0.0);

        // read coefficents
        lfsu.vread(x,xl);

        // volume evaluation
        apply_operator.template conditional<SubProblems,do_alpha_volume<> >(InvokeJacobianApplyVolume<XL,YL>(xl,yl));

        // skeleton and boundary evaluation
        if (any_child<MultiDomainGridOperatorSpace,SubProblems,do_alpha_skeleton<> >::value ||
            any_child<MultiDomainGridOperatorSpace,SubProblems,do_alpha_boundary<> >::value ||
            any_child<MultiDomainGridOperatorSpace,Couplings,do_alpha_coupling<> >::value ||
            any_child<MultiDomainGridOperatorSpace,Couplings,do_alpha_enriched_coupling<> >::value)
          {
            // local function spaces in neighbor
            LFSU lfsun(gfsu);
            LFSV lfsvn(gfsv);
            apply_operator.setlfsun(lfsun);
            apply_operator.setlfsvn(lfsvn);

            typedef typename GFSU::CouplingLocalFunctionSpace CouplingLFSU;
            CouplingLFSU couplinglfsu(gfsu);
            typedef typename GFSV::CouplingLocalFunctionSpace CouplingLFSV;
            CouplingLFSV couplinglfsv(gfsv);
            apply_operator.setcouplinglfsu(couplinglfsu);
            apply_operator.setcouplinglfsv(couplinglfsv);

            unsigned int intersection_index = 0;
            IntersectionIterator endit = gfsu.gridview().iend(*it);
            for (IntersectionIterator iit = gfsu.gridview().ibegin(*it);
                 iit!=endit; ++iit, ++intersection_index)
              {
                apply_operator.setIntersection(*iit,intersection_index);
                // skeleton term
                if (iit->neighbor())
                  {
                    // assign offset for geometry type;
                    Dune::GeometryType gtn = iit->outside()->type();
                    if (gtoffset.find(gtn)==gtoffset.end())
                      {
                        gtoffset[gtn] = offset;
                        offset += chunk;
                      }

                    // compute unique id for neighbor
                    int idn = is.index(*(iit->outside()))+gtoffset[gtn];

                    lfsun.bind(*(iit->outside()));
                    lfsvn.bind(*(iit->outside()));
                    apply_operator.setNeighborSubDomains(gfsu.gridview().indexSet().subDomains(*(iit->outside())));

                    // allocate local data container
                    XL xn(lfsun.size());
                    YL yn(lfsvn.size(),0.0);

                    XL xcoupling;
                    YL ycoupling;

                    // only bind the coupling local function spaces if anyone is going to use them
                    if (any_child<MultiDomainGridOperatorSpace,Couplings,do_alpha_enriched_coupling<> >::value)
                      {
                        couplinglfsu.bind(*iit);
                        couplinglfsv.bind(*iit);
                        couplinglfsu.vread(x,xcoupling);
                        ycoupling.resize(couplinglfsv.size(),0.0);
                      }

                    // read coefficents
                    lfsun.vread(x,xn);

                    apply_operator.template conditional<SubProblems,do_alpha_skeleton_or_boundary<> >
                      (InvokeJacobianApplySkeletonOrBoundary<XL,YL>(xl,xn,yl,yn,
                                                                    id > idn ||
                                                                    (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                                    )
                       );
                    apply_operator.template conditional<Couplings,do_alpha_coupling<> >
                      (InvokeJacobianApplyCoupling<XL,YL>(xl,xn,yl,yn));
                    apply_operator.template conditional<Couplings,do_alpha_enriched_coupling<> >
                      (InvokeJacobianApplyEnrichedCoupling<XL,YL>(xl,xn,xcoupling,yl,yn,ycoupling));
                    if (apply_operator.alphaSkeletonInvoked())
                      {
                        lfsvn.vadd(yn,y);
                        apply_operator.clearAlphaSkeletonInvoked();
                      }
                    if (apply_operator.alphaEnrichedCouplingInvoked())
                      {
                        couplinglfsv.vadd(ycoupling,y);
                        apply_operator.clearAlphaEnrichedCouplingInvoked();
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    apply_operator.template conditional<SubProblems,do_alpha_boundary<> >(InvokeJacobianApplyBoundary<XL,YL>(xl,yl));
                  }
              }
          }

        apply_operator.template conditional<SubProblems,do_alpha_volume_post_skeleton<> >(InvokeJacobianApplyVolumePostSkeleton<XL,YL>(xl,yl));

        // accumulate result (note: r needs to be cleared outside)
        lfsv.vadd(yl,y);
      }

    // set residual to zero on constrained dofs
    Dune::PDELab::copy_constrained_dofs(*pconstraintsu,x,y);
  }

  //! generic assembly of Jacobian
  /**
   * \param x Where (in the space spanned by the dofs) to evaluate the Jacobian
   * \param a Jacobian (needs to be cleared before passed to this method)
   */
  template<typename X, typename A>
  void jacobian (const X& x, A& a) const
  {
    // visit each face only once
    const int chunk=1<<28;
    int offset = 0;
    const typename GV::IndexSet& is=gfsu.gridview().indexSet();
    std::map<Dune::GeometryType,int> gtoffset;

    // make local function spaces
    typedef typename GFSU::LocalFunctionSpace LFSU;
    LFSU lfsu(gfsu);
    typedef typename GFSV::LocalFunctionSpace LFSV;
    LFSV lfsv(gfsv);

    operator_applier_all_data apply_operator(*this);
    apply_operator.setlfsu(lfsu);
    apply_operator.setlfsv(lfsv);

    // traverse grid view
    for (ElementIterator it = gfsu.gridview().template begin<0>();
         it!=gfsu.gridview().template end<0>(); ++it)
      {
        // assign offset for geometry type;
        if (gtoffset.find(it->type())==gtoffset.end())
          {
            gtoffset[it->type()] = offset;
            offset += chunk;
          }

        // compute unique id
        const typename GV::IndexSet::IndexType id = is.index(*it)+gtoffset[it->type()];

        // skip ghost and overlap
        if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
          continue;

        // bind local function spaces to element
        lfsu.bind(*it);
        lfsv.bind(*it);
        apply_operator.setElement(*it);
        apply_operator.setElementSubDomains(is.subDomains(*it));

        // allocate local data container
        typedef std::vector<typename X::ElementType> XL;
        XL xl(lfsu.size());
        typedef LocalMatrix<typename A::ElementType> AL;
        AL al(lfsv.size(),lfsu.size(),0.0);

        // read coefficents
        lfsu.vread(x,xl);

        // volume evaluation
        apply_operator.template conditional<SubProblems,do_alpha_volume<> >(InvokeJacobianVolume<XL,AL>(xl,al));

        // skeleton and boundary evaluation
        if (any_child<MultiDomainGridOperatorSpace,SubProblems,do_alpha_skeleton<> >::value ||
            any_child<MultiDomainGridOperatorSpace,SubProblems,do_alpha_boundary<> >::value ||
            any_child<MultiDomainGridOperatorSpace,Couplings,do_alpha_coupling<> >::value ||
            any_child<MultiDomainGridOperatorSpace,Couplings,do_alpha_enriched_coupling<> >::value)
          {
            // local function spaces in neighbor
            LFSU lfsun(gfsu);
            LFSV lfsvn(gfsv);
            apply_operator.setlfsun(lfsun);
            apply_operator.setlfsvn(lfsvn);

            typedef typename GFSU::CouplingLocalFunctionSpace CouplingLFSU;
            CouplingLFSU couplinglfsu(gfsu);
            typedef typename GFSV::CouplingLocalFunctionSpace CouplingLFSV;
            CouplingLFSV couplinglfsv(gfsv);
            apply_operator.setcouplinglfsu(couplinglfsu);
            apply_operator.setcouplinglfsv(couplinglfsv);

            unsigned int intersection_index = 0;
            IntersectionIterator endit = gfsu.gridview().iend(*it);
            for (IntersectionIterator iit = gfsu.gridview().ibegin(*it);
                 iit!=endit; ++iit, ++intersection_index)
              {
                apply_operator.setIntersection(*iit,intersection_index);
                // skeleton term
                if (iit->neighbor())
                  {
                    // assign offset for geometry type;
                    Dune::GeometryType gtn = iit->outside()->type();
                    if (gtoffset.find(gtn)==gtoffset.end())
                      {
                        gtoffset[gtn] = offset;
                        offset += chunk;
                      }

                    // compute unique id for neighbor
                    const typename GV::IndexSet::IndexType idn = is.index(*(iit->outside()))+gtoffset[gtn];

                    // bind local function spaces to neighbor element
                    lfsun.bind(*(iit->outside()));
                    lfsvn.bind(*(iit->outside()));
                    apply_operator.setNeighborSubDomains(gfsu.gridview().indexSet().subDomains(*(iit->outside())));

                        // allocate local data container
                    XL xn(lfsun.size());
                    AL al_sn(lfsv.size() ,lfsun.size(),0.0);
                    AL al_ns(lfsvn.size(),lfsu.size() ,0.0);
                    AL al_nn(lfsvn.size(),lfsun.size(),0.0);

                    XL xcoupling;
                    AL al_cc;
                    AL al_sc;
                    AL al_cs;
                    AL al_nc;
                    AL al_cn;

                    // only bind the coupling local function spaces if anyone is going to use them
                    if (any_child<MultiDomainGridOperatorSpace,Couplings,do_alpha_enriched_coupling<> >::value)
                      {
                        couplinglfsu.bind(*iit);
                        couplinglfsv.bind(*iit);
                        couplinglfsu.vread(x,xcoupling);
                        al_cc.resize(couplinglfsv.size(),couplinglfsu.size());
                        al_cc = 0.0;
                        al_sc.resize(lfsv.size(),couplinglfsu.size());
                        al_sc = 0.0;
                        al_cs.resize(couplinglfsv.size(),lfsu.size());
                        al_cs = 0.0;
                        al_nc.resize(lfsvn.size(),couplinglfsu.size());
                        al_nc = 0.0;
                        al_cn.resize(couplinglfsv.size(),lfsun.size());
                      }

                    // read coefficents
                    lfsun.vread(x,xn);

                    apply_operator.template conditional<SubProblems,do_alpha_skeleton_or_boundary<> >
                      (InvokeJacobianSkeletonOrBoundary<XL,AL>(xl,xn,al,al_sn,al_ns,al_nn,
                                                                    id > idn ||
                                                                    (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                                    )
                       );
                    apply_operator.template conditional<Couplings,do_alpha_coupling<> >
                      (InvokeJacobianCoupling<XL,AL>(xl,xn,al,al_sn,al_ns,al_nn));
                    apply_operator.template conditional<Couplings,do_alpha_enriched_coupling<> >
                      (InvokeJacobianEnrichedCoupling<XL,AL>(xl,xn,xcoupling,al,al_nn,al_sc,al_cs,al_nc,al_cn,al_cc));
                    if (apply_operator.alphaSkeletonInvoked())
                      {
                        etadd(lfsv,lfsun,al_sn,a);
                        etadd(lfsvn,lfsu,al_ns,a);
                        etadd(lfsvn,lfsun,al_nn,a);
                        apply_operator.clearAlphaSkeletonInvoked();
                      }
                    if (apply_operator.alphaEnrichedCouplingInvoked())
                      {
                        etadd(couplinglfsv,couplinglfsu,al_cc,a);
                        etadd(lfsv,couplinglfsu,al_sc,a);
                        etadd(couplinglfsv,lfsu,al_cs,a);
                        etadd(lfsvn,couplinglfsu,al_nc,a);
                        etadd(couplinglfsv,lfsun,al_cn,a);
                        apply_operator.clearAlphaEnrichedCouplingInvoked();
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    apply_operator.template conditional<SubProblems,do_alpha_boundary<> >(InvokeJacobianBoundary<XL,AL>(xl,al));
                  }
              }
          }

        apply_operator.template conditional<SubProblems,do_alpha_volume_post_skeleton<> >(InvokeJacobianVolumePostSkeleton<XL,AL>(xl,al));

        // accumulate result (note: a needs to be cleared outside)
        etadd(lfsv,lfsu,al,a);
      }


    typedef typename CV::const_iterator global_row_iterator;
    for (global_row_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit)
      set_trivial_row(cit->first,cit->second,a);
  }

  /** \brief Transforms a vector \f$ \boldsymbol{x} \f$ from \f$
      V\f$ to \f$ V'\f$. If postrestrict == true then
      \f$\boldsymbol{R}^T_{\boldsymbol{\tilde U}', \boldsymbol{U}'}
      \boldsymbol{S}_{\boldsymbol{\tilde V}}\f$ is applied
      instead of the full transformation.  */
  template<typename X>
  void forwardtransform(X & x, const bool postrestrict = false)
  {
    typedef typename CV::const_iterator global_col_iterator;
    for (global_col_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit){
      typedef typename global_col_iterator::value_type::first_type GlobalIndex;
      const GlobalIndex & contributor = cit->first;

      typedef typename global_col_iterator::value_type::second_type ContributedMap;
      typedef typename ContributedMap::const_iterator global_row_iterator;
      const ContributedMap & contributed = cit->second;
      global_row_iterator it  = contributed.begin();
      global_row_iterator eit = contributed.end();

      for(;it!=eit;++it)
        x[it->first] += it->second * x[contributor];
    }

    if(postrestrict)
      for (global_col_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit)
        x[cit->first]=0.;
  }

  /** \brief Transforms a vector \f$ \boldsymbol{x} \f$ from \f$
      V'\f$ to \f$ V\f$. If prerestrict == true then
      \f$\boldsymbol{S}^T_{\boldsymbol{\tilde U}}\f$ is applied
      instead of the full transformation.  */
  template<typename X>
  void backtransform(X & x, const bool prerestrict = false)
  {
    typedef typename CV::const_iterator global_col_iterator;
    for (global_col_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit){
      typedef typename global_col_iterator::value_type::first_type GlobalIndex;
      const GlobalIndex & contributor = cit->first;

      typedef typename global_col_iterator::value_type::second_type ContributedMap;
      typedef typename ContributedMap::const_iterator global_row_iterator;
      const ContributedMap & contributed = cit->second;
      global_row_iterator it  = contributed.begin();
      global_row_iterator eit = contributed.end();

      if(prerestrict)
        x[contributor] = 0.;

      for(;it!=eit;++it)
        x[contributor] += it->second * x[it->first];
    }
  }

private:

  /** \brief read local stiffness matrix for entity */
  template<typename LFSV, typename LFSU, typename GC, typename T>
  void eread (const LFSV& lfsv, const LFSU& lfsu, const GC& globalcontainer,
              LocalMatrix<T>& localcontainer) const
  {
    for (int i=0; i<lfsv.size(); i++)
      for (int j=0; j<lfsu.size(); j++)
        localcontainer(i,j) = B::access(globalcontainer,lfsv.globalIndex(i),lfsu.globalIndex(j));
  }

  /** \brief write local stiffness matrix for entity */
  template<typename LFSV, typename LFSU, typename T, typename GC>
  void ewrite (const LFSV& lfsv, const LFSU& lfsu, const LocalMatrix<T>& localcontainer, GC& globalcontainer) const
  {
    for (int i=0; i<lfsv.size(); i++)
      for (int j=0; j<lfsu.size(); j++)
        B::access(globalcontainer,lfsv.globalIndex(i),lfsu.globalIndex(j)) = localcontainer(i,j);
  }

  /** \brief write local stiffness matrix for entity */
  template<typename LFSV, typename LFSU, typename T, typename GC>
  void eadd (const LFSV& lfsv, const LFSU& lfsu, const LocalMatrix<T>& localcontainer, GC& globalcontainer) const
  {
    for (size_t i=0; i<lfsv.size(); i++)
      for (size_t j=0; j<lfsu.size(); j++)
        B::access(globalcontainer,lfsv.globalIndex(i),lfsu.globalIndex(j)) += localcontainer(i,j);
  }

  /** \brief Add local matrix \f$m\f$ to global Jacobian \f$J\f$
      and apply constraints transformation. Hence we perform: \f$
      \boldsymbol{J} := \boldsymbol{J} + \boldsymbol{S}_{
      \boldsymbol{\tilde V}} m \boldsymbol{S}^T_{
      \boldsymbol{\tilde U}} \f$*/
  template<typename LFSV, typename LFSU, typename T, typename GC>
  void etadd (const LFSV& lfsv, const LFSU& lfsu, const LocalMatrix<T>& localcontainer, GC& globalcontainer) const
  {

    for (size_t i=0; i<lfsv.size(); i++)
      for (size_t j=0; j<lfsu.size(); j++){
        typename Traits::SizeType gi = lfsv.globalIndex(i);
        typename Traits::SizeType gj = lfsu.globalIndex(j);

        // Get global constraints containers for test and ansatz space
        const CV & cv = *pconstraintsv;
        const CU & cu = *pconstraintsu;

        typedef typename CV::const_iterator global_vcol_iterator;
        typedef typename global_vcol_iterator::value_type::second_type global_vrow_type;
        typedef typename global_vrow_type::const_iterator global_vrow_iterator;

        typedef typename CU::const_iterator global_ucol_iterator;
        typedef typename global_ucol_iterator::value_type::second_type global_urow_type;
        typedef typename global_urow_type::const_iterator global_urow_iterator;

        // Check whether the global indices are constrained indices
        global_vcol_iterator gvcit = cv.find(gi);
        global_ucol_iterator gucit = cu.find(gj);

        // Set constrained_v true if gi is constrained dof
        bool constrained_v(false);
        global_vrow_iterator gvrit;
        if(gvcit!=cv.end()){
          gvrit = gvcit->second.begin();
          constrained_v = true;
        }

        T vf = 1;
        do{
          // if gi is index of constrained dof
          if(constrained_v){

            if(gvrit == gvcit->second.end())
              break;

            // otherwise set gi to an index to a contributed dof
            // and set vf to the contribution weight
            gi = gvrit->first;
            vf = gvrit->second;
          }

          // Set constrained_u true if gj is constrained dof
          bool constrained_u(false);
          global_urow_iterator gurit;
          if(gucit!=cu.end()){
            gurit = gucit->second.begin();
            constrained_u = true;
            if(gurit == gucit->second.end()){
              T t = localcontainer(i,j) * vf;
              if(t != 0.0)                 // entry might not be present in the matrix
                B::access(globalcontainer,gi,gj) += t;
            }
          }

          T uf = 1;
          do{
            // if gj is index of constrained dof
            if(constrained_u){

              if(gurit == gucit->second.end())
                break;

              // otherwise set gj to an index to a contributed dof
              // and set uf to the contribution weight
              gj = gurit->first;
              uf = gurit->second;
            }

            // add weighted local entry to global matrix
            T t = localcontainer(i,j) * uf * vf;
            if (t != 0.0)                 // entry might not be present in the matrix
              B::access(globalcontainer,gi,gj) += t;

            if(constrained_u && gurit != gucit->second.end())
              ++gurit;
            else
              break;

          }while(true);

          if(constrained_v && gvrit != gvcit->second.end())
            ++gvrit;
          else
            break;

        }while(true);

      }
  }

  /** \brief Adding matrix entry to pattern with respect to the
      constraints contributions. This assembles the entries addressed
      by etadd(..). See the documentation there for more information
      about the matrix pattern. */
  template<typename GI, typename P>
  void add_entry(P & globalpattern, GI gi, GI gj) const
  {
    const CV & cv = *pconstraintsv;
    const CU & cu = *pconstraintsu;

    typedef typename CV::const_iterator global_vcol_iterator;
    typedef typename global_vcol_iterator::value_type::second_type global_vrow_type;
    typedef typename global_vrow_type::const_iterator global_vrow_iterator;

    typedef typename CU::const_iterator global_ucol_iterator;
    typedef typename global_ucol_iterator::value_type::second_type global_urow_type;
    typedef typename global_urow_type::const_iterator global_urow_iterator;

    global_vcol_iterator gvcit = cv.find(gi);
    global_ucol_iterator gucit = cu.find(gj);

    if(gi==gj)
      globalpattern.add_link(gi,gj);

    bool constrained_v(false);
    global_vrow_iterator gvrit;
    if(gvcit!=cv.end()){
      gvrit = gvcit->second.begin();
      constrained_v = true;
      if(gvrit == gvcit->second.end())
        globalpattern.add_link(gi,gj);
    }

    do{
      if(constrained_v){
        if(gvrit == gvcit->second.end())
          break;
        gi = gvrit->first;
      }

      bool constrained_u(false);
      global_urow_iterator gurit;
      if(gucit!=cu.end()){
        gurit = gucit->second.begin();
        constrained_u = true;
        if(gurit == gucit->second.end())
          globalpattern.add_link(gi,gj);
      }

      do{
        if(constrained_u){
          if(gurit == gucit->second.end())
            break;

          gj = gurit->first;
        }

        globalpattern.add_link(gi,gj);

        if(constrained_u && gurit != gucit->second.end())
          ++gurit;
        else
          break;

      }while(true);

      if(constrained_v && gvrit != gvcit->second.end())
        ++gvrit;
      else
        break;

    }while(true);

  }

  /** \brief insert dirichlet constraints for row and assemble
      T^T_U in constrained rows
  */
  template<typename GI, typename GC, typename CG>
  void set_trivial_row (GI i, const CG & cv_i, GC& globalcontainer) const
  {
    //std::cout << "clearing row " << i << std::endl;
    // set all entries in row i to zero
    B::clear_row(i,globalcontainer);

    // set diagonal element to 1
    B::access(globalcontainer,i,i) = 1;
  }


  const GFSU& gfsu;
  const GFSV& gfsv;
  const CU* pconstraintsu;
  const CV* pconstraintsv;
  CU emptyconstraintsu;
  CV emptyconstraintsv;
  bool nonoverlapping_mode;
};

//! \} group GridFunctionSpace


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_MULTIDOMAINGRIDOPERATORSPACE_HH
