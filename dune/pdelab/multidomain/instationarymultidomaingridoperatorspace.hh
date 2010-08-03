// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_MULTIDOMAIN_INSTATIONARYMULTIDOMAINGRIDOPERATORSPACE_HH
#define DUNE_MULTIDOMAIN_INSTATIONARYMULTIDOMAINGRIDOPERATORSPACE_HH

#include <map>
#include <tuple>
#include <limits>

#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>

#include <dune/pdelab/common/geometrywrapper.hh>
//#include"../gridfunctionspace/gridfunctionspace.hh"
#include <dune/pdelab/gridfunctionspace/constraints.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include <dune/pdelab/gridoperatorspace/instationarygridoperatorspace.hh>
#include <dune/pdelab/multidomain/operatorapplier.hh>

#include <dune/pdelab/multidomain/multidomaingridoperatorspaceutilities.hh>
#include <dune/pdelab/multidomain/multidomaingridoperatorspaceinvocationhelpers.hh>


namespace Dune {

namespace PDELab {

namespace MultiDomain {


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
template<typename TReal,
         typename R,
         typename GFSU,
         typename GFSV,
         typename B,
         typename... SubProblemsAndCouplings>
class InstationaryMultiDomainGridOperatorSpace : public VariadicCompositeNode<CopyStoragePolicy,SubProblemsAndCouplings...>
{

  template<typename,typename>
  friend struct BuildVolumePattern;

  template<typename,typename>
  friend struct BuildSkeletonPattern;

  template<typename,typename>
  friend struct BuildCouplingPattern;

  friend struct PreStep;

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
  typename SubProblem<k>::Type& subProblem() {
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

  template<std::size_t k>
  typename Coupling<k>::Type& coupling() {
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
    static const typename SubProblem<k>::Type& get(const InstationaryMultiDomainGridOperatorSpace& imgos)
    {
      dune_static_assert(k < N, "subProblem index too large");
      return imgos.template subProblem<k>();
    }

    template<std::size_t k>
    static typename SubProblem<k>::Type& get(InstationaryMultiDomainGridOperatorSpace& imgos)
    {
      dune_static_assert(k < N, "subProblem index too large");
      return imgos.template subProblem<k>();
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
    static const typename Coupling<k>::Type& get(const InstationaryMultiDomainGridOperatorSpace& imgos)
    {
      dune_static_assert(k < N, "coupling index too large");
      return imgos.template coupling<k>();
    }

    template<std::size_t k>
    static typename Coupling<k>::Type& get(InstationaryMultiDomainGridOperatorSpace& imgos)
    {
      dune_static_assert(k < N, "coupling index too large");
      return imgos.template coupling<k>();
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
    static const typename BaseT::template Child<k>::Type& get(const InstationaryMultiDomainGridOperatorSpace& imgos)
    {
      dune_static_assert(k < N, "child index too large");
      return imgos.template getChild<k>();
    }

    template<std::size_t k>
    static typename BaseT::template Child<k>::Type& get(InstationaryMultiDomainGridOperatorSpace& imgos)
    {
      dune_static_assert(k < N, "child index too large");
      return imgos.template getChild<k>();
    }

  };

public:

  typedef unsigned int StageType;

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
    typedef typename B::template Matrix<InstationaryMultiDomainGridOperatorSpace,E> Type;
  private:
    MatrixContainer () {}
  };

  //! construct GridOperatorSpace
  InstationaryMultiDomainGridOperatorSpace (const TimeSteppingParameterInterface<TReal>& method_,
                                            const GFSU& gfsu_,
                                            const GFSV& gfsv_,
                                            const CU& cu,
                                            const CV& cv,
                                            SubProblemsAndCouplings&... subProblems_)
    :  BaseT(subProblems_...),
       gfsu(gfsu_),
       gfsv(gfsv_),
       pconstraintsu(&cu),
       pconstraintsv(&cv),
       method(method_)
  {
  }

  //! construct GridOperatorSpace
  InstationaryMultiDomainGridOperatorSpace (const GFSU& gfsu_,
                                            const GFSV& gfsv_,
                                            const CU& cu,
                                            const CV& cv,
                                            SubProblemsAndCouplings&... subProblems_)
    :  BaseT(subProblems_...),
       gfsu(gfsu_),
       gfsv(gfsv_),
       pconstraintsu(&cu),
       pconstraintsv(&cv),
       method(&defaultMethod),
       r0(gfsv,0.0)
  {
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

  void setMethod(const TimeSteppingParameterInterface<TReal>& method_)
  {
    method = &method_;
  }

  void preStep(const TimeSteppingParameterInterface<TReal>& method_, TReal time_, TReal dt_)
  {
    setMethod(method_);
    preStep(time_,dt_);
  }

  void preStep(TReal time_, TReal dt_)
  {
    time = time_;
    dt = dt_;
    operator_applier<InstationaryMultiDomainGridOperatorSpace> apply_function(*this);
    apply_function.template apply<AllChildren>(PreStep());
  }

  void postStep()
  {
    operator_applier<InstationaryMultiDomainGridOperatorSpace> apply_function(*this);
    apply_function.template apply<AllChildren>(PostStep());
  }

  void postStage()
  {
    operator_applier<InstationaryMultiDomainGridOperatorSpace> apply_function(*this);
    apply_function.template apply<AllChildren>(PostStage());
  }

  TReal suggestTimestep(TReal dt) const
  {
    operator_applier<const InstationaryMultiDomainGridOperatorSpace> apply_function(*this);
    SuggestTimestep<TReal> suggested_dt(dt);
    apply_function.template apply<AllChildren>(suggested_dt);
    return suggested_dt.value();
  }

  typedef operator_applier<
    const InstationaryMultiDomainGridOperatorSpace,
    data::ElementData,
    data::NeighborData,
    data::IntersectionReference,
    data::SkeletonInvocationTracker
    > operator_applier_all_data;

  typedef operator_applier<
    InstationaryMultiDomainGridOperatorSpace,
    data::ElementData,
    data::NeighborData,
    data::IntersectionReference,
    data::SkeletonInvocationTracker
    > non_const_operator_applier_all_data;

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
      const InstationaryMultiDomainGridOperatorSpace,
      data::ElementData,
      data::NeighborData
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

        if (method->implicit())
          apply_operator.template conditional<SubProblems,do_pattern_volume<SpatialOperator> >(BuildVolumePattern<P,SpatialOperator>(globalpattern));
        apply_operator.template conditional<SubProblems,do_pattern_volume<TemporalOperator> >(BuildVolumePattern<P,TemporalOperator>(globalpattern));


        // skeleton and boundary pattern
        if (!(any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_pattern_skeleton<TemporalOperator> >::value ||
              any_child<InstationaryMultiDomainGridOperatorSpace,Couplings,do_pattern_coupling<CouplingOperator> >::value ||
              (any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_pattern_skeleton<SpatialOperator> >::value && method->implicit()))) continue;

        // local function spaces in neighbor
        LFSU lfsun(gfsu);
        LFSV lfsvn(gfsv);
        apply_operator.setlfsun(lfsun);
        apply_operator.setlfsvn(lfsvn);

        IntersectionIterator endit = gfsu.gridview().iend(*it);
        for (IntersectionIterator iit = gfsu.gridview().ibegin(*it); iit!=endit; ++iit)
          {
            // skip if there is no neighbor
            if (!iit->neighbor()) continue;

            // bind local function spaces to neighbor element
            lfsun.bind(*(iit->outside()));
            lfsvn.bind(*(iit->outside()));
            apply_operator.setNeighborSubDomains(gfsu.gridview().indexSet().subDomains(*(iit->outside())));

            // get pattern
            if (method->implicit())
              {
                apply_operator.template conditional<SubProblems,do_pattern_skeleton<SpatialOperator> >(BuildSkeletonPattern<P,SpatialOperator>(globalpattern));
                apply_operator.template conditional<Couplings,do_pattern_coupling<CouplingOperator> >(BuildCouplingPattern<P,CouplingOperator>(globalpattern));
              }
            apply_operator.template conditional<SubProblems,do_pattern_skeleton<TemporalOperator> >(BuildSkeletonPattern<P,TemporalOperator>(globalpattern));
          }
      }
  }


  template<typename X>
  void preStage(StageType stage_, const std::vector<X*>& x)
  {
    stage = stage_;
    if (x.size()!=stage)
      DUNE_THROW(Exception,"wrong number of solutions in InstationaryGridOperatorSpace");
    if (stage<1 || stage>method->s())
      DUNE_THROW(Exception,"invalid stage number in InstationaryGridOperatorSpace");

    // visit each face only once
    const int chunk=1<<28;
    int offset = 0;
    const typename GV::IndexSet& is=gfsu.gridview().indexSet();
    std::map<Dune::GeometryType,int> gtoffset;

    // extract coefficients of time stepping scheme
    std::vector<TReal> a(stage);
    for (size_t i=0; i<stage; ++i) a[i] = method->a(stage,i);
    std::vector<TReal> b(stage);
    for (size_t i=0; i<stage; ++i) b[i] = method->b(stage,i);
    std::vector<TReal> d(stage);
    for (size_t i=0; i<stage; ++i) d[i] = method->d(i);

    // make local function spaces
    typedef typename GFSU::LocalFunctionSpace LFSU;
    LFSU lfsu(gfsu);
    typedef typename GFSV::LocalFunctionSpace LFSV;
    LFSV lfsv(gfsv);

    non_const_operator_applier_all_data apply_operator(*this);
    apply_operator.setlfsu(lfsu);
    apply_operator.setlfsv(lfsv);

    // prepare local operators for stage
    apply_operator.template apply<AllChildren>(PreStage<TReal,StageType>(time+method->d(stage)*dt,stage));

    // clear constant part residual before assembling
    r0 = 0.0;

    const bool needsSkeleton =
      any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_alpha_skeleton<SpatialOperator> >::value ||
      any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_alpha_boundary<SpatialOperator> >::value ||
      any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_lambda_boundary<SpatialOperator> >::value ||
      any_child<InstationaryMultiDomainGridOperatorSpace,Couplings,do_alpha_coupling<CouplingOperator> >::value;

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


        // loop over all previous time steps
        for(StageType i = 0; i < stage; ++i)
          {
            // set time in local operators for evaluation
            apply_operator.template apply<AllChildren>(SetTime<TReal>(time+d[i]*dt));

            // allocate local data container
            typedef std::vector<typename X::ElementType> XL;
            XL xl(lfsu.size());
            typedef std::vector<typename R::ElementType> RL;
            RL rl_a(lfsv.size(),0.0);
            RL rl_m(lfsv.size(),0.0);

            // read coefficents
            lfsu.vread(*x[i],xl);

            const bool doM = std::abs(a[i]) > 1e-6;
            const bool doA = std::abs(b[i]) > 1e-6;

            // volume evaluation
            if (doA)
              {
                apply_operator.template conditional<SubProblems,do_alpha_volume<SpatialOperator> >(InvokeAlphaVolume<XL,RL,SpatialOperator>(xl,rl_a));
                apply_operator.template conditional<SubProblems,do_lambda_volume<SpatialOperator> >(InvokeLambdaVolume<RL,SpatialOperator>(rl_a));
              }
            if (doM)
              {
                apply_operator.template conditional<SubProblems,do_alpha_volume<TemporalOperator> >(InvokeAlphaVolume<XL,RL,TemporalOperator>(xl,rl_m));
              }

            // skip if no intersection iterator is needed
            if (needsSkeleton && doA)
              {
                // local function spaces in neighbor
                LFSU lfsun(gfsu);
                LFSV lfsvn(gfsv);
                apply_operator.setlfsun(lfsun);
                apply_operator.setlfsvn(lfsvn);

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

                        // read coefficents
                        lfsun.vread(*x[i],xn);

                        // unique vist of intersection
                        apply_operator.template conditional<SubProblems,do_alpha_skeleton_or_boundary<SpatialOperator> >
                          (InvokeAlphaSkeletonOrBoundary<XL,RL,SpatialOperator>(xl,xn,rl_a,rn,
                                                                                id > idn ||
                                                                                (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                                                )
                           );
                        apply_operator.template conditional<Couplings,do_alpha_coupling<CouplingOperator> >
                          (InvokeAlphaCoupling<XL,RL,CouplingOperator>(xl,xn,rl_a,rn));
                        if (apply_operator.alphaSkeletonInvoked())
                          {
                            for(typename RL::iterator it = rn.begin(); it != rn.end(); ++it)
                              (*it) *= b[i] * dt;
                            lfsvn.vadd(rn,r0);
                            apply_operator.clearAlphaSkeletonInvoked();
                          }
                      }

                    // boundary term
                    if (iit->boundary())
                      {
                        apply_operator.template conditional<SubProblems,do_alpha_boundary<SpatialOperator> >(InvokeAlphaBoundary<XL,RL,SpatialOperator>(xl,rl_a));
                        apply_operator.template conditional<SubProblems,do_lambda_boundary<SpatialOperator> >(InvokeLambdaBoundary<RL,SpatialOperator>(rl_a));
                      }
                  }
              }

            if (doA)
              {
                apply_operator.template conditional<SubProblems,do_alpha_volume_post_skeleton<SpatialOperator> >(InvokeAlphaVolumePostSkeleton<XL,RL,SpatialOperator>(xl,rl_a));
                apply_operator.template conditional<SubProblems,do_lambda_volume_post_skeleton<SpatialOperator> >(InvokeLambdaVolumePostSkeleton<RL,SpatialOperator>(rl_a));

                // accumulate result (note: r needs to be cleared outside)
                for (typename RL::iterator it = rl_a.begin(); it != rl_a.end(); ++it)
                  (*it) *= b[i] * dt;
                lfsv.vadd(rl_a,r0);
              }
            if (doM)
              {
                for (typename RL::iterator it = rl_m.begin(); it != rl_m.end(); ++it)
                  (*it) *= a[i];
                lfsv.vadd(rl_m,r0);
              }
          }
      }
  }

  template<typename X, typename A>
  void explicit_jacobian_residual(StageType stage_, const std::vector<X*>& x, A& mat, R& alpha, R& beta)
  {
    // process arguments
    stage = stage_;
    if (x.size()!=stage+1)
      DUNE_THROW(Exception,"wrong number of solutions in InstationaryGridOperatorSpace");
    if (stage<1 || stage>method->s())
      DUNE_THROW(Exception,"invalid stage number in InstationaryGridOperatorSpace");
    if (method->implicit())
      DUNE_THROW(Exception,"explicit mode called with implicit scheme");

    // visit each face only once
    const int chunk=1<<28;
    int offset = 0;
    const typename GV::IndexSet& is=gfsu.gridview().indexSet();
    std::map<Dune::GeometryType,int> gtoffset;

    // extract coefficients of time stepping scheme
    std::vector<TReal> a(stage);
    for (size_t i=0; i<stage; ++i) a[i] = method->a(stage,i);
    std::vector<TReal> b(stage);
    for (size_t i=0; i<stage; ++i) b[i] = method->b(stage,i);
    std::vector<TReal> d(stage);
    for (size_t i=0; i<stage; ++i) d[i] = method->d(i);
    const TReal b_rr = method->b(stage,stage);
    const TReal d_r = method->d(stage);

    // make local function spaces
    typedef typename GFSU::LocalFunctionSpace LFSU;
    LFSU lfsu(gfsu);
    typedef typename GFSV::LocalFunctionSpace LFSV;
    LFSV lfsv(gfsv);

    non_const_operator_applier_all_data apply_operator(*this);
    apply_operator.setlfsu(lfsu);
    apply_operator.setlfsv(lfsv);

    // prepare local operators for stage
    apply_operator.template apply<AllChildren>(PreStage<TReal,StageType>(time+method->d(stage)*dt,stage));

    // clear constant part residual before assembling
    r0 = 0.0;

    const bool needsSkeleton =
      any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_alpha_skeleton<SpatialOperator> >::value ||
      any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_alpha_boundary<SpatialOperator> >::value ||
      any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_lambda_boundary<SpatialOperator> >::value ||
      any_child<InstationaryMultiDomainGridOperatorSpace,Couplings,do_alpha_coupling<CouplingOperator> >::value;

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

        // residual part
        // loop over all previous time steps
        for(StageType i = 0; i < stage; ++i)
          {
            // set time in local operators for evaluation
            apply_operator.template apply<AllChildren>(SetTime<TReal>(time+d[i]*dt));

            // allocate local data container
            typedef std::vector<typename X::ElementType> XL;
            XL xl(lfsu.size());
            typedef std::vector<typename R::ElementType> RL;
            RL rl_a(lfsv.size(),0.0);
            RL rl_m(lfsv.size(),0.0);

            // read coefficents
            lfsu.vread(*x[i],xl);

            const bool doM = std::abs(a[i]) > 1e-6;
            const bool doA = std::abs(b[i]) > 1e-6;

            // volume evaluation
            if (doA)
              {
                apply_operator.template conditional<SubProblems,do_alpha_volume<SpatialOperator> >(InvokeAlphaVolume<XL,RL,SpatialOperator>(xl,rl_a));
                apply_operator.template conditional<SubProblems,do_lambda_volume<SpatialOperator> >(InvokeLambdaVolume<RL,SpatialOperator>(rl_a));
              }
            if (doM)
              {
                apply_operator.template conditional<SubProblems,do_alpha_volume<TemporalOperator> >(InvokeAlphaVolume<XL,RL,TemporalOperator>(xl,rl_m));
              }

            // skip if no intersection iterator is needed
            if (needsSkeleton && doA)
              {
                // local function spaces in neighbor
                LFSU lfsun(gfsu);
                LFSV lfsvn(gfsv);
                apply_operator.setlfsun(lfsun);
                apply_operator.setlfsvn(lfsvn);

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

                        // read coefficents
                        lfsun.vread(*x[i],xn);

                        // unique vist of intersection
                        apply_operator.template conditional<SubProblems,do_alpha_skeleton_or_boundary<SpatialOperator> >
                          (InvokeAlphaSkeletonOrBoundary<XL,RL,SpatialOperator>(xl,xn,rl_a,rn,
                                                                                id > idn ||
                                                                                (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                                                )
                           );
                        apply_operator.template conditional<Couplings,do_alpha_coupling<CouplingOperator> >
                          (InvokeAlphaCoupling<XL,RL,CouplingOperator>(xl,xn,rl_a,rn));
                        if (apply_operator.alphaSkeletonInvoked())
                          {
                            for(typename RL::iterator it = rn.begin(); it != rn.end(); ++it)
                              (*it) *= -b[i];
                            lfsvn.vadd(rn,beta);
                            apply_operator.clearAlphaSkeletonInvoked();
                          }
                      }

                    // boundary term
                    if (iit->boundary())
                      {
                        apply_operator.template conditional<SubProblems,do_alpha_boundary<SpatialOperator> >(InvokeAlphaBoundary<XL,RL,SpatialOperator>(xl,rl_a));
                        apply_operator.template conditional<SubProblems,do_lambda_boundary<SpatialOperator> >(InvokeLambdaBoundary<RL,SpatialOperator>(rl_a));
                      }
                  }
              }

            if (doA)
              {
                apply_operator.template conditional<SubProblems,do_alpha_volume_post_skeleton<SpatialOperator> >(InvokeAlphaVolumePostSkeleton<XL,RL,SpatialOperator>(xl,rl_a));
                apply_operator.template conditional<SubProblems,do_lambda_volume_post_skeleton<SpatialOperator> >(InvokeLambdaVolumePostSkeleton<RL,SpatialOperator>(rl_a));

                // accumulate result (note: r needs to be cleared outside)
                for (typename RL::iterator it = rl_a.begin(); it != rl_a.end(); ++it)
                  (*it) *= -b[i];
                lfsv.vadd(rl_a,beta);
              }
            if (doM)
              {
                for (typename RL::iterator it = rl_m.begin(); it != rl_m.end(); ++it)
                  (*it) *= -a[i];
                lfsv.vadd(rl_m,alpha);
              }
          }

        // Jacobian part
        // Note:
        // - we are explicit; there is no spatial part here
        // - temporal part has only alpha_volume

        // allocate local data container
        typedef std::vector<typename X::ElementType> XL;
        XL xl(lfsu.size());
        typedef LocalMatrix<typename A::ElementType> AL;
        AL ml(lfsv.size(),lfsu.size(),0.0);

        // set time in local operator for evaluation
        apply_operator.template apply<AllChildren>(SetTime<TReal>(time+d_r*dt));

        // read coefficents; this is only a dummy since Jacobian should not depend on solution !
        // but of course it is required to give this parameter
        lfsu.vread(*x[stage],xl);

        // compute local jacobian
        apply_operator.template conditional<SubProblems,do_alpha_volume<TemporalOperator> >(InvokeJacobianVolume<XL,AL,TemporalOperator>(xl,ml));

        // accumulate to global matrix
        etadd(lfsv,lfsu,ml,mat);
      }

    // set trivial conditions for constrained degrees of freedom
    typedef typename CV::const_iterator global_row_iterator;
    for (global_row_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit)
      set_trivial_row(cit->first,cit->second,mat);

    // set residual to zero on constrained dofs of spatial part (which is scaled by dt)
    Dune::PDELab::constrain_residual(*pconstraintsv,beta);

    // copy solution on constrained dofs from solution of stage to temporal part (which is not scaled)
    // this makes the boundary conditions appear in the solution !
    Dune::PDELab::copy_constrained_dofs(*pconstraintsu,*x[stage],alpha);
  }

  //! generic evaluation of residual
  /**
   * \param r residual (needs to be cleared before this method is called)
   */
  template<typename X>
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

    const TReal b_rr = method->b(stage,stage);
    const TReal d_r = method->d(stage);
    const bool implicit = method->implicit();

    // set time in local operators for evaluation
    apply_operator.template apply<AllChildren>(SetTime<TReal>(time+d_r*dt));

    // copy constant part of residual
    r = r0;

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
        RL rl_a(lfsv.size(),0.0);
        RL rl_m(lfsv.size(),0.0);

        // read coefficents
        lfsu.vread(x,xl);

        // volume evaluation
        if (implicit)
          {
            apply_operator.template conditional<SubProblems,do_alpha_volume<SpatialOperator> >(InvokeAlphaVolume<XL,RL,SpatialOperator>(xl,rl_a));
            apply_operator.template conditional<SubProblems,do_lambda_volume<SpatialOperator> >(InvokeLambdaVolume<RL,SpatialOperator>(rl_a));
          }
        apply_operator.template conditional<SubProblems,do_alpha_volume<TemporalOperator> >(InvokeAlphaVolume<XL,RL,TemporalOperator>(xl,rl_m));


        // skip if no intersection iterator is needed
        if (implicit &&
            (any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_alpha_skeleton<SpatialOperator> >::value ||
             any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_alpha_boundary<SpatialOperator> >::value ||
             any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_lambda_boundary<SpatialOperator> >::value ||
             any_child<InstationaryMultiDomainGridOperatorSpace,Couplings,do_alpha_coupling<CouplingOperator> >::value)
            )
          {
            // local function spaces in neighbor
            LFSU lfsun(gfsu);
            LFSV lfsvn(gfsv);
            apply_operator.setlfsun(lfsun);
            apply_operator.setlfsvn(lfsvn);

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

                    // read coefficents
                    lfsun.vread(x,xn);

                    // unique vist of intersection
                    apply_operator.template conditional<SubProblems,do_alpha_skeleton_or_boundary<SpatialOperator> >
                      (InvokeAlphaSkeletonOrBoundary<XL,RL,SpatialOperator>(xl,xn,rl_a,rn,
                                                                            id > idn ||
                                                                            (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                                            )
                       );
                    apply_operator.template conditional<Couplings,do_alpha_coupling<CouplingOperator> >
                      (InvokeAlphaCoupling<XL,RL,CouplingOperator>(xl,xn,rl_a,rn));
                    if (apply_operator.alphaSkeletonInvoked())
                      {
                        for(typename RL::iterator it = rn.begin(); it != rn.end(); ++it)
                          (*it) *= b_rr * dt;
                        lfsvn.vadd(rn,r);
                        apply_operator.clearAlphaSkeletonInvoked();
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    apply_operator.template conditional<SubProblems,do_alpha_boundary<SpatialOperator> >(InvokeAlphaBoundary<XL,RL,SpatialOperator>(xl,rl_a));
                    apply_operator.template conditional<SubProblems,do_lambda_boundary<SpatialOperator> >(InvokeLambdaBoundary<RL,SpatialOperator>(rl_a));
                  }
              }
          }

        if (implicit)
          {
            apply_operator.template conditional<SubProblems,do_alpha_volume_post_skeleton<SpatialOperator> >(InvokeAlphaVolumePostSkeleton<XL,RL,SpatialOperator>(xl,rl_a));
            apply_operator.template conditional<SubProblems,do_lambda_volume_post_skeleton<SpatialOperator> >(InvokeLambdaVolumePostSkeleton<RL,SpatialOperator>(rl_a));

            // accumulate result (note: r needs to be cleared outside)
            for (typename RL::iterator it = rl_a.begin(); it != rl_a.end(); ++it)
                   (*it) *= b_rr * dt;
            lfsv.vadd(rl_a,r);
          }

        // scheme is normalized !
        lfsv.vadd(rl_m,r);
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

    const TReal b_rr = method->b(stage,stage);
    const TReal d_r = method->d(stage);
    const bool implicit = method->implicit();

    // set time in local operators for evaluation
    apply_operator.template apply<AllChildren>(SetTime<TReal>(time+d_r*dt));

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
        YL yl_a(lfsv.size(),0.0);
        YL yl_m(lfsv.size(),0.0);


        // read coefficents
        lfsu.vread(x,xl);

        // volume evaluation
        if (implicit)
          apply_operator.template conditional<SubProblems,do_alpha_volume<SpatialOperator> >(InvokeJacobianApplyVolume<XL,YL,SpatialOperator>(xl,yl_a));
        apply_operator.template conditional<SubProblems,do_alpha_volume<TemporalOperator> >(InvokeJacobianApplyVolume<XL,YL,TemporalOperator>(xl,yl_m));

        // skeleton and boundary evaluation
        if (implicit && (any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_alpha_skeleton<SpatialOperator> >::value ||
                         any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_alpha_boundary<SpatialOperator> >::value ||
                         any_child<InstationaryMultiDomainGridOperatorSpace,Couplings,do_alpha_coupling<CouplingOperator> >::value))
          {
            // local function spaces in neighbor
            LFSU lfsun(gfsu);
            LFSV lfsvn(gfsv);
            apply_operator.setlfsun(lfsun);
            apply_operator.setlfsvn(lfsvn);

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
                    YL yn_a(lfsvn.size(),0.0);

                    // read coefficents
                    lfsun.vread(x,xn);

                    apply_operator.template conditional<SubProblems,do_alpha_skeleton_or_boundary<SpatialOperator> >
                      (InvokeJacobianApplySkeletonOrBoundary<XL,YL,SpatialOperator>(xl,xn,yl_a,yn_a,
                                                                                    id > idn ||
                                                                                    (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                                                    )
                       );
                    apply_operator.template conditional<Couplings,do_alpha_coupling<CouplingOperator> >
                      (InvokeJacobianCoupling<XL,YL,CouplingOperator>(xl,xn,yl_a,yn_a));
                    if (apply_operator.alphaSkeletonInvoked())
                      {
                        for(typename YL::iterator it = yn_a.begin(); it != yn_a.end(); ++it)
                          (*it) *= b_rr*dt;
                        lfsvn.vadd(yn_a,y);
                        apply_operator.clearAlphaSkeletonInvoked();
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    apply_operator.template conditional<SubProblems,do_alpha_boundary<SpatialOperator> >(InvokeJacobianApplyBoundary<XL,YL,SpatialOperator>(xl,yl_a));
                  }
              }
          }

        if (implicit)
          {
            apply_operator.template conditional<SubProblems,do_alpha_volume_post_skeleton<SpatialOperator> >(InvokeJacobianApplyVolumePostSkeleton<XL,YL,SpatialOperator>(xl,yl_a));
            for(typename YL::iterator it = yl_a.begin(); it != yl_a.end(); ++it)
              (*it) *= b_rr*dt;
            // accumulate result (note: r needs to be cleared outside)
            lfsv.vadd(yl_a,y);
          }
        lfsv.vadd(yl_m,y);
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

    const TReal b_rr = method->b(stage,stage);
    const TReal d_r = method->d(stage);
    const bool implicit = method->implicit();

    // set time in local operators for evaluation
    apply_operator.template apply<AllChildren>(SetTime<TReal>(time+d_r*dt));

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
        AL ml(lfsv.size(),lfsu.size(),0.0);

        // read coefficents
        lfsu.vread(x,xl);

        // volume evaluation
        if (implicit)
          apply_operator.template conditional<SubProblems,do_alpha_volume<SpatialOperator> >(InvokeJacobianVolume<XL,AL,SpatialOperator>(xl,al));
        apply_operator.template conditional<SubProblems,do_alpha_volume<TemporalOperator> >(InvokeJacobianVolume<XL,AL,TemporalOperator>(xl,ml));


        // skeleton and boundary evaluation
        if (implicit && (any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_alpha_skeleton<SpatialOperator> >::value ||
                         any_child<InstationaryMultiDomainGridOperatorSpace,SubProblems,do_alpha_boundary<SpatialOperator> >::value ||
                         any_child<InstationaryMultiDomainGridOperatorSpace,Couplings,do_alpha_coupling<CouplingOperator> >::value))
          {
            // local function spaces in neighbor
            LFSU lfsun(gfsu);
            LFSV lfsvn(gfsv);
            apply_operator.setlfsun(lfsun);
            apply_operator.setlfsvn(lfsvn);

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

                    // read coefficents
                    lfsun.vread(x,xn);

                    apply_operator.template conditional<SubProblems,do_alpha_skeleton_or_boundary<SpatialOperator> >
                      (InvokeJacobianSkeletonOrBoundary<XL,AL,SpatialOperator>(xl,xn,al,al_sn,al_ns,al_nn,
                                                                               id > idn ||
                                                                               (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                                               )
                       );
                    apply_operator.template conditional<Couplings,do_alpha_coupling<CouplingOperator> >
                      (InvokeJacobianCoupling<XL,AL,CouplingOperator>(xl,xn,al,al_sn,al_ns,al_nn));
                    if (apply_operator.alphaSkeletonInvoked())
                      {
                        al_sn *= b_rr*dt;
                        etadd(lfsv,lfsun,al_sn,a);
                        al_ns *= b_rr*dt;
                        etadd(lfsvn,lfsu,al_ns,a);
                        al_nn *= b_rr*dt;
                        etadd(lfsvn,lfsun,al_nn,a);
                        apply_operator.clearAlphaSkeletonInvoked();
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    apply_operator.template conditional<SubProblems,do_alpha_boundary<SpatialOperator> >(InvokeJacobianBoundary<XL,AL,SpatialOperator>(xl,al));
                  }
              }
          }

        if (implicit) {
          apply_operator.template conditional<SubProblems,do_alpha_volume_post_skeleton<SpatialOperator> >(InvokeJacobianVolumePostSkeleton<XL,AL,SpatialOperator>(xl,al));
          al *= b_rr*dt;
          etadd(lfsv,lfsu,al,a);
        }

        // accumulate result (note: a needs to be cleared outside)
        etadd(lfsv,lfsu,ml,a);
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
  const TimeSteppingParameterInterface<TReal>* method;
  TReal time, dt;
  StageType stage;
  R r0;
  ImplicitEulerParameter<TReal> defaultMethod;
};

//! \} group GridFunctionSpace


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_INSTATIONARYMULTIDOMAINGRIDOPERATORSPACE_HH
