// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_MULTIDOMAIN_MULTIDOMAINGRIDOPERATORSPACE_HH
#define DUNE_MULTIDOMAIN_MULTIDOMAINGRIDOPERATORSPACE_HH

#include<map>
#include<tuple>

#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>

//#include"../common/geometrywrapper.hh"
//#include"../gridfunctionspace/gridfunctionspace.hh"
//#include"../gridfunctionspace/constraints.hh"
//#include"localmatrix.hh"
//#include"gridoperatorspaceutilities.hh"

#include <dune/pdelab/multidomain/typemap.hh>


namespace Dune {

namespace PDELab {

namespace MultiDomain {

//namespace {

template<typename T>
struct tag
{
  static const bool value = false;
};

template<>
struct tag<int>
{
  static const bool value = true;
};

template<typename T, std::size_t i, std::size_t j>
struct SubProblemEntry
{
  typedef T type;
  typedef T key;
  static const std::size_t globalPos = i;
  static const std::size_t localPos = j;
};

template<std::size_t i, std::size_t j, typename... ProblemsAndCouplings>
struct extract_problems_helper;

template<std::size_t i, std::size_t j, bool is_subproblem, typename T, typename... ProblemsAndCouplings>
struct extract_problem_switch;

template<std::size_t i, std::size_t j, typename T, typename... ProblemsAndCouplings>
struct extract_problem_switch<i,j,false,T,ProblemsAndCouplings...>
{
  template<typename... Problems>
  struct result
  {
    typedef typename extract_problems_helper<i+1,j,ProblemsAndCouplings...>::template result<Problems...>::type type;
    typedef typename extract_problems_helper<i+1,j,ProblemsAndCouplings...>::template result<Problems...>::map_type map_type;
  };
};

template<std::size_t i, std::size_t j, typename T, typename... ProblemsAndCouplings>
struct extract_problem_switch<i,j,true,T,ProblemsAndCouplings...>
{
  template<typename... Problems>
  struct result
  {
    typedef typename extract_problems_helper<i+1,j+1,ProblemsAndCouplings...>::template result<Problems...,SubProblemEntry<T,i,j> >::type type;
    typedef typename extract_problems_helper<i+1,j+1,ProblemsAndCouplings...>::template result<Problems...,SubProblemEntry<T,i,j> >::map_type map_type;
  };
};

template<std::size_t i, std::size_t j, typename T, typename... ProblemsAndCouplings>
struct extract_problems_helper<i,j,T,ProblemsAndCouplings...>
{
  template<typename... Problems>
  struct result
  {
    typedef typename extract_problem_switch<i,j,tag<T>::value,T,ProblemsAndCouplings...>::template result<Problems...>::type type;
    typedef typename extract_problem_switch<i,j,tag<T>::value,T,ProblemsAndCouplings...>::template result<Problems...>::map_type map_type;
  };
};

template<std::size_t i, std::size_t j>
struct extract_problems_helper<i,j>
{
  template<typename... Problems>
  struct result
  {
    typedef std::tuple<Problems...> type;
    typedef embedded_key_map<Problems...> map_type;
  };
};

template<typename... ProblemsAndCouplings>
struct extract_problems
{
  typedef typename extract_problems_helper<0,0,ProblemsAndCouplings...>::template result<>::type type;
  typedef typename extract_problems_helper<0,0,ProblemsAndCouplings...>::template result<>::map_type map_type;
};

template<typename T1, typename T2>
struct lookup_subproblem
{
  typedef T1 type;
};

template<typename T, std::size_t i, std::size_t j, typename SubProblems>
struct CouplingEntry
{
  typedef T type;
  static const std::size_t globalPos = i;
  static const std::size_t localPos = j;
  typedef typename get_map_entry<typename T::Traits::FirstSubProblem,SubProblems>::type FirstSubProblem;
  typedef typename get_map_entry<typename T::Traits::SecondSubProblem,SubProblems>::type SecondSubProblem;

};

template<typename SubProblems, std::size_t i, std::size_t j, typename... ProblemsAndCouplings>
struct extract_couplings_helper;

template<typename SubProblems, std::size_t i, std::size_t j, bool is_coupling, typename T, typename... ProblemsAndCouplings>
struct extract_coupling_switch;

template<typename SubProblems, std::size_t i, std::size_t j, typename T, typename... ProblemsAndCouplings>
struct extract_coupling_switch<SubProblems,i,j,false,T,ProblemsAndCouplings...>
{
  template<typename... Couplings>
  struct result
  {
    typedef typename extract_couplings_helper<SubProblems,i+1,j,ProblemsAndCouplings...>::template result<Couplings...>::type type;
  };
};

template<typename SubProblems, std::size_t i, std::size_t j, typename T, typename... ProblemsAndCouplings>
struct extract_coupling_switch<SubProblems,i,j,true,T,ProblemsAndCouplings...>
{
  template<typename... Couplings>
  struct result
  {
    typedef typename extract_couplings_helper<SubProblems,i+1,j+1,ProblemsAndCouplings...>::template result<Couplings...,CouplingEntry<T,i,j,SubProblems> >::type type;
  };
};

template<typename SubProblems, std::size_t i, std::size_t j, typename T, typename... ProblemsAndCouplings>
struct extract_couplings_helper<SubProblems,i,j,T,ProblemsAndCouplings...>
{
  template<typename... Couplings>
  struct result
  {
    typedef typename extract_coupling_switch<SubProblems,i,j,!tag<T>::value,T,ProblemsAndCouplings...>::template result<Couplings...>::type type;
  };
};

template<typename SubProblems,std::size_t i, std::size_t j>
struct extract_couplings_helper<SubProblems,i,j>
{
  template<typename... Couplings>
  struct result
  {
    typedef std::tuple<Couplings...> type;
  };
};

template<typename... ProblemsAndCouplings>
struct extract_couplings
{
  typedef typename extract_couplings_helper<typename extract_problems<ProblemsAndCouplings...>::map_type,
                                            0,
                                            0,
                                            ProblemsAndCouplings...>::template result<>::type type;
};


//} // anonymous namespace

template<typename MDGOS, template<typename> class Condition, template<bool,bool> class BooleanOp, bool start_value, std::size_t i, std::size_t n>
struct child_condition
{
  static const bool value = BooleanOp<Condition<typename MDGOS::template Child<i>::Type>::value,
                                      child_condition<MDGOS,Condition,BooleanOp,start_value,i+1,n>::value
                                     >::value;
};

template<typename MDGOS, template<typename> class Condition, template<bool,bool> class BooleanOp, bool start_value, std::size_t n>
struct child_condition<MDGOS,Condition,BooleanOp,start_value,n,n>
{
  static const bool value = start_value;
};

template<bool a, bool b>
struct and_
{
  static const bool value = a && b;
};

template<bool a, bool b>
struct or_
{
  static const bool value = a || b;
};

template<typename MDGOS, template<typename> class Condition>
struct all_childs : public child_condition<MDGOS,Condition,and_,true,0,MDGOS::CHILDREN> {};

template<typename MDGOS, template<typename> class Condition>
struct any_child : public child_condition<MDGOS,Condition,or_,false,0,MDGOS::CHILDREN> {};

template<typename Applier, typename Operator, std::size_t i, std::size_t n>
struct apply_operator_helper
{
  static void apply(Applier& applier, Operator& op)
  {
    op(applier,mdgos.template getChild<i>());
    apply_operator_helper<Applier,Operator,i+1,n>::apply(applier,op);
  }
};

// end of recursion
template<typename Applier, typename Operator, std::size_t n>
struct apply_operator_helper<Applier, Operator, n,n>
{
  static void apply(Applier& applier, Operator& op)
  {
  }
};


template<typename Applier, typename Operator, bool do_apply>
struct conditional_apply;

template<typename Applier, typename Operator>
struct conditional_apply<Applier,Operator,true>
{
  static void apply(Applier& applier, Operator& op)
  {
    op(applier,mdgos.template getChild<i>());
  }
};

template<typename Applier, typename Operator>
struct conditional_apply<Applier,Operator,false>
{
  static void apply(Applier& applier, Operator& op)
  {
  }
};

template<typename Applier, template<typename> class Condition, typename Operator, std::size_t i, std::size_t n>
struct conditional_apply_operator_helper
{
  static void apply(Applier& applier, Operator& op)
  {
    conditional_apply<Applier,Operator,Condition<typename MDGOS::template Child<i>::Type>::value>::apply(applier,op);
    apply_operator_helper<Applier,Operator,i+1,n>::apply(applier,op);
  }
};

// end of recursion
template<typename Applier, template<typename> class Condition, typename Operator, std::size_t n>
struct conditional_apply_operator_helper<Applier, Condition, Operator, n,n>
{
  static void apply(Applier& applier, Operator& op)
  {
  }
};


template<typename MDGOS, typename LFSU, typename LFSV>
class operator_applier
{

  typedef typename MDGOS::Traits::TrialGridFunctionSpace::Traits::Grid Grid;
  typedef typename Grid::template Codim<0>::SubDomainSet ElementSubDomainSet;

public:

  operator_applier(MDGOS& mdgos) :
    _mdgos(mdgos),
    _alphaSkeletonInvoked(false)
  {}

  template<typename Operator>
  void operator()(Operator op)
  {
    apply_operator_helper<MDGOS,Operator,0,MDGOS::CHILDREN>::apply(_mdgos,op);
  }

  template<template<typename> class Condition, typename Operator)
  void conditional_apply(Operator op)
  {
    conditional_apply_operator_helper<MDGOS,Condition,Operator,0,MDGOS::CHILDREN>::apply(_mdgos,op);
  }

  MDGOS& gos() {
    return _mdgos;
  }

  LFSU& lfsu() {
    return *_plfsu;
  }

  void set(LFSU& lfsu) {
    _plfsu = &lfsu;
  }

  LFSV& lfsv() {
    return *_plfsv;
  }

  void set(LFSV& lfsv) {
    _plfsv = &lfsv;
  }

  LFSU& lfsun() {
    return *_plfsun;
  }

  void setn(LFSU& lfsun) {
    _plfsun = &lfsun;
  }

  LFSV& lfsvn() {
    return *_plfsvn;
  }

  void setn(LFSV& lfsvn) {
    _plfsvn = &lfsvn;
  }

  ElementSubDomainSet& elementSubDomains() {
    return _elementSet;
  }

  void setElementSubDomains(const ElementSubDomainSet& elementSet) {
    _elementSet = elementSet;
  }

  ElementSubDomainSet& neighborSubDomains() {
    return _neighborSet;
  }

  void setNeighborSubDomains(const NeighborSubDomainSet& intersectionSet) {
    _neighborSet = neighborSet;
  }

  void setAlphaSkeletonInvoked() {
    _alphaSkeletonInvoked = true;
  }

  void clearAlphaSkeletonInvoked() {
    _alphaSkeletonInvoked = false;
  }

  bool alphaSkeletonInvoked() const {
    return _alphaSkeletonInvoked;
  }


private:

  MDGOS& _mdgos;
  LFSU* _plfsu;
  LFSV* _plfsv;
  LFSU* _plfsun;
  LFSV* _plfsvn;
  ElementSubdomainSet _elementSet;
  ElementSubDomainSet _neighborSet;
  bool _alphaSkeletonInvoked;

};

template<typename T>
struct do_pattern_skeleton
{
  static const bool value = T::Traits::LocalOperator::doPatternSkeleton;
};

template<typename T>
struct do_pattern_volume
{
  static const bool value = T::Traits::LocalOperator::doPatternVolume;
};

template<typename T>
struct do_alpha_volume
{
  static const bool value = T::Traits::LocalOperator::doAlphaVolume;
};

template<typename T>
struct do_alpha_skeleton
{
  static const bool value = T::Traits::LocalOperator::doAlphaSkeleton;
};

template<typename T>
struct do_alpha_boundary
{
  static const bool value = T::Traits::LocalOperator::doAlphaBoundary;
};

template<typename T>
struct do_alpha_skeleton_or_boundary
{
  static const bool value = T::Traits::LocalOperator::doAlphaSkeleton || T::Traits::LocalOperator::doAlphaBoundary;
};

template<typename T>
struct do_lambda_volume
{
  static const bool value = T::Traits::LocalOperator::doLambdaVolume;
};

template<typename T>
struct do_lambda_boundary
{
  static const bool value = T::Traits::LocalOperator::doLambdaBoundary;
};

#if 0

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
class MultiDomainGridOperatorSpace : public VariadicCompositeNode<SubProblemsAndCouplings...>
{

  typedef typename extract_problems<SubProblemsAndCouplings...>::type SubProblemList;
  typedef typename extract_problems<SubProblemsAndCouplings...>::map_type SubProblemMap;
  typedef typename extract_couplings<SubProblemsAndCouplings...>::type CouplingList;

  static const std::size_t subProblemCount = std::tuple_size<SubProblemList>::value;
  static const std::size_t couplingCount = std::tuple_size<SubProblemList>::value;

  template<std::size_t k>
  struct SubProblem {
    typedef typename std::tuple_element<k,SubProblemList>::type::type Type;
  };

  template<std::size_t k>
  const SubProblem<k>::Type& subProblem() const {
    return this->getChild<std::tuple_element<k,SubProblemList>::type::globalIndex>();
  }

  template<std::size_t k>
  struct Coupling {
    typedef typename std::tuple_element<k,CouplingList>::type::type Type;
  };

  template<std::size_t k>
  const Coupling<k>::Type& coupling() const {
    return this->getChild<std::tuple_element<k,CouplingList>::type::globalIndex>();
  }

  // extract useful types
  typedef typename GFSU::Traits::GridViewType GV;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;

  template<typename P>
  struct BuildVolumePattern
  {
    BuildVolumePattern(P& gp) : globalpattern(gp) {}

    template<typename Data, typename Child>
    void operator()(Data& data, Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSU lfsu(data.lfsu());
      LFSV lfsv(data.lfsv());
      LocalSparsityPattern localpattern;
      child.localOperator().pattern_volume(lfsu,lfsv,localpattern);

      // translate local to global indices and add to global pattern
      for (size_t k=0; k<localpattern.size(); ++k)
        add_entry(globalpattern,
                  lfsv.globalIndex(localpattern[k].i()),
                      lfsu.globalIndex(localpattern[k].j())
                  );
    }

    P& globalpattern;
  };

  template<typename P>
  struct BuildSkeletonPattern
  {
    BuildSkeletonPattern(P& gp) : globalpattern(gp) {}

    template<typename Data, typename Child>
    void operator()(Data& data, Child& child)
    {
      // skip subdomain borders that do not coincide with the overall domain border
      if (!(child.appliesTo(data.elementSubDomains()) && child.appliesTo(data.neighborSubDomains())))
        return;
      typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSU lfsu(data.lfsu());
      LFSV lfsv(data.lfsv());
      LFSU lfsun(data.lfsun());
      LFSV lfsvn(data.lfsvn());
      LocalSparsityPattern localpattern_sn, localpattern_ns;
      child.localOperator().pattern_skeleton(lfsu,lfsv,lfsun,lfsvn,localpattern_sn,localpattern_ns);

      // translate local to global indices and add to global pattern
      for (size_t k=0; k<localpattern_sn.size(); ++k)
        add_entry(globalpattern,
                  lfsv.globalIndex(localpattern_sn[k].i()),
                  lfsun.globalIndex(localpattern_sn[k].j())
                  );

      for (size_t k=0; k<localpattern_ns.size(); ++k)
        add_entry(globalpattern,
                  lfsvn.globalIndex(localpattern_ns[k].i()),
                  lfsu.globalIndex(localpattern_ns[k].j())
                  );
    }

    P& globalpattern;
  };

  template<typename XL, typename RL>
  struct InvokeAlphaVolume
  {

    InvokeAlphaVolume(const XL& xl, RL& rl) :
      x(xl),
      r(rl)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSU lfsu(data.lfsu());
      LFSV lfsv(data.lfsv());
      child.localOperator().alpha_volume(ElementGeometry<typename Data::Entity>(data.element()),lfsu,x,lfsv,r);
    }

    const XL& x;
    RL& r;
  };

  template<typename RL>
  struct InvokeLambdaVolume
  {

    InvokeLambdaVolume(RL& rl) :
      r(rl)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSV lfsv(data.lfsv());
      child.localOperator().lambda_volume(ElementGeometry<typename Data::Entity>(data.element()),lfsv,r);
    }

    RL& r;
  };


  template<typename XL, typename RL>
  struct InvokeAlphaSkeletonOrBoundary
  {

    InvokeAlphaSkeletonOrBoundary(const XL& xl, const XL& xn, RL& rl, RL& rn, bool applyOneSided) :
      _xl(xl),
      _xn(xn),
      _rl(rl),
      _rn(rn),
      _applyOneSided(applyOneSided)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::LocalOperator LOP;
      typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSU lfsu(data.lfsu());
      LFSV lfsv(data.lfsv());
      if (child.appliesTo(data.neighborSubDomains()))
        {
          if (applyOneSided || LOP::doSkeletonTwoSided)
            {
              LFSU lfsun(data.lfsun());
              LFSV lfsvn(data.lfsvn());
              LocalAssemblerCallSwitch<Child,LOP::doAlphaSkeleton>::
                alpha_skeleton(child.localOperator(),
                               IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                 data.intersection_index),
                               lfsu,_xl,lfsv,lfsun,_xn,lfsvn,_rl,_rn);
              if(LOP::doAlphaSkeleton)
                data.setAlphaSkeletonInvoked();
            }
        }
      else
        {
          LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
            alpha_boundary(child.localOperator(),
                           IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                             data.intersection_index),
                           lfsu,_xl,lfsv,_rl);
          LocalAssemblerCallSwitch<LOP,LOP::doLambdaBoundary>::
            lambda_boundary(child.localOperator(),
                            IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                              data.intersection_index),
                            lfsv,_rl);
        }
    }

    const XL& _xl;
    const XL& _xn;
    RL& _rl;
    RL& _rn;
    const bool _applyOneSided;
  };


  template<typename XL, typename RL>
  struct InvokeAlphaBoundary
  {

    InvokeAlphaBoundary(const XL& xl, RL& rl) :
      x(xl),
      r(rl)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSU lfsu(data.lfsu());
      LFSV lfsv(data.lfsv());
      child.localOperator().alpha_boundary(IntersectionGeometry<typename Data::Intersection>(data.intersection(),data.intersectionIndex()),lfsu,x,lfsv,r);
    }

    const XL& x;
    RL& r;
  };

  template<typename RL>
  struct InvokeLambdaBoundary
  {

    InvokeLambdaBoundary(RL& rl) :
      r(rl)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSV lfsv(data.lfsv());
      child.localOperator().lambda_volume(IntersctionGeometry<typename Data::Intersection>(data.intersection(),data.intersectionIndex()),lfsv,r);
    }

    RL& r;
  };


  template<typename XL, typename RL>
  struct InvokeAlphaVolumePostSkeleton
  {

    InvokeAlphaVolumePostSkeleton(const XL& xl, RL& rl) :
      x(xl),
      r(rl)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSU lfsu(data.lfsu());
      LFSV lfsv(data.lfsv());
      child.localOperator().alpha_volume_post_skeleton(ElementGeometry<typename Data::Entity>(data.element()),lfsu,x,lfsv,r);
    }

    const XL& x;
    RL& r;
  };

  template<typename RL>
  struct InvokeLambdaVolumePostSkeleton
  {

    InvokeLambdaVolumePostSkeleton(RL& rl) :
      r(rl)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSV lfsv(data.lfsv());
      child.localOperator().lambda_volume_post_skeleton(ElementGeometry<typename Data::Entity>(data.element()),lfsv,r);
    }

    RL& r;
  };

public:
  typedef GridOperatorSpaceTraits<GFSU,GFSV,B,CU,CV> Traits;

  template<typename E>
  struct MatrixContainer
  {
    //! \brief define Type as the Type of a Matrix of E's
    typedef typename B::template Matrix<GridOperatorSpace,E> Type;
  private:
    MatrixContainer () {}
  };

  //! construct GridOperatorSpace
  GridOperatorSpace (const GFSU& gfsu_, const GFSV& gfsv_, const LA& la_)
    : gfsu(gfsu_), gfsv(gfsv_), la(la_)
  {
    pconstraintsu = &emptyconstraintsu;
    pconstraintsv = &emptyconstraintsv;
  }

  //! construct GridOperatorSpace, with constraints
  GridOperatorSpace (const GFSU& gfsu_, const CU& cu,
                     const GFSV& gfsv_, const CV& cv,
                     const LA& la_)
    : gfsu(gfsu_), gfsv(gfsv_), la(la_)
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

    operator_applier<MultiDomainGridOperatorSpace> apply_operator(*this);
    apply_operator.set(lfsu);
    apply_operator.set(lfsv);

    for (ElementIterator it = gfsu.gridview().template begin<0>();
         it!=gfsu.gridview().template end<0>(); ++it)
      {
        // bind local function spaces to element
        lfsu.bind(*it);
        lfsv.bind(*it);
        apply_operator.setElementSubDomains(gfsu.gridview().indexSet().subDomains(*it));

        apply_operator.conditional<do_pattern_volume>(BuildVolumePattern<P>(globalpattern));

        // skeleton and boundary pattern
        if (!any_child<MultiDomainGridOperatorSpace,do_pattern_skeleton>::value) continue;

        // local function spaces in neighbor
        LFSU lfsun(gfsu);
        LFSV lfsvn(gfsv);
        apply_operator.setn(lfsun);
        apply_operator.setn(lfsvn);

        IntersectionIterator endit = gfsu.gridview().iend(*it);
        for (IntersectionIterator iit = gfsu.gridview().ibegin(*it); iit!=endit; ++iit)
          {
            // skip if there is no neighbor
            if (!iit->neighbor()) continue;

            // bind local function spaces to neighbor element
            lfsun.bind(*(iit->outside()));
            lfsvn.bind(*(iit->outside()));
            apply_operator.setNeighborSubDomains(gfsu.gridview().indexSet().subDomains(*it));

            // get pattern
            apply_operator.conditional<do_pattern_skeleton>(BuildSkeletonPattern<P>(globalpattern));
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

    operator_applier<MultiDomainGridOperatorSpace> apply_operator(*this);
    apply_operator.set(lfsu);
    apply_operator.set(lfsv);

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
        apply_operator.setElementSubDomains(is.subDomains(*it));

        // allocate local data container
        typedef std::vector<typename X::ElementType> XL;
        XL xl(lfsu.size());
        typedef std::vector<typename R::ElementType> RL;
        RL rl(lfsv.size(),0.0);

        // read coefficents
        lfsu.vread(x,xl);

        // volume evaluation
        apply_operator.conditional<do_alpha_volume>(InvokeAlphaVolume<XL,RL>(xl,rl));
        apply_operator.conditional<do_lambda_volume>(InvokeLambdaVolume<RL>(rl));

        // skip if no intersection iterator is needed
        if (any_child<MultiDomainGridOperatorSpace,do_alpha_skeleton>::value ||
            any_child<MultiDomainGridOperatorSpace,do_alpha_boundary>::value ||
            any_child<MultiDomainGridOperatorSpace,do_lambda_skeleton>::value ||
            any_child<MultiDomainGridOperatorSpace,do_lambda_boundary>::value)
          {
            // local function spaces in neighbor
            LFSU lfsun(gfsu);
            LFSV lfsvn(gfsv);
            apply_operator.setn(lfsun);
            apply_operator.setn(lfsvn);

            // traverse intersections
            unsigned int intersection_index = 0;
            IntersectionIterator endit = gfsu.gridview().iend(*it);
            for (IntersectionIterator iit = gfsu.gridview().ibegin(*it);
                 iit!=endit; ++iit, ++intersection_index)
              {
                apply_operator.setIntersectionIndex(intersection_index);
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

                    // allocate local data container
                    XL xn(lfsun.size());
                    RL rn(lfsvn.size(),0.0);

                    // read coefficents
                    lfsun.vread(x,xn);

                    // unique vist of intersection
                    apply_operator.conditional<do_alpha_skeleton_or_boundary>
                      (InvokeAlphaSkeletonOrBoundary<XL,RL>(xl,xn,rl,rn,
                                                            id > idn ||
                                                            (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                            )
                       );
                    if (data.alphaSkeletonInvoked())
                      {
                        lfsvn.vadd(rn,r);
                        data.clearAlphaSkeletonInvoked();
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    apply_operator.conditional<do_alpha_boundary>(InvokeAlphaBoundary<XL,RL>(xl,rl));
                    apply_operator.conditional<do_lambda_boundary>(InvokeLambdaBoundary<RL>(rl));
                  }
              }
          }

        apply_operator.conditional<do_alpha_volume_post_skeleton>(InvokeAlphaVolumePostSkeleton<XL,RL>(xl,rl));
        apply_operator.conditional<do_lambda_volume_post_skeleton>(InvokeLambdaVolumePostSkeleton<RL>(rl));

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

        // allocate local data container
        std::vector<typename X::ElementType> xl(lfsu.size());
        std::vector<typename Y::ElementType> yl(lfsv.size(),0.0);

        // read coefficents
        lfsu.vread(x,xl);

        // volume evaluation
        LocalAssemblerCallSwitch<LA,LA::doAlphaVolume>::
          jacobian_apply_volume(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,yl);

        // skeleton and boundary evaluation
        if (LA::doAlphaSkeleton||LA::doAlphaBoundary)
          {
            // local function spaces in neighbor
            LFSU lfsun(gfsu);
            LFSV lfsvn(gfsv);

            unsigned int intersection_index = 0;
            IntersectionIterator endit = gfsu.gridview().iend(*it);
            for (IntersectionIterator iit = gfsu.gridview().ibegin(*it);
                 iit!=endit; ++iit, ++intersection_index)
              {
                // skeleton term
                if (iit->neighbor() && LA::doAlphaSkeleton )
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

                    // unique vist of intersection
                    if (LA::doSkeletonTwoSided || id>idn ||
                        (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity))
                      {
                        // bind local function spaces to neighbor element
                        lfsun.bind(*(iit->outside()));
                        lfsvn.bind(*(iit->outside()));

                        // allocate local data container
                        std::vector<typename X::ElementType> xn(lfsun.size());
                        std::vector<typename Y::ElementType> yn(lfsvn.size(),0.0);

                        // read coefficents
                        lfsun.vread(x,xn);

                        // skeleton evaluation
                        LocalAssemblerCallSwitch<LA,LA::doAlphaSkeleton>::
                          jacobian_apply_skeleton(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsu,xl,lfsv,lfsun,xn,lfsvn,yl,yn);

                        // accumulate result (note: r needs to be cleared outside)
                        lfsvn.vadd(yn,y);
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    LocalAssemblerCallSwitch<LA,LA::doAlphaBoundary>::
                      jacobian_apply_boundary(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsu,xl,lfsv,yl);
                  }
              }
          }

        LocalAssemblerCallSwitch<LA,LA::doAlphaVolumePostSkeleton>::
          jacobian_apply_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,yl);

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
        //            std::cout << "[" << gfsu.gridview().comm().rank() << "] " << " element: " << id << std::endl;

        // skip ghost and overlap
        if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
          continue;

        // bind local function spaces to element
        lfsu.bind(*it);
        lfsv.bind(*it);

        // allocate local data container
        std::vector<typename X::ElementType> xl(lfsu.size());
        LocalMatrix<typename A::ElementType> al(lfsv.size(),lfsu.size(),0.0);

        // read coefficents
        lfsu.vread(x,xl);

        // volume evaluation
        LocalAssemblerCallSwitch<LA,LA::doAlphaVolume>::
          jacobian_volume(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,al);

        // skeleton and boundary evaluation
        if (LA::doAlphaSkeleton||LA::doAlphaBoundary)
          {
            // local function spaces in neighbor
            LFSU lfsun(gfsu);
            LFSV lfsvn(gfsv);

            unsigned int intersection_index = 0;
            IntersectionIterator endit = gfsu.gridview().iend(*it);
            for (IntersectionIterator iit = gfsu.gridview().ibegin(*it);
                 iit!=endit; ++iit, ++intersection_index)
              {
                // skeleton term
                if (iit->neighbor() && LA::doAlphaSkeleton )
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

                    // unique vist of intersection
                    if (LA::doSkeletonTwoSided || id>idn ||
                        (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity) )
                      {
                        // bind local function spaces to neighbor element
                        lfsun.bind(*(iit->outside()));
                        lfsvn.bind(*(iit->outside()));

                        // allocate local data container
                        std::vector<typename X::ElementType> xn(lfsun.size());
                        LocalMatrix<typename A::ElementType> al_sn(lfsv.size() ,lfsun.size(),0.0);
                        LocalMatrix<typename A::ElementType> al_ns(lfsvn.size(),lfsu.size() ,0.0);
                        LocalMatrix<typename A::ElementType> al_nn(lfsvn.size(),lfsun.size(),0.0);

                        // read coefficents
                        lfsun.vread(x,xn);

                        // skeleton evaluation
                        LocalAssemblerCallSwitch<LA,LA::doAlphaSkeleton>::
                          jacobian_skeleton(la,IntersectionGeometry<Intersection>(*iit,intersection_index),
                                            lfsu,xl,lfsv,lfsun,xn,lfsvn,al,al_sn,al_ns,al_nn);

                        // accumulate result
                        etadd(lfsv,lfsun,al_sn,a);
                        etadd(lfsvn,lfsu,al_ns,a);
                        etadd(lfsvn,lfsun,al_nn,a);
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    LocalAssemblerCallSwitch<LA,LA::doAlphaBoundary>::
                      jacobian_boundary(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsu,xl,lfsv,al);
                  }
              }
          }

        LocalAssemblerCallSwitch<LA,LA::doAlphaVolumePostSkeleton>::
          jacobian_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,al);

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
  const LA& la;
  const CU* pconstraintsu;
  const CV* pconstraintsv;
  CU emptyconstraintsu;
  CV emptyconstraintsv;
};

//! \} group GridFunctionSpace

#endif

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_MULTIDOMAINGRIDOPERATORSPACE_HH