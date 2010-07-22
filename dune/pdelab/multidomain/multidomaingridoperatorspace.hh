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
template<typename GFSU, typename GFSV,
         typename B,
         typename... SubProblemsAndCouplings>
class MultiDomainGridOperatorSpace : public VariadicCompositeNode<CopyStoragePolicy,SubProblemsAndCouplings...>
{
  typedef VariadicCompositeNode<CopyStoragePolicy,SubProblemsAndCouplings...> BaseT;

  typedef typename extract_problems<SubProblemsAndCouplings...>::type SubProblemList;
  typedef typename extract_problems<SubProblemsAndCouplings...>::map_type SubProblemMap;
  typedef typename extract_couplings<SubProblemsAndCouplings...>::type CouplingList;

  static const std::size_t subProblemCount = std::tuple_size<SubProblemList>::value;
  static const std::size_t couplingCount = std::tuple_size<SubProblemList>::value;

  typedef typename GFSU::template ConstraintsContainer<double>::Type CU;
  typedef typename GFSV::template ConstraintsContainer<double>::Type CV;

  template<std::size_t k>
  struct SubProblem {
    typedef typename std::tuple_element<k,SubProblemList>::type::type Type;
  };

  template<std::size_t k>
  const typename SubProblem<k>::Type& subProblem() const {
    return this->template getChild<std::tuple_element<k,SubProblemList>::type::globalIndex>();
  }

  template<std::size_t k>
  struct Coupling {
    typedef typename std::tuple_element<k,CouplingList>::type::type Type;
  };

  template<std::size_t k>
  const typename Coupling<k>::Type& coupling() const {
    return this->template getChild<std::tuple_element<k,CouplingList>::type::globalIndex>();
  }

  // extract useful types
  typedef typename GFSU::Traits::GridType Grid;
  typedef typename Grid::LeafGridView GV;
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
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      LocalSparsityPattern localpattern;
      child.localOperator().pattern_volume(lfsu,lfsv,localpattern);

      // translate local to global indices and add to global pattern
      for (size_t k=0; k<localpattern.size(); ++k)
        // TODO: Figure out if we really need the localIndex() call
        data.gos().add_entry(globalpattern,
                             lfsv.globalIndex(lfsv.localIndex(localpattern[k].i())),
                             lfsu.globalIndex(lfsu.localIndex(localpattern[k].j()))
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
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      LFSU lfsun(data.lfsun(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsvn(data.lfsvn(),child,child.testGridFunctionSpaceConstraints());
      lfsun.bind();
      lfsvn.bind();
      LocalSparsityPattern localpattern_sn, localpattern_ns;
      child.localOperator().pattern_skeleton(lfsu,lfsv,lfsun,lfsvn,localpattern_sn,localpattern_ns);

      // translate local to global indices and add to global pattern
      for (size_t k=0; k<localpattern_sn.size(); ++k)
        data.gos().add_entry(globalpattern,
                             data.lfsv().globalIndex(data.lfsv().localIndex(localpattern_sn[k].i())),
                             data.lfsun().globalIndex(data.lfsun().localIndex(localpattern_sn[k].j()))
                             );

      for (size_t k=0; k<localpattern_ns.size(); ++k)
        data.gos().add_entry(globalpattern,
                             data.lfsvn().globalIndex(data.lfsvn().localIndex(localpattern_ns[k].i())),
                             data.lfsu().globalIndex(data.lfsu().localIndex(localpattern_ns[k].j()))
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
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      child.localOperator().alpha_volume(ElementGeometry<typename Data::Element>(data.element()),lfsu,x,lfsv,r);
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
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsv.bind();
      child.localOperator().lambda_volume(ElementGeometry<typename Data::Element>(data.element()),lfsv,r);
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
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      if (child.appliesTo(data.neighborSubDomains()))
        {
          if (_applyOneSided || LOP::doSkeletonTwoSided)
            {
              LFSU lfsun(data.lfsun(),child,child.trialGridFunctionSpaceConstraints());
              LFSV lfsvn(data.lfsvn(),child,child.testGridFunctionSpaceConstraints());
              lfsun.bind();
              lfsvn.bind();
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
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      child.localOperator().alpha_boundary(IntersectionGeometry<typename Data::Intersection>(data.intersection(),data.intersectionIndex()),lfsu,x,lfsv,r);
    }

    const XL& x;
    RL& r;
  };

  template<typename Operator, typename XL, typename RL>
  struct InvokeAlphaCoupling
  {

    InvokeAlphaCoupling(const XL& xl, const XL& xn, RL& rl, RL& rn) :
      _local_x(xl),
      _remote_x(xn),
      _local_r(rl),
      _remote_r(rn)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains(),data.neighborSubDomains()))
        return;
      typedef typename Operator::template ExtractType<Child>::Type LOP;
      typedef typename Child::Traits::LocalSubProblem LocalSubProblem;
      typedef typename Child::Traits::RemoteSubProblem RemoteSubProblem;
      typedef typename LocalSubProblem::TrialLocalFunctionSpace LocalLFSU;
      typedef typename LocalSubProblem::TestLocalFunctionSpace LocalLFSV;
      typedef typename RemoteSubProblem::TrialLocalFunctionSpace RemoteLFSU;
      typedef typename RemoteSubProblem::TestLocalFunctionSpace RemoteLFSV;
      const LocalSubProblem& localSubProblem = child.localSubProblem();
      const RemoteSubProblem& remoteSubProblem = child.remoteSubProblem();
      LocalLFSU local_lfsu(data.lfsu(),localSubProblem,localSubProblem.trialGridFunctionSpaceConstraints());
      LocalLFSV local_lfsv(data.lfsv(),localSubProblem,localSubProblem.testGridFunctionSpaceConstraints());
      local_lfsu.bind();
      local_lfsv.bind();
      RemoteLFSU remote_lfsu(data.lfsun(),remoteSubProblem,remoteSubProblem.trialGridFunctionSpaceConstraints());
      RemoteLFSV remote_lfsv(data.lfsvn(),remoteSubProblem,remoteSubProblem.testGridFunctionSpaceConstraints());
      remote_lfsu.bind();
      remote_lfsv.bind();
      Operator::extract(child).alpha_coupling(IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                                data.intersection_index),
                                              local_lfsu,
                                              _local_x,
                                              local_lfsv,
                                              remote_lfsu,
                                              _remote_x,
                                              remote_lfsv,
                                              _local_r,
                                              _remote_r);
      if(LOP::doAlphaSkeleton)
        data.setAlphaSkeletonInvoked();
    }

    const XL& _local_x;
    const XL& _remote_x;
    RL& _local_r;
    RL& _remote_r;
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
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsv.bind();
      child.localOperator().lambda_boundary(IntersectionGeometry<typename Data::Intersection>(data.intersection(),data.intersectionIndex()),lfsv,r);
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
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      child.localOperator().alpha_volume_post_skeleton(ElementGeometry<typename Data::Element>(data.element()),lfsu,x,lfsv,r);
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
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsv.bind();
      child.localOperator().lambda_volume_post_skeleton(ElementGeometry<typename Data::Element>(data.element()),lfsv,r);
    }

    RL& r;
  };

  template<typename XL, typename YL>
  struct InvokeJacobianApplyVolume
  {

    InvokeJacobianApplyVolume(const XL& xl, YL& yl) :
      x(xl),
      y(yl)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      child.localOperator().jacobian_apply_volume(ElementGeometry<typename Data::Element>(data.element()),lfsu,x,lfsv,y);
    }

    const XL& x;
    YL& y;
  };

  template<typename XL, typename YL>
  struct InvokeJacobianApplySkeletonOrBoundary
  {

    InvokeJacobianApplySkeletonOrBoundary(const XL& xl, const XL& xn, YL& yl, YL& yn, bool applyOneSided) :
      _xl(xl),
      _xn(xn),
      _yl(yl),
      _yn(yn),
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
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      if (child.appliesTo(data.neighborSubDomains()))
        {
          if (_applyOneSided || LOP::doSkeletonTwoSided)
            {
              LFSU lfsun(data.lfsun(),child,child.trialGridFunctionSpaceConstraints());
              LFSV lfsvn(data.lfsvn(),child,child.testGridFunctionSpaceConstraints());
              lfsun.bind();
              lfsvn.bind();
              LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
                jacobian_apply_skeleton(child.localOperator(),
                                        IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                          data.intersection_index),
                                        lfsu,_xl,lfsv,lfsun,_xn,lfsvn,_yl,_yn);
              if(LOP::doAlphaSkeleton)
                data.setAlphaSkeletonInvoked();
            }
        }
      else
        {
          LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
            jacobian_apply_boundary(child.localOperator(),
                                    IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                      data.intersection_index),
                                    lfsu,_xl,lfsv,_yl);
        }
    }

    const XL& _xl;
    const XL& _xn;
    YL& _yl;
    YL& _yn;
    const bool _applyOneSided;
  };

  template<typename XL, typename YL>
  struct InvokeJacobianApplyBoundary
  {

    InvokeJacobianApplyBoundary(const XL& xl, YL& yl) :
      x(xl),
      y(yl)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      child.localOperator().jacobian_apply_boundary(IntersectionGeometry<typename Data::Intersection>(data.intersection(),data.intersectionIndex()),lfsu,x,lfsv,y);
    }

    const XL& x;
    YL& y;
  };


  template<typename Operator, typename XL, typename YL>
  struct InvokeJacobianApplyCoupling
  {

    InvokeJacobianApplyCoupling(const XL& xl, const XL& xn, YL& yl, YL& yn) :
      _local_x(xl),
      _remote_x(xn),
      _local_y(yl),
      _remote_y(yn)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains(),data.neighborSubDomains()))
        return;
      typedef typename Operator::template ExtractType<Child>::Type LOP;
      typedef typename Child::Traits::LocalSubProblem LocalSubProblem;
      typedef typename Child::Traits::RemoteSubProblem RemoteSubProblem;
      typedef typename LocalSubProblem::TrialLocalFunctionSpace LocalLFSU;
      typedef typename LocalSubProblem::TestLocalFunctionSpace LocalLFSV;
      typedef typename RemoteSubProblem::TrialLocalFunctionSpace RemoteLFSU;
      typedef typename RemoteSubProblem::TestLocalFunctionSpace RemoteLFSV;
      const LocalSubProblem& localSubProblem = child.localSubProblem();
      const RemoteSubProblem& remoteSubProblem = child.remoteSubProblem();
      LocalLFSU local_lfsu(data.lfsu(),localSubProblem,localSubProblem.trialGridFunctionSpaceConstraints());
      LocalLFSV local_lfsv(data.lfsv(),localSubProblem,localSubProblem.testGridFunctionSpaceConstraints());
      local_lfsu.bind();
      local_lfsv.bind();
      RemoteLFSU remote_lfsu(data.lfsun(),remoteSubProblem,remoteSubProblem.trialGridFunctionSpaceConstraints());
      RemoteLFSV remote_lfsv(data.lfsvn(),remoteSubProblem,remoteSubProblem.testGridFunctionSpaceConstraints());
      remote_lfsu.bind();
      remote_lfsv.bind();
      Operator::extract(child).jacobian_apply_coupling(IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                                         data.intersection_index),
                                                       local_lfsu,
                                                       _local_x,
                                                       local_lfsv,
                                                       remote_lfsu,
                                                       _remote_x,
                                                       remote_lfsv,
                                                       _local_y,
                                                       _remote_y);
      if(LOP::doAlphaSkeleton)
        data.setAlphaSkeletonInvoked();
    }

    const XL& _local_x;
    const XL& _remote_x;
    YL& _local_y;
    YL& _remote_y;
  };


  template<typename XL, typename YL>
  struct InvokeJacobianApplyVolumePostSkeleton
  {

    InvokeJacobianApplyVolumePostSkeleton(const XL& xl, YL& yl) :
      x(xl),
      y(yl)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      child.localOperator().jacobian_apply_volume_post_skeleton(ElementGeometry<typename Data::Element>(data.element()),lfsu,x,lfsv,y);
    }

    const XL& x;
    YL& y;
  };


  template<typename XL, typename AL>
  struct InvokeJacobianVolume
  {

    InvokeJacobianVolume(const XL& xl, AL& al) :
      x(xl),
      a(al)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::LocalTrialFunctionSpace LFSU;
      typedef typename Child::Traits::LocalTestFunctionSpace LFSV;
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      child.localOperator().jacobian_volume(ElementGeometry<typename Data::Element>(data.element()),lfsu,x,lfsv,a);
    }

    const XL& x;
    AL& a;
  };

  template<typename XL, typename AL>
  struct InvokeJacobianSkeletonOrBoundary
  {

    InvokeJacobianSkeletonOrBoundary(const XL& xl, const XL& xn, AL& al, AL& al_sn, AL& al_ns, AL& al_nn, bool applyOneSided) :
      _xl(xl),
      _xn(xn),
      _al(al),
      _al_sn(al_sn),
      _al_ns(al_ns),
      _al_nn(al_nn),
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
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      if (child.appliesTo(data.neighborSubDomains()))
        {
          if (_applyOneSided || LOP::doSkeletonTwoSided)
            {
              LFSU lfsun(data.lfsun(),child,child.trialGridFunctionSpaceConstraints());
              LFSV lfsvn(data.lfsvn(),child,child.testGridFunctionSpaceConstraints());
              lfsun.bind();
              lfsvn.bind();
              LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
                jacobian_skeleton(child.localOperator(),
                                  IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                    data.intersectionIndex()),
                                  lfsu,_xl,lfsv,lfsun,_xn,lfsvn,_al,_al_sn,_al_ns,_al_nn);
              if(LOP::doAlphaSkeleton)
                data.setAlphaSkeletonInvoked();
            }
        }
      else
        {
          LocalAssemblerCallSwitch<LOP,LOP::doAlphaBoundary>::
            jacobian_boundary(child.localOperator(),
                              IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                data.intersectionIndex()),
                              lfsu,_xl,lfsv,_al);
        }
    }

    const XL& _xl;
    const XL& _xn;
    AL& _al;
    AL& _al_sn;
    AL& _al_ns;
    AL& _al_nn;
    const bool _applyOneSided;
  };

  template<typename XL, typename AL>
  struct InvokeJacobianBoundary
  {

    InvokeJacobianBoundary(const XL& xl, AL& al) :
      x(xl),
      a(al)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      child.localOperator().jacobian_boundary(IntersectionGeometry<typename Data::Intersection>(data.intersection(),data.intersectionIndex()),lfsu,x,lfsv,a);
    }

    const XL& x;
    AL& a;
  };

  template<typename Operator, typename XL, typename AL>
  struct InvokeJacobianCoupling
  {

    InvokeJacobianCoupling(const XL& xl, const XL& xn, AL& al, AL& al_sn, AL& al_ns, AL& al_nn) :
      _local_x(xl),
      _remote_x(xn),
      _local_a(al),
      _local_to_remote_a(al_sn),
      _remote_to_local_a(al_ns),
      _remote_a(al_nn)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains(),data.neighborSubDomains()))
        return;
      typedef typename Operator::template ExtractType<Child>::Type LOP;
      typedef typename Child::Traits::LocalSubProblem LocalSubProblem;
      typedef typename Child::Traits::RemoteSubProblem RemoteSubProblem;
      typedef typename LocalSubProblem::TrialLocalFunctionSpace LocalLFSU;
      typedef typename LocalSubProblem::TestLocalFunctionSpace LocalLFSV;
      typedef typename RemoteSubProblem::TrialLocalFunctionSpace RemoteLFSU;
      typedef typename RemoteSubProblem::TestLocalFunctionSpace RemoteLFSV;
      const LocalSubProblem& localSubProblem = child.localSubProblem();
      const RemoteSubProblem& remoteSubProblem = child.remoteSubProblem();
      LocalLFSU local_lfsu(data.lfsu(),localSubProblem,localSubProblem.trialGridFunctionSpaceConstraints());
      LocalLFSV local_lfsv(data.lfsv(),localSubProblem,localSubProblem.testGridFunctionSpaceConstraints());
      local_lfsu.bind();
      local_lfsv.bind();
      RemoteLFSU remote_lfsu(data.lfsun(),remoteSubProblem,remoteSubProblem.trialGridFunctionSpaceConstraints());
      RemoteLFSV remote_lfsv(data.lfsvn(),remoteSubProblem,remoteSubProblem.testGridFunctionSpaceConstraints());
      remote_lfsu.bind();
      remote_lfsv.bind();
      Operator::extract(child).jacobian_apply_coupling(IntersectionGeometry<typename Data::Intersection>(data.intersection(),
                                                                                                         data.intersection_index),
                                                       local_lfsu,
                                                       _local_x,
                                                       local_lfsv,
                                                       remote_lfsu,
                                                       _remote_x,
                                                       remote_lfsv,
                                                       _local_a,
                                                       _local_to_remote_a,
                                                       _remote_to_local_a,
                                                       _remote_a);
      if(LOP::doAlphaSkeleton)
        data.setAlphaSkeletonInvoked();
    }

    const XL& _local_x;
    const XL& _remote_x;
    AL& _local_a;
    AL& _local_to_remote_a;
    AL& _remote_to_local_a;
    AL& _remote_a;
  };


  template<typename XL, typename AL>
  struct InvokeJacobianVolumePostSkeleton
  {

    InvokeJacobianVolumePostSkeleton(const XL& xl, AL& al) :
      x(xl),
      a(al)
    {}

    template<typename Data, typename Child>
    void operator()(Data& data, const Child& child)
    {
      if (!child.appliesTo(data.elementSubDomains()))
        return;
      typedef typename Child::Traits::TrialLocalFunctionSpace LFSU;
      typedef typename Child::Traits::TestLocalFunctionSpace LFSV;
      LFSU lfsu(data.lfsu(),child,child.trialGridFunctionSpaceConstraints());
      LFSV lfsv(data.lfsv(),child,child.testGridFunctionSpaceConstraints());
      lfsu.bind();
      lfsv.bind();
      child.localOperator().jacobian_volume_post_skeleton(ElementGeometry<typename Data::Element>(data.element()),lfsu,x,lfsv,a);
    }

    const XL& x;
    AL& a;
  };


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
    data::IntersectionReference,
    data::SkeletonInvocationTracker
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

        apply_operator.template conditional<do_pattern_volume<> >(BuildVolumePattern<P>(globalpattern));

        // skeleton and boundary pattern
        if (!any_child<MultiDomainGridOperatorSpace,do_pattern_skeleton<> >::value) continue;

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
            apply_operator.template conditional<do_pattern_skeleton<> >(BuildSkeletonPattern<P>(globalpattern));
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
        apply_operator.template conditional<do_alpha_volume<> >(InvokeAlphaVolume<XL,RL>(xl,rl));
        apply_operator.template conditional<do_lambda_volume<> >(InvokeLambdaVolume<RL>(rl));

        // skip if no intersection iterator is needed
        if (any_child<MultiDomainGridOperatorSpace,do_alpha_skeleton<> >::value ||
            any_child<MultiDomainGridOperatorSpace,do_alpha_boundary<> >::value ||
            any_child<MultiDomainGridOperatorSpace,do_lambda_boundary<> >::value ||
            any_child<MultiDomainGridOperatorSpace,do_alpha_coupling<> >::value)
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
                    apply_operator.template conditional<do_alpha_skeleton_or_boundary<> >
                      (InvokeAlphaSkeletonOrBoundary<XL,RL>(xl,xn,rl,rn,
                                                            id > idn ||
                                                            (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                            )
                       );
                    apply_operator.template conditional<do_alpha_coupling<> >
                      (InvokeAlphaCoupling<CouplingOperator,XL,RL>(xl,xn,rl,rn));
                    if (apply_operator.alphaSkeletonInvoked())
                      {
                        lfsvn.vadd(rn,r);
                        apply_operator.clearAlphaSkeletonInvoked();
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    apply_operator.template conditional<do_alpha_boundary<> >(InvokeAlphaBoundary<XL,RL>(xl,rl));
                    apply_operator.template conditional<do_lambda_boundary<> >(InvokeLambdaBoundary<RL>(rl));
                  }
              }
          }

        apply_operator.template conditional<do_alpha_volume_post_skeleton<> >(InvokeAlphaVolumePostSkeleton<XL,RL>(xl,rl));
        apply_operator.template conditional<do_lambda_volume_post_skeleton<> >(InvokeLambdaVolumePostSkeleton<RL>(rl));

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
        apply_operator.template conditional<do_alpha_volume<> >(InvokeJacobianApplyVolume<XL,YL>(xl,yl));

        // skeleton and boundary evaluation
        if (any_child<MultiDomainGridOperatorSpace,do_alpha_skeleton<> >::value ||
            any_child<MultiDomainGridOperatorSpace,do_alpha_boundary<> >::value ||
            any_child<MultiDomainGridOperatorSpace,do_alpha_coupling<> >::value)
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
                    YL yn(lfsvn.size(),0.0);

                    // read coefficents
                    lfsun.vread(x,xn);

                    apply_operator.template conditional<do_alpha_skeleton_or_boundary<> >
                      (InvokeJacobianApplySkeletonOrBoundary<XL,YL>(xl,xn,yl,yn,
                                                                    id > idn ||
                                                                    (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                                    )
                       );
                    apply_operator.template conditional<do_alpha_coupling<> >
                      (InvokeJacobianApplyCoupling<CouplingOperator,XL,YL>(xl,xn,yl,yn));
                    if (apply_operator.alphaSkeletonInvoked())
                      {
                        lfsvn.vadd(yn,y);
                        apply_operator.clearAlphaSkeletonInvoked();
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    apply_operator.template conditional<do_alpha_boundary<> >(InvokeJacobianApplyBoundary<XL,YL>(xl,yl));
                  }
              }
          }

        apply_operator.template conditional<do_alpha_volume_post_skeleton<> >(InvokeJacobianApplyVolumePostSkeleton<XL,YL>(xl,yl));

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
        apply_operator.template conditional<do_alpha_volume<> >(InvokeJacobianVolume<XL,AL>(xl,al));

        // skeleton and boundary evaluation
        if (any_child<MultiDomainGridOperatorSpace,do_alpha_skeleton<> >::value ||
            any_child<MultiDomainGridOperatorSpace,do_alpha_boundary<> >::value ||
            any_child<MultiDomainGridOperatorSpace,do_alpha_coupling<> >::value)
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

                    apply_operator.template conditional<do_alpha_skeleton_or_boundary<> >
                      (InvokeJacobianSkeletonOrBoundary<XL,AL>(xl,xn,al,al_sn,al_ns,al_nn,
                                                                    id > idn ||
                                                                    (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity)
                                                                    )
                       );
                    apply_operator.template conditional<do_alpha_coupling<> >
                      (InvokeJacobianCoupling<CouplingOperator,XL,AL>(xl,al,al_sn,al_ns,al_nn));
                    if (apply_operator.alphaSkeletonInvoked())
                      {
                        etadd(lfsv,lfsun,al_sn,a);
                        etadd(lfsvn,lfsu,al_ns,a);
                        etadd(lfsvn,lfsun,al_nn,a);
                        apply_operator.clearAlphaSkeletonInvoked();
                      }
                  }

                // boundary term
                if (iit->boundary())
                  {
                    apply_operator.template conditional<do_alpha_boundary<> >(InvokeJacobianBoundary<XL,AL>(xl,al));
                  }
              }
          }

        apply_operator.template conditional<do_alpha_volume_post_skeleton<> >(InvokeJacobianVolumePostSkeleton<XL,AL>(xl,al));

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
