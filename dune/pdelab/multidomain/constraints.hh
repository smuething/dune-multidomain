#ifndef DUNE_PDELAB_MULTIDOMAIN_CONSTRAINTS_HH
#define DUNE_PDELAB_MULTIDOMAIN_CONSTRAINTS_HH

#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/multidomain/operatorflagtests.hh>
#include <dune/pdelab/multidomain/constraintsfunctors.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {


#ifndef DOXYGEN

template<typename MDLFS, typename Constraints>
struct BuildConstraintsDescriptor;

struct NoConstraintsParameters:
    public Dune::PDELab::NoConstraintsParameters
{

  template<typename IG, typename X>
  bool isDirichlet(const IG& ig, const X& x) const
  {
    return false;
  }

  template<typename IG, typename X>
  bool isNeumann(const IG& ig, const X& x) const
  {
    return false;
  }

};

template<typename Parameters>
class ParameterHolder
{

public:

  static const bool noParameters = is_same<Parameters,NoConstraintsParameters>::value;

  const Parameters& parameters() const
  {
    return _parameters;
  }

protected:

  ParameterHolder()
  {}

  ParameterHolder(const Parameters& parameters)
    : _parameters(parameters)
  {}

private:
  //! Save object if type is NoConstraintsParameters, const reference otherwise
  typename std::conditional<noParameters,Parameters,const Parameters&>::type _parameters;

};

#endif // DOXYGEN

template<typename Parameters>
class ConstrainMultiDomainGridFunctionSpace
  : public ParameterHolder<Parameters>
{

public:

  ConstrainMultiDomainGridFunctionSpace()
  {}

  ConstrainMultiDomainGridFunctionSpace(const Parameters& parameters)
    : ParameterHolder<Parameters>(parameters)
  {}

};

template<typename SubProblem, typename Parameters>
class ConstrainSubProblem
  : public ParameterHolder<Parameters>
{

public:

  ConstrainSubProblem(const SubProblem& subProblem)
    : _subProblem(subProblem)
  {}

  ConstrainSubProblem(const SubProblem& subProblem, const Parameters& parameters)
    : ParameterHolder<Parameters>(parameters)
    , _subProblem(subProblem)
  {}

  const SubProblem& subProblem() const
  {
    return _subProblem;
  }

private:

  const SubProblem& _subProblem;

};


ConstrainMultiDomainGridFunctionSpace<NoConstraintsParameters>
constrainMultiDomainGridFunctionSpace()
{
  return ConstrainMultiDomainGridFunctionSpace<NoConstraintsParameters>();
}

template<typename Parameters>
ConstrainMultiDomainGridFunctionSpace<Parameters>
constrainMultiDomainGridFunctionSpace(const Parameters& parameters)
{
  return ConstrainMultiDomainGridFunctionSpace<Parameters>(parameters);
}

template<typename SubProblem>
ConstrainSubProblem<SubProblem,NoConstraintsParameters>
constrainSubProblem(const SubProblem& subProblem)
{
  return ConstrainSubProblem<SubProblem,NoConstraintsParameters>(subProblem);
}

template<typename SubProblem, typename Parameters>
ConstrainSubProblem<SubProblem,Parameters>
constrainSubProblem(const SubProblem& subProblem, const Parameters& parameters)
{
  return ConstrainSubProblem<SubProblem,Parameters>(subProblem,parameters);
}



template<typename MDLFS, typename Parameters>
struct MultiDomainGridFunctionSpaceConstraints
  : public Dune::PDELab::TypeTree::LeafNode
  , public ParameterHolder<Parameters>
{

  typedef ParameterHolder<Parameters> BaseT;

  static const bool doVolume = any_child<MDLFS,do_constraints_volume<LFSConstraintsExtractor> >::value;

  static const bool doSkeleton = any_child<MDLFS,do_constraints_skeleton<LFSConstraintsExtractor> >::value;

  static const bool doBoundary = any_child<MDLFS,do_constraints_boundary<LFSConstraintsExtractor> >::value;

  static const bool doProcessor = any_child<MDLFS,do_constraints_processor<LFSConstraintsExtractor> >::value;

  template<typename EG>
  bool appliesTo(const EG& eg) const
  {
    return true;
  }

  void rebind_lfs_s() const
  {}

  void rebind_lfs_n() const
  {}

  const MDLFS& lfs_s() const
  {
    return _mdlfs_s;
  }

  const MDLFS& lfs_n() const
  {
    return _mdlfs_n;
  }

  MultiDomainGridFunctionSpaceConstraints(const MDLFS& mdlfs_s,
                                          const MDLFS& mdlfs_n)
    : _mdlfs_s(mdlfs_s)
    , _mdlfs_n(mdlfs_n)
  {}

  MultiDomainGridFunctionSpaceConstraints(const MDLFS& mdlfs_s,
                                          const MDLFS& mdlfs_n,
                                          const Parameters& parameters)
    : BaseT(parameters)
    , _mdlfs_s(mdlfs_s)
    , _mdlfs_n(mdlfs_n)
  {}

private:

  const MDLFS& _mdlfs_s;
  const MDLFS& _mdlfs_n;

};

template<typename MDLFS, typename Parameters>
struct BuildConstraintsDescriptor<MDLFS,ConstrainMultiDomainGridFunctionSpace<Parameters> >
{
  typedef MultiDomainGridFunctionSpaceConstraints<MDLFS,Parameters> type;
};

template<typename MDLFS, typename Parameters>
shared_ptr<MultiDomainGridFunctionSpaceConstraints<MDLFS,Parameters> >
buildConstraintsDescriptor(const MDLFS& mdlfs_s, const MDLFS& mdlfs_n,
                           const ConstrainMultiDomainGridFunctionSpace<Parameters>& d)
{
  return make_shared<MultiDomainGridFunctionSpaceConstraints<MDLFS,Parameters> >(mdlfs_s,mdlfs_n,d.parameters());
}


template<typename MDLFS, typename SubProblem, typename Parameters>
struct SubProblemConstraints
  : public Dune::PDELab::TypeTree::LeafNode
  , public ParameterHolder<Parameters>
{

  typedef ParameterHolder<Parameters> BaseT;

  typedef typename SubProblem::Traits::MultiDomainTrialGridFunctionSpace TrialGFS;

  typedef typename std::conditional<is_same<typename MDLFS::Traits::GridFunctionSpaceType,
                                            TrialGFS
                                            >::value,
                                    typename SubProblem::Traits::TrialLocalFunctionSpace,
                                    typename SubProblem::Traits::TestLocalFunctionSpace
                                    >::type SubProblemLocalFunctionSpace;

  static const bool doVolume = any_child<SubProblemLocalFunctionSpace,do_constraints_volume<LFSConstraintsExtractor> >::value;

  static const bool doSkeleton = any_child<SubProblemLocalFunctionSpace,do_constraints_skeleton<LFSConstraintsExtractor> >::value;

  static const bool doBoundary = any_child<SubProblemLocalFunctionSpace,do_constraints_boundary<LFSConstraintsExtractor> >::value;

  static const bool doProcessor = any_child<SubProblemLocalFunctionSpace,do_constraints_processor<LFSConstraintsExtractor> >::value;

  template<typename EG>
  bool appliesTo(const EG& eg) const
  {
    return _subProblem.appliesTo(eg);
  }

  void rebind_lfs_s()
  {
    _splfs_s.bind();
  }

  void rebind_lfs_n()
  {
    _splfs_n.bind();
  }

  const SubProblemLocalFunctionSpace& lfs_s() const
  {
    return _splfs_s;
  }

  const SubProblemLocalFunctionSpace& lfs_n() const
  {
    return _splfs_n;
  }

  SubProblemConstraints(const MDLFS& mdlfs_s,
                        const MDLFS& mdlfs_n,
                        const SubProblem& subProblem)
    : _subProblem(subProblem)
    , _splfs_s(mdlfs_s,subProblem)
    , _splfs_n(mdlfs_n,subProblem)
  {}

  SubProblemConstraints(const MDLFS& mdlfs_s,
                        const MDLFS& mdlfs_n,
                        const SubProblem& subProblem,
                        const Parameters& parameters)
    : BaseT(parameters)
    , _subProblem(subProblem)
    , _splfs_s(mdlfs_s,subProblem)
    , _splfs_n(mdlfs_n,subProblem)
  {}

private:

  const SubProblem& _subProblem;
  SubProblemLocalFunctionSpace _splfs_s;
  SubProblemLocalFunctionSpace _splfs_n;

};

template<typename MDLFS, typename SubProblem, typename Parameters>
struct BuildConstraintsDescriptor<MDLFS,ConstrainSubProblem<SubProblem,Parameters> >
{
  typedef SubProblemConstraints<MDLFS,SubProblem,Parameters> type;
};

template<typename MDLFS, typename SubProblem, typename Parameters>
shared_ptr<SubProblemConstraints<MDLFS,SubProblem,Parameters> >
buildConstraintsDescriptor(const MDLFS& mdlfs_s, const MDLFS& mdlfs_n,
                           const ConstrainSubProblem<SubProblem,Parameters>& d)
{
  return make_shared<SubProblemConstraints<MDLFS,SubProblem,Parameters> >(mdlfs_s,
                                                                          mdlfs_n,
                                                                          d.subProblem(),
                                                                          d.parameters());
}


template<typename CA, typename... ConstraintsDescriptors>
class ConstraintsAssemblerEngine
  : public Dune::PDELab::TypeTree::VariadicCompositeNode<ConstraintsDescriptors...>
  , public Dune::PDELab::LocalAssemblerEngineBase
{

  typedef Dune::PDELab::TypeTree::VariadicCompositeNode<ConstraintsDescriptors...> NodeT;
  typedef typename CA::CG CG;

public:

  bool requireSkeleton() const
  {
    return
      requireVSkeleton() ||
      requireVBoundary() ||
      requireVProcessor();
  }

  bool requireSkeletonTwoSided() const
  {
    return requireVBoundary();
  }

  bool requireVVolume() const
  {
    return any_child<ConstraintsAssemblerEngine,do_constraints_volume<> >::value;
  }

  bool requireVSkeleton() const
  {
    return
      any_child<ConstraintsAssemblerEngine,do_constraints_skeleton<> >::value ||
      any_child<ConstraintsAssemblerEngine,do_constraints_boundary<> >::value;
  }

  bool requireVBoundary() const
  {
    return any_child<ConstraintsAssemblerEngine,do_constraints_boundary<> >::value;
  }

  bool requireVProcessor() const
  {
    return any_child<ConstraintsAssemblerEngine,do_constraints_processor<> >::value;
  }


  template<typename EG, typename LFSV>
  void onBindLFSV(const EG& eg, const LFSV& lfsv)
  {
    applyVisitor(visitor<functors::rebind_subproblem_lfs_s>::add_data());
  }

  template<typename IG, typename LFSV_N>
  void onBindLFSVOutside(const IG& ig, const LFSV_N& lfsv_n)
  {
    applyVisitor(visitor<functors::rebind_subproblem_lfs_n>::add_data());
  }


  template<typename EG, typename LFSV>
  void assembleVVolume(const EG& eg, const LFSV& lfsv)
  {
    typedef visitor<functors::volume_constraints,do_constraints_volume<> > Visitor;
    applyVisitor(Visitor::add_data(wrap_eg(eg),wrap_lfsv_s(lfsv),wrap_cg(*cg)));
  }

  template<typename IG, typename LFSV_S, typename LFSV_N>
  void assembleVSkeleton(const IG& ig,
                         const LFSV_S& lfsv_s,
                         const LFSV_N& lfsv_n)
  {
    typedef visitor<functors::skeleton_or_boundary_constraints,do_constraints_skeleton_or_boundary<> > Visitor;
    applyVisitor(Visitor::add_data(wrap_ig(ig),wrap_lfsv_s(lfsv_s),wrap_lfsv_n(lfsv_n),wrap_cg(*cg)));
  }

  template<typename IG, typename LFSV>
  void assembleVBoundary(const IG& ig, const LFSV& lfsv)
  {
    typedef visitor<functors::boundary_constraints,do_constraints_boundary<> > Visitor;
    applyVisitor(Visitor::add_data(wrap_ig(ig),wrap_lfsv_s(lfsv),wrap_cg(*cg)));
  }

  template<typename IG, typename LFSV>
  void assembleVProcessor(const IG& ig, const LFSV& lfsv)
  {
    typedef visitor<functors::processor_constraints,do_constraints_processor<> > Visitor;
    applyVisitor(Visitor::add_data(wrap_ig(ig),wrap_lfsv_s(lfsv),wrap_cg(*cg)));
  }


  void setConstraintsContainer(CG& cg_)
  {
    cg = &cg_;
  }

  ConstraintsAssemblerEngine(shared_ptr<ConstraintsDescriptors>... constraintsSpecifications)
    : NodeT(constraintsSpecifications...)
  {}

private:

  template<typename Visitor>
  void applyVisitor(Visitor&& visitor)
  {
    Dune::PDELab::TypeTree::applyToTree(*this,visitor);
  }

  CG* cg;

};

template<typename GFS, typename RF, typename... ConstraintsSpecifications>
class ConstraintsAssembler
{

public:

  typedef GlobalAssembler<GFS,GFS> Assembler;
  typedef ConstraintsAssemblerEngine<ConstraintsAssembler,
                                     typename BuildConstraintsDescriptor<typename Assembler::LFSV,
                                                                         ConstraintsSpecifications
                                                                         >::type...
                                     > Engine;
  typedef typename GFS::template ConstraintsContainer<RF>::Type CG;


  ConstraintsAssembler(const GFS& gfs,
                       const ConstraintsSpecifications&... specifications)
    : _assembler(gfs,gfs)
    , _engine(buildConstraintsDescriptor(_assembler.lfsv_s(),
                                         _assembler.lfsv_n(),
                                         specifications)...)
  {}

  void assemble(CG& cg)
  {
    _engine.setConstraintsContainer(cg);
    _assembler.assemble(_engine);
  }

private:

  Assembler _assembler;
  Engine _engine;

};


template<typename RF, typename GFS, typename... ConstraintsSpecifications>
ConstraintsAssembler<GFS,RF,ConstraintsSpecifications...>
constraints(const GFS& gfs, const ConstraintsSpecifications&... specifications)
{
  return ConstraintsAssembler<GFS,RF,ConstraintsSpecifications...>(gfs,specifications...);
}


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_CONSTRAINTS_HH
