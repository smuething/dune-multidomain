#ifndef DUNE_PDELAB_MULTIDOMAIN_CONSTRAINTSFUNCTORS_HH
#define DUNE_PDELAB_MULTIDOMAIN_CONSTRAINTSFUNCTORS_HH

#include <dune/pdelab/constraints/constraints.hh>
#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/multidomain/operatorflagtests.hh>
#include <dune/pdelab/multidomain/datawrappers.hh>

namespace Dune {
namespace PDELab {
namespace MultiDomain {

namespace functors {

  template<typename Data>
  struct rebind_subproblem_lfs_s
    : public data_accessor<Data>
  {

    template<typename Descriptor>
    void operator()(Descriptor& descriptor)
    {
      descriptor.rebind_lfs_s();
    }

  };

  template<typename Data>
  struct rebind_subproblem_lfs_n
    : public data_accessor<Data>
  {

    template<typename Descriptor>
    void operator()(Descriptor& descriptor)
    {
      descriptor.rebind_lfs_n();
    }

  };

  template<typename Data>
  struct volume_constraints
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Descriptor>
    void operator()(Descriptor& descriptor)
    {
      if (!descriptor.appliesTo(data().eg()))
        return;
      typedef VolumeConstraints<typename Data::EG,typename Data::CG> Visitor;
      Dune::PDELab::TypeTree::applyToTree(descriptor.lfs_s(),Visitor(data().eg(),data().cg()));
    }

  };

  template<typename Data>
  struct skeleton_or_processor_or_boundary_constraints
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Descriptor>
    void operator()(Descriptor& descriptor)
    {
      if (!descriptor.appliesTo(data().ig().insideElement()))
        return;
      if (!descriptor.appliesTo(data().ig().outsideElement()))
        {
          // handle subproblem boundary
          typedef BoundaryConstraints<typename Data::IG,typename Data::CG> Visitor;
          Dune::PDELab::TypeTree::applyToTreePair(descriptor.parameters(),descriptor.lfs_s(),
                                                  Visitor(data().ig(),data().cg()));
        }
    }

  };


  template<typename Data>
  struct boundary_constraints
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Descriptor>
    void operator()(Descriptor& descriptor)
    {
      if (!descriptor.appliesTo(data().ig().insideElement())
          return;
      typedef BoundaryConstraints<typename Data::IG,typename Data::CG> Visitor;
      Dune::PDELab::TypeTree::applyToTreePair(descriptor.parameters(),descriptor.lfs_s(),
                                              Visitor(data().ig(),data().cg()));
    }

  };



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
    apply(subProblemConstraints,BindSubProblemLFS_S());
  }

  template<typename IG, typename LFSV_N>
  void onBindLFSVOutside(const IG& ig, const LFSV_N& lfsv_n)
  {
    typedef visitor<functors::bind:
    TypeTree::applyToTree(subProblemConstraints,bind_subproblem_lfs_n());
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

} // namespace functors

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_CONSTRAINTSFUNCTORS_HH
