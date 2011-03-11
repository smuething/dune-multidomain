#ifndef DUNE_PDELAB_MULTIDOMAIN_CONSTRAINTSFUNCTORS_HH
#define DUNE_PDELAB_MULTIDOMAIN_CONSTRAINTSFUNCTORS_HH

#include <dune/pdelab/constraints/constraints.hh>
#include <dune/pdelab/common/typetree.hh>
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
  struct skeleton_or_boundary_constraints
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Descriptor>
    void operator()(Descriptor& descriptor)
    {
      if (!descriptor.appliesTo(data().ig().insideElement()))
        return;
      if (descriptor.appliesTo(data().ig().outsideElement()) &&
          data().ig().oneSidedDirection()) // Constraints are always applied one-sided
        {
          typedef SkeletonConstraints<typename Data::IG,typename Data::CG> Visitor;
          Dune::PDELab::TypeTree::applyToTreePair(descriptor.lfs_s(),
                                                  descriptor.lfs_n(),
                                                  Visitor(data().ig(),data().cg()));
        }
      else
        {
          typedef BoundaryConstraints<typename Data::IG, typename Data::CG> Visitor;
          Dune::PDELab::TypeTree::applyToTreePair(descriptor.parameters(),
                                                  descriptor.lfs_s(),
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
      if (!descriptor.appliesTo(data().ig().insideElement()))
          return;
      typedef BoundaryConstraints<typename Data::IG,typename Data::CG> Visitor;
      Dune::PDELab::TypeTree::applyToTreePair(descriptor.parameters(),
                                              descriptor.lfs_s(),
                                              Visitor(data().ig(),data().cg()));
    }

  };

  template<typename Data>
  struct processor_constraints
    : public data_accessor<Data>
  {

    using data_accessor<Data>::data;

    template<typename Descriptor>
    void operator()(Descriptor& descriptor)
    {
      if (!descriptor.appliesTo(data().ig().insideElement()))
        return;
      typedef ProcessorConstraints<typename Data::IG, typename Data::CG> Visitor;
      Dune::PDELab::TypeTree::applyToTree(descriptor.lfs_s(),
                                          Visitor(data().ig(),data().cg()));
    }

  };

} // namespace functors

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_CONSTRAINTSFUNCTORS_HH
