// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_POLICY_HH
#define DUNE_PDELAB_MULTIDOMAIN_POLICY_HH

namespace Dune {

namespace PDELab {

namespace MultiDomain {


struct NoConstraintsCachingPolicy
{

  static const bool cache_trial_constraints = false;
  static const bool cache_test_constraints = false;
  static const bool cache_coupling_trial_constraints = false;
  static const bool cache_coupling_test_constraints = false;

};


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_POLICY_HH
