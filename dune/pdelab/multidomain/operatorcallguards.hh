#ifndef DUNE_PDELAB_MULTIDOMAIN_OPERATORCALLGUARDS_HH
#define DUNE_PDELAB_MULTIDOMAIN_OPERATORCALLGUARDS_HH

namespace Dune {

namespace PDELab {

namespace guarded_call {


#define OPERATORCALLGUARD_DEFINE_GUARD(GUARD,METHOD) \
template<typename Target, typename... Args> \
inline typename enable_if<remove_reference<Target>::type::GUARD>::type \
METHOD(Target&& target, Args&&... args) \
{ \
  target.METHOD(std::forward<Args>(args)...); \
} \
\
template<typename Target, typename... Args> \
inline typename enable_if<!remove_reference<Target>::type::GUARD>::type \
METHOD(Target&& target, Args&&... args) \
{}

OPERATORCALLGUARD_DEFINE_GUARD(doPatternSkeleton,pattern_skeleton)
OPERATORCALLGUARD_DEFINE_GUARD(doPatternBoundary,pattern_boundary)

OPERATORCALLGUARD_DEFINE_GUARD(doAlphaSkeleton,alpha_skeleton)
OPERATORCALLGUARD_DEFINE_GUARD(doAlphaBoundary,alpha_boundary)

OPERATORCALLGUARD_DEFINE_GUARD(doLambdaSkeleton,lambda_skeleton)
OPERATORCALLGUARD_DEFINE_GUARD(doLambdaBoundary,lambda_boundary)

OPERATORCALLGUARD_DEFINE_GUARD(doAlphaSkeleton,jacobian_skeleton)
OPERATORCALLGUARD_DEFINE_GUARD(doAlphaBoundary,jacobian_boundary)


OPERATORCALLGUARD_DEFINE_GUARD(doPatternEnrichedCouplingToSubProblems,pattern_enriched_coupling_first)
OPERATORCALLGUARD_DEFINE_GUARD(doPatternEnrichedCouplingToSubProblems,pattern_enriched_coupling_second)
OPERATORCALLGUARD_DEFINE_GUARD(doPatternEnrichedCoupling,pattern_enriched_coupling)

OPERATORCALLGUARD_DEFINE_GUARD(doAlphaEnrichedCouplingToSubProblems,alpha_enriched_coupling_first)
OPERATORCALLGUARD_DEFINE_GUARD(doAlphaEnrichedCouplingToSubProblems,alpha_enriched_coupling_second)
OPERATORCALLGUARD_DEFINE_GUARD(doAlphaEnrichedCoupling,alpha_enriched_coupling)

OPERATORCALLGUARD_DEFINE_GUARD(doLambdaEnrichedCouplingToSubProblems,lambda_enriched_coupling_first)
OPERATORCALLGUARD_DEFINE_GUARD(doLambdaEnrichedCouplingToSubProblems,lambda_enriched_coupling_second)
OPERATORCALLGUARD_DEFINE_GUARD(doLambdaEnrichedCoupling,lambda_enriched_coupling)

OPERATORCALLGUARD_DEFINE_GUARD(doAlphaEnrichedCouplingToSubProblems,jacobian_enriched_coupling_first)
OPERATORCALLGUARD_DEFINE_GUARD(doAlphaEnrichedCouplingToSubProblems,jacobian_enriched_coupling_second)
OPERATORCALLGUARD_DEFINE_GUARD(doAlphaEnrichedCoupling,jacobian_enriched_coupling)


} // namespace guarded_call
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_OPERATORCALLGUARDS_HH
