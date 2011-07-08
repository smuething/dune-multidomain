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

} // namespace guarded_call
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_OPERATORCALLGUARDS_HH
