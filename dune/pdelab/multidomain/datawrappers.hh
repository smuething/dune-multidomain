// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_DATAWRAPPERS_HH
#define DUNE_PDELAB_MULTIDOMAIN_DATAWRAPPERS_HH

namespace Dune {

namespace PDELab {

namespace MultiDomain {

#ifndef DOXYGEN // nothing in this file is public API

#define DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(TYPE_NAME,VARIABLE_NAME) \
template<typename TYPE_NAME##_> \
struct VARIABLE_NAME##_wrapper \
{ \
\
  typedef TYPE_NAME##_ TYPE_NAME; \
\
  VARIABLE_NAME##_wrapper(TYPE_NAME & VARIABLE_NAME) \
    : _##VARIABLE_NAME(VARIABLE_NAME) \
  {} \
\
  TYPE_NAME & VARIABLE_NAME() \
  { \
    return _##VARIABLE_NAME; \
  } \
\
  const TYPE_NAME & VARIABLE_NAME() const \
  { \
    return _##VARIABLE_NAME; \
  } \
\
private: \
\
  TYPE_NAME & _##VARIABLE_NAME; \
\
}; \
\
template<typename TYPE_NAME> \
VARIABLE_NAME##_wrapper<TYPE_NAME> wrap_##VARIABLE_NAME(TYPE_NAME & VARIABLE_NAME) \
{ \
  return VARIABLE_NAME##_wrapper<TYPE_NAME>(VARIABLE_NAME); \
} \
// end DATA_WRAPPER macro


#define DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(TYPE_NAME,VARIABLE_NAME) \
template<typename TYPE_NAME##_> \
struct VARIABLE_NAME##_wrapper \
{ \
\
  typedef TYPE_NAME##_ TYPE_NAME; \
\
  VARIABLE_NAME##_wrapper(const TYPE_NAME & VARIABLE_NAME) \
    : _##VARIABLE_NAME(VARIABLE_NAME) \
  {} \
\
  const TYPE_NAME & VARIABLE_NAME() const \
  { \
    return _##VARIABLE_NAME; \
  } \
\
private: \
\
  const TYPE_NAME & _##VARIABLE_NAME; \
\
}; \
\
template<typename TYPE_NAME> \
VARIABLE_NAME##_wrapper<TYPE_NAME> wrap_##VARIABLE_NAME(const TYPE_NAME & VARIABLE_NAME) \
{ \
  return VARIABLE_NAME##_wrapper<TYPE_NAME>(VARIABLE_NAME); \
} \
// end CONST_DATA_WRAPPER macro


#define DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_CONTAINER(VARIABLE_NAME) \
template<typename T> \
struct VARIABLE_NAME##_container \
{ \
\
  VARIABLE_NAME##_container(const T & VARIABLE_NAME) \
    : _##VARIABLE_NAME(VARIABLE_NAME) \
  {} \
\
  T & VARIABLE_NAME() \
  { \
    return _##VARIABLE_NAME; \
  } \
\
  const T & VARIABLE_NAME() const \
  { \
    return _##VARIABLE_NAME; \
  } \
\
private: \
\
  T _##VARIABLE_NAME; \
\
}; \
\
template<typename T> \
VARIABLE_NAME##_container<T> store_##VARIABLE_NAME(const T & VARIABLE_NAME) \
{ \
  return VARIABLE_NAME##_container<T>(VARIABLE_NAME); \
} \
// end DATA_CONTAINER macro

DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(Operator,operator_type)
DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(EG,eg)
DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(IG,ig)
DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(LFSU_S_Cache,lfsu_s_cache)
DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(LFSV_S_Cache,lfsv_s_cache)
DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(LFSU_N_Cache,lfsu_n_cache)
DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(LFSV_N_Cache,lfsv_n_cache)
DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(LFSU_C_Cache,lfsu_c_cache)
DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(LFSV_C_Cache,lfsv_c_cache)
DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(X,x)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(R,r)
DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(X_S,x_s)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(R_S,r_s)
DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(X_N,x_n)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(R_N,r_n)
DUNE_PDELAB_MULTIDOMAIN_CREATE_CONST_DATA_WRAPPER(X_C,x_c)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(R_C,r_c)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(A,a)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(A_SS,a_ss)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(A_SN,a_sn)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(A_NS,a_ns)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(A_NN,a_nn)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(A_SC,a_sc)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(A_CS,a_cs)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(A_NC,a_nc)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(A_CN,a_cn)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(A_CC,a_cc)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(Pattern_SS,pattern_ss)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(Pattern_SN,pattern_sn)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(Pattern_SC,pattern_sc)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(Pattern_nn,pattern_nn)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(Pattern_NS,pattern_ns)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(Pattern_NC,pattern_nc)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(Pattern_CC,pattern_cc)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(Pattern_CS,pattern_cs)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(Pattern_CN,pattern_cn)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(CG,cg)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_CONTAINER(time)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_CONTAINER(weight)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_CONTAINER(dt)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_CONTAINER(stages)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_CONTAINER(stage)

#endif // DOXYGEN

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_DATAWRAPPERS_HH
