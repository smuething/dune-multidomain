#ifndef DUNE_PDELAB_MULTIDOMAIN_VISITOR_HH
#define DUNE_PDELAB_MULTIDOMAIN_VISITOR_HH

namespace Dune {
namespace PDELab {
namespace MultiDomain {


template<typename data_container>
struct data_accessor
{

  typedef data_container Data;

  data_container& data()
  {
    return static_cast<data_container&>(*this);
  }

  const data_container& data() const
  {
    return static_cast<const data_container&>(*this);
  }

};


template<template<typename data> class visitor_template, typename condition>
struct visitor
{

  template<typename... data_wrappers>
  struct data_enriched_visitor
    : public visitor_template<data_enriched_visitor<data_wrappers...> >
    , public data_wrappers...
  {
    data_enriched_visitor(data_wrappers... wrappers)
      : data_wrappers(wrappers)...
    {}
  };

  template<typename... data_wrappers>
  static data_enriched_visitor<data_wrappers...>
  add_data(data_wrappers... wrappers)
  {
    return data_enriched_visitor<data_wrappers...>(wrappers...);
  }

};

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
  return VARIABLE_NAME##_wrapper<T>(VARIABLE_NAME); \
} \
// end DATA_CONTAINER macro

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_VISITOR_HH