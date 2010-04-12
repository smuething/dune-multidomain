#ifndef DUNE_MULTIDOMAIN_VARIADICCOMPOSITENODE_HH
#define DUNE_MULTIDOMAIN_VARIADICCOMPOSITENODE_HH

#include <tuple>

namespace Dune {

namespace PDELab {

namespace MultiDomain {

/**
 * TMP for deriving storage types in VariadicCompositeNode
 */
template<typename... OArgs>
struct replace;

template<typename OHead, typename... OTail>
struct replace<OHead,OTail...> {
  template<template<typename T> class Transform, typename... SArgs>
  struct with {
    typedef typename replace<OTail...>::template with<Transform,SArgs...,typename Transform<OHead>::type>::type type;
  };
};

/**
 * End of recursion - export the transformed argument list as a tuple
 */
template<>
struct replace<> {
  template<template<typename T> class Transform, typename... SArgs>
  struct with {
    typedef typename Transform<void>::template container<SArgs...>::type type;
  };
};

/**
 * wrapper to simplify calling the TMP
 */
template<template<typename T> class Transform, typename... Args>
struct transform {
  typedef typename replace<Args...>::template with<Transform>::type type;
};


/**
 * VariadicCompositeNode
 *
 */
template<typename P, typename... Children>
class VariadicCompositeNode
{

  template<typename T>
  struct forwarder {
    typedef typename P::template Storage<T>::Type type;

    template<typename... Args>
    struct container {
      typedef std::tuple<Args...> type;
    };
  };

  typedef std::tuple<Children...> OT;
  typedef typename transform<forwarder,Children...>::type ST;

public:

  enum { isLeaf = false };
  enum { isPower = false /**< */ };
  enum { isComposite = true /**< */ };
  enum { CHILDREN = sizeof...(Children) };

  VariadicCompositeNode (Children&&... c_) : c(P::convert(c_)...) {}

  template<int i>
  struct Child
  {
    typedef typename std::tuple_element<i,OT>::type Type;
  };

  template<int i>
  typename Child<i>::Type& getChild ()
  {
    return P::get(std::get<i>(c));
  }

  template<int i>
  const typename Child<i>::Type& getChild () const
  {
    return P::get(std::get<i>(c));
  }

private:
  ST c;
};

} // namespace MultiDomain

} // namespace PDELab

} // namespace Dune

#endif // DUNE_MULTIDOMAIN_VARIADICCOMPOSITENODE_HH
