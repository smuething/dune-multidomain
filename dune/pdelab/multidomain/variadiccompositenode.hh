#ifndef DUNE_MULTIDOMAIN_VARIADICCOMPOSITENODE_HH
#define DUNE_MULTIDOMAIN_VARIADICCOMPOSITENODE_HH

#include <tuple>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


/**
 * VariadicCompositeNode
 *
 */
template<typename P, typename... Children>
class VariadicCompositeNode
{

  struct forwarder {

    template<typename T>
    struct transform {
      typedef typename P::template Storage<T>::Type type;
    };

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

  VariadicCompositeNode (Children&... c_) : c(P::convert(c_)...) {}

  VariadicCompositeNode () {}

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
