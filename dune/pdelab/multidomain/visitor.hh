#ifndef DUNE_PDELAB_MULTIDOMAIN_VISITOR_HH
#define DUNE_PDELAB_MULTIDOMAIN_VISITOR_HH

namespace Dune {
namespace PDELab {
namespace MultiDomain {

#include <dune/pdelab/common/typetree.hh>

template<typename data_container>
struct data_accessor
{

  data_container& data()
  {
    return static_cast<data_container&>(*this);
  }

  const data_container& data() const
  {
    return static_cast<const data_container&>(*this);
  }

};


struct match_all_condition
{
  template<typename T>
  struct test
  {
    static const bool value = true;
  };
};


template<template<typename data> class functor_template, typename condition = match_all_condition>
struct visitor
{

  template<typename... data_wrappers>
  struct data_enriched_visitor
    : public functor_template<data_enriched_visitor<data_wrappers...> >
    , public Dune::PDELab::TypeTree::TreeVisitor
    , public Dune::PDELab::TypeTree::DynamicTraversal
    , public data_wrappers...
  {

    data_enriched_visitor(data_wrappers... wrappers)
      : data_wrappers(wrappers)...
    {}

    template<typename Node, typename TreePath>
    typename enable_if<condition::template test<typename remove_reference<Node>::type>::value == true>::type
    leaf(Node&& node, TreePath treePath)
    {
      (*this)(std::forward<Node>(node));
    }

    template<typename Node, typename TreePath>
    typename enable_if<condition::template test<typename remove_reference<Node>::type>::value == false>::type
    leaf(const Node& node, TreePath treePath)
    {}

  };

  template<typename... data_wrappers>
  static data_enriched_visitor<data_wrappers...>
  add_data(data_wrappers... wrappers)
  {
    return data_enriched_visitor<data_wrappers...>(wrappers...);
  }

};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_VISITOR_HH
