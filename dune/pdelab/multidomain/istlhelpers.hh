#ifndef DUNE_MULTIDOMAIN_ISTLHELPERS_HH
#define DUNE_MULTIDOMAIN_ISTLHELPERS_HH

#include <dune/pdelab/backend/istlvectorbackend.hh>

namespace Dune {

namespace PDELab {

  namespace MultiDomain {

    struct CouplingGridFunctionSpaceTag;

  }

  namespace istl {

    template<typename E, typename Node>
    struct vector_descriptor_helper<E,Node,MultiDomain::CouplingGridFunctionSpaceTag>
    {
      typedef leaf_vector_descriptor<E,typename Node::Traits::Backend> type;
    };

  }

} // namespace PDELab

} // namespace Dune

#endif // DUNE_MULTIDOMAIN_ISTLHELPERS_HH
