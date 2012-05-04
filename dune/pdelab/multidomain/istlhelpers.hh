#ifndef DUNE_MULTIDOMAIN_ISTLHELPERS_HH
#define DUNE_MULTIDOMAIN_ISTLHELPERS_HH

#include <dune/pdelab/backend/istlvectorbackend.hh>

namespace Dune {

namespace PDELab {

  namespace MultiDomain {

    struct CouplingGridFunctionSpaceTag;

  }

template<typename E, typename Node>
struct extract_istl_vector_helper<E,Node,MultiDomain::CouplingGridFunctionSpaceTag>
{
  typedef BlockVector<FieldVector<E,Node::Traits::Backend::blockSize> > type;
};

} // namespace PDELab

} // namespace Dune

#endif // DUNE_MULTIDOMAIN_ISTLHELPERS_HH
