#include <dune/pdelab/multidomain/visitor.hh>
#include <iostream>

template<typename data_container>
struct Visitor
  : public Dune::PDELab::MultiDomain::data_accessor<data_container>
{

  typedef Dune::PDELab::MultiDomain::data_accessor<data_container> DataAccessor;
  typedef typename DataAccessor::Data Data;
  using DataAccessor::data;

  void foo()
  {
    typename Data::LFSU_S& lfsu = data().lfsu_s();
    std::cout << data().x() << " " << lfsu << std::endl;
  }

};

DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(X,x)
DUNE_PDELAB_MULTIDOMAIN_CREATE_DATA_WRAPPER(LFSU_S,lfsu_s)

int main()
{
  int x = 3;
  double lfsu = 3.4;
  Dune::PDELab::MultiDomain::visitor<Visitor,int>::add_data(wrap_x(x),wrap_lfsu_s(lfsu)).foo();
  return 0;
}
