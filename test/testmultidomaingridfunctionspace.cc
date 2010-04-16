#include "config.h"

#include <dune/grid/sgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/multidomain/subproblemgridfunctionspace.hh>

int main(int argc, char** argv) {
  const int dim = 2;
  typedef Dune::SGrid<dim,dim> BaseGrid;
  const int s[2] = {4,4};
  const double h[2] = {1.0,1.0};
  BaseGrid baseGrid(s,h);
  typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,4> > Grid;
  Grid grid(baseGrid,false);
  typedef Grid::SubDomainGrid SubDomainGrid;
  SubDomainGrid& sdg0 = grid.subDomain(0);
  SubDomainGrid& sdg1 = grid.subDomain(1);
  typedef Grid::ctype ctype;
  typedef Grid::LeafGridView MDGV;
  typedef SubDomainGrid::LeafGridView SDGV;
  MDGV mdgv = grid.leafView();
  SDGV sdgv0 = sdg0.leafView();
  SDGV sdgv1 = sdg1.leafView();
  grid.startSubDomainMarking();
  for (MDGV::Codim<0>::Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
    {
      if (it->geometry().center()[0] > 0.5)
        grid.addToSubDomain(0,*it);
      if (it->geometry().center()[1] > 0.5)
        grid.addToSubDomain(1,*it);
    }
  grid.preUpdateSubDomains();
  grid.updateSubDomains();
  grid.postUpdateSubDomains();
  typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM;
  FEM fem;
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<MDGV,FEM,CON,VBE> MDGFS;
  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM,CON,VBE> SDGFS;
  MDGFS mdgfs(mdgv,fem);
  SDGFS sdgfs0(sdgv0,fem);
  SDGFS sdgfs1(sdgv1,fem);

  typedef BaseGrid::LeafGridView BGV;
  BGV bgv = baseGrid.leafView();
  typedef Dune::PDELab::GridFunctionSpace<BGV,FEM,CON,VBE> BGFS;
  BGFS bgfs(bgv,fem);

  typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,MDGFS,SDGFS,SDGFS> MultiGFS;
  MultiGFS multigfs(grid,mdgfs,sdgfs0,sdgfs1);
  MultiGFS::LocalFunctionSpace lfs(multigfs);
  MDGV::Codim<0>::Iterator it = mdgv.begin<0>();
  lfs.bind(*it);
  lfs.debug();
  const MultiGFS::LocalFunctionSpace::Child<0>::Type& c0lfs = lfs.getChild<0>();
  const MultiGFS::LocalFunctionSpace::Child<1>::Type& c1lfs = lfs.getChild<1>();

  MDGFS::LocalFunctionSpace mdlfs(mdgfs);
  mdlfs.bind(*it);
  mdlfs.debug();
  typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::GridFunctionSpaceLexicographicMapper,MDGFS,MDGFS> CGFS;
  CGFS cgfs(mdgfs,mdgfs);
  CGFS::LocalFunctionSpace clfs(cgfs);
  clfs.bind(*it);
  clfs.debug();

  typedef MDGV::IndexSet::SubDomainSet SDS;

  SDS sds1, sds2, sds3;
  sds1.add(0);sds1.add(1);sds1.add(2);
  sds2.add(1);sds2.add(2);
  sds3.add(0);sds3.add(1);

  Dune::PDELab::MultiDomain::IncludesSubDomains<SDS> t1(1,2);
  Dune::PDELab::MultiDomain::EqualsSubDomains<SDS> t2(1);

  std::cout << "t1(sds1) == " << (t1(sds1) ? "true" : "false") << std::endl;
  std::cout << "t1(sds2) == " << (t1(sds2) ? "true" : "false") << std::endl;
  std::cout << "t1(sds3) == " << (t1(sds3) ? "true" : "false") << std::endl;

  std::cout << "t2(sds1) == " << (t2(sds1) ? "true" : "false") << std::endl;
  std::cout << "t2(sds2) == " << (t2(sds2) ? "true" : "false") << std::endl;
  std::cout << "t2(sds3) == " << (t2(sds3) ? "true" : "false") << std::endl;

  Dune::PDELab::MultiDomain::SubProblemGridFunctionSpace<MultiGFS,decltype(t2),0,2> spgfs(multigfs);

}
