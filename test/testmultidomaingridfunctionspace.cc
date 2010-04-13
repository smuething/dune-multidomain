#include "config.h"

#include <dune/grid/sgrid.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>

int main(int argc, char** argv) {
  const int dim = 2;
  typedef Dune::SGrid<dim,dim> BaseGrid;
  const int s[2] = {4,4};
  const double h[2] = {1.0,1.0};
  BaseGrid baseGrid(s,h);
  typedef Dune::MultiDomainGrid<BaseGrid,Dune::mdgrid::FewSubDomainsTraits<BaseGrid::dimension,2> > Grid;
  Grid grid(baseGrid,false);
  typedef Grid::SubDomainGrid SubDomainGrid;
  SubDomainGrid& sdg = grid.subDomain(0);
  typedef Grid::ctype ctype;
  typedef Grid::LeafGridView MDGV;
  typedef SubDomainGrid::LeafGridView SDGV;
  MDGV mdgv = grid.leafView();
  SDGV sdgv = sdg.leafView();
  typedef Dune::PDELab::Q1LocalFiniteElementMap<ctype,double,dim> FEM;
  FEM fem;
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<MDGV,FEM,CON,VBE> MDGFS;
  typedef Dune::PDELab::GridFunctionSpace<SDGV,FEM,CON,VBE> SDGFS;
  MDGFS mdgfs(mdgv,fem);
  SDGFS sdgfs(sdgv,fem);

  typedef BaseGrid::LeafGridView BGV;
  BGV bgv = baseGrid.leafView();
  typedef Dune::PDELab::GridFunctionSpace<BGV,FEM,CON,VBE> BGFS;
  BGFS bgfs(bgv,fem);

  typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<Grid,MDGFS,SDGFS/*,BGFS*/> MultiGFS;
  MultiGFS multigfs(grid,mdgfs,sdgfs/*,bgfs*/);
}
