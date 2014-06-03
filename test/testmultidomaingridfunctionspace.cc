#include "config.h"

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/sgrid.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>

int main(int argc, char** argv) {

  try {

    Dune::MPIHelper::instance(argc,argv);

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
    MDGV mdgv = grid.leafGridView();
    SDGV sdgv0 = sdg0.leafGridView();
    SDGV sdgv1 = sdg1.leafGridView();
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
    typedef Dune::PDELab::QkLocalFiniteElementMap<MDGV,ctype,double,1> MDFEM;
    typedef Dune::PDELab::QkLocalFiniteElementMap<SDGV,ctype,double,1> SDFEM;
    MDFEM mdfem(mdgv);
    SDFEM sdfem0(sdgv0);
    SDFEM sdfem1(sdgv1);
    typedef Dune::PDELab::NoConstraints CON;
    typedef Dune::PDELab::ISTLVectorBackend<> VBE;
    typedef Dune::PDELab::GridFunctionSpace<MDGV,MDFEM,CON,VBE> MDGFS;
    typedef Dune::PDELab::GridFunctionSpace<SDGV,SDFEM,CON,VBE> SDGFS;
    MDGFS mdgfs(mdgv,mdfem);
    SDGFS sdgfs0(sdgv0,sdfem0);
    SDGFS sdgfs1(sdgv1,sdfem1);

    typedef BaseGrid::LeafGridView BGV;
    BGV bgv = baseGrid.leafGridView();
    typedef Dune::PDELab::QkLocalFiniteElementMap<BGV,ctype,double,1> BFEM;
    BFEM bfem(bgv);
    typedef Dune::PDELab::GridFunctionSpace<BGV,BFEM,CON,VBE> BGFS;
    BGFS bgfs(bgv,bfem);

    typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<
      Grid,
      VBE,
      Dune::PDELab::LexicographicOrderingTag,
      MDGFS,
      SDGFS,
      SDGFS
      > MultiGFS;
    MultiGFS multigfs(grid,mdgfs,sdgfs0,sdgfs1);
    Dune::PDELab::LocalFunctionSpace<MultiGFS> lfs(multigfs);
    MDGV::Codim<0>::Iterator it = mdgv.begin<0>();
    lfs.bind(*it);
    const Dune::PDELab::LocalFunctionSpace<MultiGFS>::Child<0>::Type& c0lfs = lfs.child<0>();
    const Dune::PDELab::LocalFunctionSpace<MultiGFS>::Child<1>::Type& c1lfs = lfs.child<1>();

    std::cout << c0lfs.size() << " " << c1lfs.size() << std::endl;

    MDGFS mdgfs1(mdgv,mdfem);
    MDGFS mdgfs2(mdgv,mdfem);
    typedef Dune::PDELab::CompositeGridFunctionSpace<
      VBE,
      Dune::PDELab::LexicographicOrderingTag,
      MDGFS,MDGFS
      > CGFS;
    CGFS cgfs(mdgfs1,mdgfs2);
    Dune::PDELab::LocalFunctionSpace<CGFS> clfs(cgfs);
    clfs.bind(*it);

    typedef MDGV::IndexSet::SubDomainSet SDS;

    SDS sds1, sds2, sds3;
    sds1.add(0);sds1.add(1);sds1.add(2);
    sds2.add(1);sds2.add(2);
    sds3.add(0);sds3.add(1);

    // Dune::PDELab::MultiDomain::IncludesSubDomains<SDS> t1(1,2);
    // Dune::PDELab::MultiDomain::EqualsSubDomains<SDS> t2(0);

    // std::cout << "t1(sds1) == " << (t1(sds1) ? "true" : "false") << std::endl;
    // std::cout << "t1(sds2) == " << (t1(sds2) ? "true" : "false") << std::endl;
    // std::cout << "t1(sds3) == " << (t1(sds3) ? "true" : "false") << std::endl;

    // std::cout << "t2(sds1) == " << (t2(sds1) ? "true" : "false") << std::endl;
    // std::cout << "t2(sds2) == " << (t2(sds2) ? "true" : "false") << std::endl;
    // std::cout << "t2(sds3) == " << (t2(sds3) ? "true" : "false") << std::endl;

    // Dune::PDELab::MultiDomain::SubProblemGridFunctionSpace<MultiGFS,decltype(t2),0,2> spgfs(multigfs);

    // Dune::PDELab::MultiDomain::SubProblemLocalFunctionSpace<MultiGFS::LocalFunctionSpace,decltype(t2),0,1> splfs(lfs,t2);

  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }

}
