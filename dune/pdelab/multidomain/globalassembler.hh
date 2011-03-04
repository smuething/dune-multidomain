// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_GLOBALASSEMBLER_HH
#define DUNE_PDELAB_MULTIDOMAIN_GLOBALASSEMBLER_HH

#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>

#include <dune/pdelab/multidomain/entitywrappers.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


template<typename GridOperator, bool nonoverlapping_mode = false>
class GlobalAssembler
{

  typedef typename GridOperator::Traits::TrialGridFunctionSpace GFSU;
  typedef typename GridOperator::Traits::TestGridFunctionSpace GFSV;

  typedef typename GFSU::Traits::GridViewType GV;
  typedef typename GV::IndexSet IndexSet;
  typedef typename GV::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;

  typedef LocalFunctionSpace<GFSU> LFSU;
  typedef LocalFunctionSpace<GFSV> LFSV;
  typedef CouplingLocalFunctionSpace<GFSU> CouplingLFSU;
  typedef CouplingLocalFunctionSpace<GFSV> CouplingLFSV;


public:

  template<typename LocalAssemblerEngine>
  void assemble(LocalAssemblerEngine& engine)
  {

    GV gv = gfsu.gridview();
    const IndexSet& is = gv.indexSet();

    engine.preAssembly();

    // make local function spaces

    typedef Dune::PDELab::MultiDomain::ElementWrapper<GV> ElementWrapper;
    typedef Dune::PDELab::MultiDomain::BoundaryIntersectionWrapper<GV> BoundaryIntersectionWrapper;
    typedef Dune::PDELab::MultiDomain::SkeletonIntersectionWrapper<GV> SkeletonIntersectionWrapper;

    MultiGeomUniqueIDMapper<GV> cell_mapper(gv);

    const ElementIterator endit = gv.template end<0>();

    // traverse grid view
    for (ElementIterator it = gv.template begin<0>();
         it!=endit; ++it)
      {

        typename GV::IndexSet::IndexType ids = cell_mapper.map(*it);

        // skip ghost and overlap
        if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
          continue;

        // bind local function spaces to element
        lfsu.bind(*it);
        lfsv.bind(*it);

        ElementWrapper elementWrapper(*it,is.subDomains(*it));

        engine.onBindLFSUV(elementWrapper,lfsu,lfsv);
        engine.loadCoefficientsLFSUInside(lfsu);
        engine.onBindLFSV(elementWrapper,lfsv);

        engine.assembleUVVolume(elementWrapper,lfsu,lfsv);
        engine.assembleVVolume(elementWrapper,lfsv);

        // skip if no intersection iterator is needed
        if (engine.requireIntersections())
          {
            // traverse intersections
            unsigned int intersection_index = 0;
            IntersectionIterator endiit = gv.iend(*it);
            for (IntersectionIterator iit = gv.ibegin(*it);
                 iit!=endiit; ++iit, ++intersection_index)
              {
                // skeleton term
                if (iit->neighbor())
                  {
                    const typename GV::IndexSet::IndexType idn = cell_mapper.map(*(iit->outside()));
                    if (ids < idn && !engine.requireIntersectionsTwoSided())
                      continue;

                    lfsun.bind(*(iit->outside()));
                    lfsvn.bind(*(iit->outside()));

                    SkeletonIntersectionWrapper skeletonIntersectionWrapper(*iit,intersection_index,
                                                                            elementWrapper,
                                                                            is.subDomains(*(iit->outside())));

                    engine.onBindLFSUVOutside(skeletonIntersectionWrapper,lfsun,lfsvn);
                    engine.loadCoefficientsLFSUOutside(lfsun);
                    engine.onBindLFSVOutside(skeletonIntersectionWrapper,lfsvn);

                    engine.assembleUVSkeleton(skeletonIntersectionWrapper,lfsu,lfsv,lfsun,lfsvn);
                    engine.assembleVSkeleton(skeletonIntersectionWrapper,lfsv,lfsvn);

                    bool unbindLFSVCoupling = false;
                    if (engine.requireVEnrichedCoupling() || engine.requireUVEnrichedCoupling())
                      {
                        unbindLFSVCoupling = true;
                        couplinglfsv.bind(*iit);
                        engine.onBindLFSVCoupling(skeletonIntersectionWrapper,couplinglfsv);
                      }

                    if (engine.requireVEnrichedCoupling())
                      {
                      engine.assembleVEnrichedCoupling(skeletonIntersectionWrapper,lfsv,lfsvn,couplinglfsv);
                      }

                    if (engine.requireUVEnrichedCoupling())
                      {
                        couplinglfsu.bind(*iit);
                        engine.onBindLFSUVCoupling(skeletonIntersectionWrapper,couplinglfsu,couplinglfsv);
                        engine.loadCoefficientsLFSUCoupling(couplinglfsu);
                        engine.assembleUVEnrichedCoupling(skeletonIntersectionWrapper,lfsu,lfsv,lfsun,lfsvn,couplinglfsu,couplinglfsv);
                        engine.onUnbindLFSUVCoupling(skeletonIntersectionWrapper,couplinglfsu,couplinglfsv);
                      }

                    if (unbindLFSVCoupling)
                      {
                        engine.onUnbindLFSVCoupling(skeletonIntersectionWrapper,couplinglfsv);
                      }

                    engine.onUnbindLFSUVOutside(skeletonIntersectionWrapper,lfsun,lfsvn);
                    engine.onUnbindLFSVOutside(skeletonIntersectionWrapper,lfsvn);
                  }

                // boundary term
                if (iit->boundary())
                  {
                    BoundaryIntersectionWrapper boundaryIntersectionWrapper(*iit,intersection_index,elementWrapper);
                    engine.assembleUVBoundary(boundaryIntersectionWrapper,lfsu,lfsv);
                    engine.assembleVBoundary(boundaryIntersectionWrapper,lfsv);
                  }
              }
          }

        engine.assembleUVVolumePostSkeleton(elementWrapper,lfsu,lfsv);
        engine.assembleVVolumePostSkeleton(elementWrapper,lfsv);

        engine.onUnbindLFSUV(elementWrapper,lfsu,lfsv);
        engine.onUnbindLFSV(elementWrapper,lfsv);

      }

    engine.postAssembly();
  }

  GlobalAssembler(const GFSU& gfsu_, const GFSV& gfsv_)
    : gfsu(gfsu_)
    , gfsv(gfsv_)
    , lfsu(gfsu)
    , lfsv(gfsv)
    , lfsun(gfsu)
    , lfsvn(gfsv)
    , couplinglfsu(gfsu)
    , couplinglfsv(gfsv)
  {}

  const GFSU& trialGridFunctionSpace() const
  {
    return gfsu;
  }

  const GFSV& testGridFunctionSpace() const
  {
    return gfsv;
  }

private:

  const GFSU& gfsu;
  const GFSV& gfsv;

  LFSU lfsu;
  LFSV lfsv;

  LFSU lfsun;
  LFSV lfsvn;

  CouplingLFSU couplinglfsu;
  CouplingLFSV couplinglfsv;

};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_GLOBALASSEMBLER_HH
