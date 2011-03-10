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

  typedef Dune::PDELab::MultiDomain::ElementWrapper<GV> ElementWrapper;
  typedef Dune::PDELab::MultiDomain::BoundaryIntersectionWrapper<GV> BoundaryIntersectionWrapper;
  typedef Dune::PDELab::MultiDomain::SkeletonIntersectionWrapper<GV> SkeletonIntersectionWrapper;

public:

  template<typename LocalAssemblerEngine>
  void assemble(LocalAssemblerEngine& engine)
  {

    // extract relevant requirements
    const bool require_intersections = engine.requireSkeleton();
    const bool require_two_sided = engine.requireSkeletonTwoSided();
    const bool require_uv_volume = engine.requireUVVolume();
    const bool require_uv_skeleton = engine.requireUVSkeleton();
    const bool require_v_skeleton = engine.requireVSkeleton();
    const bool require_uv_boundary = engine.requireUVBoundary();
    const bool require_v_boundary = engine.requireVBoundary();
    const bool require_uv_processor = engine.requireUVProcessor();
    const bool require_uv_enriched_coupling = engine.requireUVEnrichedCoupling();
    const bool require_v_enriched_coupling = engine.requireVEnrichedCoupling();
    const bool require_uv_volume_post_skeleton = engine.requireUVVolumePostSkeleton();

    // do we need to bind the trial function space in the neighbor?
    const bool require_lfsu_n =
      require_uv_skeleton ||
      require_uv_enriched_coupling;

    // do we need to bind the coupling test function space?
    const bool require_lfsv_coupling =
      require_uv_enriched_coupling ||
      require_v_enriched_coupling;

    // do we need to bind the local trial function space?
    const bool require_lfsu_s =
      require_uv_volume ||
      require_uv_boundary ||
      require_uv_processor ||
      require_uv_volume_post_skeleton ||
      require_lfsu_n ||
      require_uv_enriched_coupling;

    // do we need to visit skeleton intersections?
    const bool require_skeleton =
      require_uv_skeleton ||
      require_v_skeleton ||
      require_uv_boundary ||
      require_v_boundary ||
      require_uv_enriched_coupling ||
      require_v_enriched_coupling;


    GV gv = gfsu.gridview();
    const IndexSet& is = gv.indexSet();

    engine.preAssembly();

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

        ElementWrapper elementWrapper(*it,is.subDomains(*it));

        if (engine.assembleCell(elementWrapper))
          continue;

        // bind lfsv_s
        lfsv_s.bind(*it);
        engine.onBindLFSV(elementWrapper,lfsv_s);

        if (require_lfsu_s)
          {
            // bind lfsu_s
            lfsu_s.bind(*it);
            engine.onBindLFSUV(elementWrapper,lfsu_s,lfsv_s);
            engine.loadCoefficientsLFSUInside(lfsu_s);
          }

        // volume assembly
        engine.assembleUVVolume(elementWrapper,lfsu_s,lfsv_s);
        engine.assembleVVolume(elementWrapper,lfsv_s);

        // skip if intersections are not needed
        if (require_intersections)
          {
            // traverse intersections
            unsigned int intersection_index = 0;
            IntersectionIterator endiit = gv.iend(*it);
            for (IntersectionIterator iit = gv.ibegin(*it);
                 iit!=endiit; ++iit, ++intersection_index)
              {
                switch (IntersectionType::get(*iit))
                  {

                  case IntersectionType::skeleton:
                  case IntersectionType::periodic:

                    if (require_skeleton)
                      {
                        const typename GV::IndexSet::IndexType idn = cell_mapper.map(*(iit->outside()));

                        if (!require_two_sided && ids < idn)
                          break;

                        SkeletonIntersectionWrapper skeletonIntersectionWrapper(*iit,intersection_index,
                                                                                elementWrapper,
                                                                                is.subDomains(*(iit->outside())),
                                                                                ids < idn);

                        // bind lfsv_n
                        lfsv_n.bind(*(iit->outside()));
                        engine.onBindLFSVOutside(skeletonIntersectionWrapper,lfsv_n);

                        if (require_lfsu_n)
                          {
                            // bind lfsu_n
                            lfsu_n.bind(*it);
                            engine.onBindLFSUVOutside(skeletonIntersectionWrapper,lfsu_n,lfsv_n);
                            engine.loadCoefficientsLFSUOutside(lfsu_n);
                          }

                        engine.assembleUVSkeleton(skeletonIntersectionWrapper,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
                        engine.assembleVSkeleton(skeletonIntersectionWrapper,lfsv_s,lfsv_n);

                        if (require_lfsv_coupling)
                          {
                            // bind lfsv_c
                            lfsv_c.bind(*iit);
                            engine.onBindLFSVCoupling(skeletonIntersectionWrapper,lfsv_c);

                            if (require_uv_enriched_coupling)
                              {
                                // bind_lfsu_c
                                lfsu_c.bind(*iit);
                                engine.onBindLFSUVCoupling(skeletonIntersectionWrapper,lfsu_c,lfsv_c);
                                engine.loadCoefficientsLFSUCoupling(lfsu_c);
                              }

                            engine.assembleUVEnrichedCoupling(skeletonIntersectionWrapper,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_c,lfsv_c);
                            engine.assembleVEnrichedCoupling(skeletonIntersectionWrapper,lfsv_c,lfsv_n,lfsv_c);

                            engine.onUnbindLFSUVCoupling(skeletonIntersectionWrapper,lfsu_c,lfsv_c);
                            engine.onUnbindLFSVCoupling(skeletonIntersectionWrapper,lfsv_c);
                          }

                        engine.onUnbindLFSUVOutside(skeletonIntersectionWrapper,lfsu_n,lfsv_n);
                        engine.onUnbindLFSVOutside(skeletonIntersectionWrapper,lfsv_n);
                      }
                    break;


                  case IntersectionType::boundary:

                    if (require_uv_boundary)
                      {
                        BoundaryIntersectionWrapper boundaryIntersectionWrapper(*iit,intersection_index,elementWrapper);
                        engine.assembleUVBoundary(boundaryIntersectionWrapper,lfsu_s,lfsv_s);
                        engine.assembleVBoundary(boundaryIntersectionWrapper,lfsv_s);
                      }
                    break;


                  case IntersectionType::processor:

                    if (require_uv_processor)
                      {
                        BoundaryIntersectionWrapper boundaryIntersectionWrapper(*iit,intersection_index,elementWrapper);
                        engine.assembleUVProcessor(boundaryIntersectionWrapper,lfsu_s,lfsv_s);
                        engine.assembleVProcessor(boundaryIntersectionWrapper,lfsv_s);
                      }
                    break;

                  } // switch
              } // loop
          }  // if

        engine.assembleUVVolumePostSkeleton(elementWrapper,lfsu_s,lfsv_s);
        engine.assembleVVolumePostSkeleton(elementWrapper,lfsv_s);

        engine.onUnbindLFSUV(elementWrapper,lfsu_s,lfsv_s);
        engine.onUnbindLFSV(elementWrapper,lfsv_s);

      }

    engine.postAssembly();
  }

  GlobalAssembler(const GFSU& gfsu_, const GFSV& gfsv_)
    : gfsu(gfsu_)
    , gfsv(gfsv_)
    , lfsu_s(gfsu)
    , lfsv_s(gfsv)
    , lfsu_n(gfsu)
    , lfsv_n(gfsv)
    , lfsu_c(gfsu)
    , lfsv_c(gfsv)
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

  LFSU lfsu_s;
  LFSV lfsv_s;

  LFSU lfsu_n;
  LFSV lfsv_n;

  CouplingLFSU lfsu_c;
  CouplingLFSV lfsv_c;

};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_GLOBALASSEMBLER_HH
