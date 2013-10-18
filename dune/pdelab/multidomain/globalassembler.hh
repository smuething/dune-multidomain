// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_GLOBALASSEMBLER_HH
#define DUNE_PDELAB_MULTIDOMAIN_GLOBALASSEMBLER_HH

#include<dune/common/exceptions.hh>
#include<dune/geometry/type.hh>

#include <dune/pdelab/common/elementmapper.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/multidomain/couplinglocalfunctionspace.hh>
#include <dune/pdelab/multidomain/entitywrappers.hh>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


template<typename GFSU_, typename GFSV_, bool nonoverlapping_mode = false>
class GlobalAssembler
{

  typedef GFSU_ GFSU;
  typedef GFSV_ GFSV;

  typedef typename GFSU::Traits::GridViewType GV;
  typedef typename GV::IndexSet IndexSet;
  typedef typename GV::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;

  typedef Dune::PDELab::MultiDomain::ElementWrapper<GV> ElementWrapper;
  typedef Dune::PDELab::MultiDomain::BoundaryIntersectionWrapper<GV> BoundaryIntersectionWrapper;
  typedef Dune::PDELab::MultiDomain::SkeletonIntersectionWrapper<GV> SkeletonIntersectionWrapper;

public:

  typedef LocalFunctionSpace<GFSU> LFSU;
  typedef LocalFunctionSpace<GFSV> LFSV;
  typedef CouplingLocalFunctionSpace<GFSU> CouplingLFSU;
  typedef CouplingLocalFunctionSpace<GFSV> CouplingLFSV;

  template<typename LocalAssemblerEngine>
  void assemble(LocalAssemblerEngine& engine)
  {

    typedef typename LocalAssemblerEngine::Traits::Spaces Spaces;

    typedef typename LocalAssemblerEngine::Traits::TrialGridFunctionSpaceConstraints CU;
    typedef typename LocalAssemblerEngine::Traits::TestGridFunctionSpaceConstraints CV;

    const CU& cu = engine.trialConstraints();
    const CV& cv = engine.testConstraints();

    typedef typename Spaces::LFSU_Cache LFSUCache;
    typedef typename Spaces::LFSV_Cache LFSVCache;

    typedef typename Spaces::LFSU_C_Cache CouplingLFSUCache;
    typedef typename Spaces::LFSU_C_Cache CouplingLFSVCache;

    LFSUCache lfsu_s_cache(_lfsu_s,cu);
    LFSVCache lfsv_s_cache(_lfsv_s,cv);
    LFSUCache lfsu_n_cache(_lfsu_n,cu);
    LFSVCache lfsv_n_cache(_lfsv_n,cv);
    CouplingLFSUCache lfsu_c_cache(_lfsu_c,cu);
    CouplingLFSVCache lfsv_c_cache(_lfsv_c,cv);

    // extract relevant requirements
    const bool require_intersections = engine.requireSkeleton();
    const bool require_two_sided = engine.requireSkeletonTwoSided();
    const bool require_uv_volume = engine.requireUVVolume();
    const bool require_uv_skeleton = engine.requireUVSkeleton();
    const bool require_v_skeleton = engine.requireVSkeleton();
    const bool require_uv_boundary = engine.requireUVBoundary();
    const bool require_v_boundary = engine.requireVBoundary();
    const bool require_uv_processor = engine.requireUVProcessor();
    const bool require_v_processor = engine.requireVProcessor();
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


    GV gv = gfsu.gridView();
    const IndexSet& is = gv.indexSet();

    engine.preAssembly();

    ElementMapper<GV> cell_mapper(gv);

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
        _lfsv_s.bind(*it);
        lfsv_s_cache.update();

        engine.onBindLFSV(elementWrapper,lfsv_s_cache);

        if (require_lfsu_s)
          {
            // bind lfsu_s
            _lfsu_s.bind(*it);
            lfsu_s_cache.update();

            engine.onBindLFSUV(elementWrapper,lfsu_s_cache,lfsv_s_cache);
            engine.loadCoefficientsLFSUInside(lfsu_s_cache);
          }

        // volume assembly
        engine.assembleUVVolume(elementWrapper,lfsu_s_cache,lfsv_s_cache);
        engine.assembleVVolume(elementWrapper,lfsv_s_cache);

        // skip if intersections are not needed
        if (require_intersections)
          {
            // traverse intersections
            unsigned int intersection_index = 0;
            const IntersectionIterator endiit = gv.iend(*it);
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
                                                                                ids > idn);

                        // bind lfsv_n
                        _lfsv_n.bind(*(iit->outside()));
                        lfsv_n_cache.update();

                        engine.onBindLFSVOutside(skeletonIntersectionWrapper,
                                                 lfsv_s_cache,
                                                 lfsv_n_cache);

                        if (require_lfsu_n)
                          {
                            // bind lfsu_n
                            _lfsu_n.bind(*(iit->outside()));
                            lfsu_n_cache.update();

                            engine.onBindLFSUVOutside(skeletonIntersectionWrapper,
                                                      lfsu_s_cache,lfsv_s_cache,
                                                      lfsu_n_cache,lfsv_n_cache);
                            engine.loadCoefficientsLFSUOutside(lfsu_n_cache);
                          }

                        engine.assembleUVSkeleton(skeletonIntersectionWrapper,
                                                  lfsu_s_cache,lfsv_s_cache,
                                                  lfsu_n_cache,lfsv_n_cache);
                        engine.assembleVSkeleton(skeletonIntersectionWrapper,
                                                 lfsv_s_cache,lfsv_n_cache);

                        if (require_lfsv_coupling)
                          {
                            // bind lfsv_c
                            _lfsv_c.bind(*iit);
                            lfsv_c_cache.update();

                            engine.onBindLFSVCoupling(skeletonIntersectionWrapper,
                                                      lfsv_s_cache,
                                                      lfsv_n_cache,
                                                      lfsv_c_cache);

                            if (require_uv_enriched_coupling)
                              {
                                // bind_lfsu_c
                                _lfsu_c.bind(*iit);
                                lfsu_c_cache.update(),

                                engine.onBindLFSUVCoupling(skeletonIntersectionWrapper,
                                                           lfsu_s_cache,lfsv_s_cache,
                                                           lfsu_n_cache,lfsv_n_cache,
                                                           lfsu_c_cache,lfsv_c_cache);
                                engine.loadCoefficientsLFSUCoupling(lfsu_c_cache);
                              }

                            engine.assembleUVEnrichedCoupling(skeletonIntersectionWrapper,
                                                              lfsu_s_cache,lfsv_s_cache,
                                                              lfsu_n_cache,lfsv_n_cache,
                                                              lfsu_c_cache,lfsv_c_cache);
                            engine.assembleVEnrichedCoupling(skeletonIntersectionWrapper,
                                                             lfsv_s_cache,
                                                             lfsv_n_cache,
                                                             lfsv_c_cache);

                            engine.onUnbindLFSUVCoupling(skeletonIntersectionWrapper,
                                                         lfsu_s_cache,lfsv_s_cache,
                                                         lfsu_n_cache,lfsv_n_cache,
                                                         lfsu_c_cache,lfsv_c_cache);
                            engine.onUnbindLFSVCoupling(skeletonIntersectionWrapper,
                                                        lfsv_s_cache,
                                                        lfsv_n_cache,
                                                        lfsv_c_cache);
                          }

                        engine.onUnbindLFSUVOutside(skeletonIntersectionWrapper,
                                                    lfsu_s_cache,lfsv_s_cache,
                                                    lfsu_n_cache,lfsv_n_cache);
                        engine.onUnbindLFSVOutside(skeletonIntersectionWrapper,
                                                   lfsv_s_cache,
                                                   lfsv_n_cache);
                      }
                    break;


                  case IntersectionType::boundary:

                    if (require_uv_boundary || require_v_boundary)
                      {
                        BoundaryIntersectionWrapper boundaryIntersectionWrapper(*iit,intersection_index,elementWrapper);
                        engine.assembleUVBoundary(boundaryIntersectionWrapper,lfsu_s_cache,lfsv_s_cache);
                        engine.assembleVBoundary(boundaryIntersectionWrapper,lfsv_s_cache);
                      }
                    break;


                  case IntersectionType::processor:

                    if (require_uv_processor || require_v_processor)
                      {
                        BoundaryIntersectionWrapper boundaryIntersectionWrapper(*iit,intersection_index,elementWrapper);
                        engine.assembleUVProcessor(boundaryIntersectionWrapper,lfsu_s_cache,lfsv_s_cache);
                        engine.assembleVProcessor(boundaryIntersectionWrapper,lfsv_s_cache);
                      }
                    break;

                  } // switch
              } // loop
          }  // if

        engine.assembleUVVolumePostSkeleton(elementWrapper,lfsu_s_cache,lfsv_s_cache);
        engine.assembleVVolumePostSkeleton(elementWrapper,lfsv_s_cache);

        engine.onUnbindLFSUV(elementWrapper,lfsu_s_cache,lfsv_s_cache);
        engine.onUnbindLFSV(elementWrapper,lfsv_s_cache);

      }

    engine.postAssembly(gfsu,gfsv);
  }

  GlobalAssembler(const GFSU& gfsu_, const GFSV& gfsv_)
    : gfsu(gfsu_)
    , gfsv(gfsv_)
    , _lfsu_s(gfsu)
    , _lfsv_s(gfsv)
    , _lfsu_n(gfsu)
    , _lfsv_n(gfsv)
    , _lfsu_c(gfsu)
    , _lfsv_c(gfsv)
  {}

  const GFSU& trialGridFunctionSpace() const
  {
    return gfsu;
  }

  const GFSV& testGridFunctionSpace() const
  {
    return gfsv;
  }

  // ************************************************************************************
  // Access to internal local function spaces for creating SubProblemLocalFunctionSpaces
  // ************************************************************************************

  const LFSU& lfsu_s() const
  {
    return _lfsu_s;
  }

  const LFSV& lfsv_s() const
  {
    return _lfsv_s;
  }

  const LFSU& lfsu_n() const
  {
    return _lfsu_n;
  }

  const LFSV& lfsv_n() const
  {
    return _lfsv_n;
  }

  const LFSU& lfsu_c() const
  {
    return _lfsu_c;
  }

  const LFSV& lfsv_c() const
  {
    return _lfsv_c;
  }


private:

  const GFSU& gfsu;
  const GFSV& gfsv;

  LFSU _lfsu_s;
  LFSV _lfsv_s;

  LFSU _lfsu_n;
  LFSV _lfsv_n;

  CouplingLFSU _lfsu_c;
  CouplingLFSV _lfsv_c;

};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_GLOBALASSEMBLER_HH
