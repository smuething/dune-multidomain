// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTIDOMAIN_GLOBALASSEMBLER_HH
#define DUNE_PDELAB_MULTIDOMAIN_GLOBALASSEMBLER_HH

#include<dune/common/exceptions.hh>
#include<dune/geometry/type.hh>

#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
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
        _lfsv_s.bind(*it);
        engine.onBindLFSV(elementWrapper,_lfsv_s);

        if (require_lfsu_s)
          {
            // bind lfsu_s
            _lfsu_s.bind(*it);
            engine.onBindLFSUV(elementWrapper,_lfsu_s,_lfsv_s);
            engine.loadCoefficientsLFSUInside(_lfsu_s);
          }

        // volume assembly
        engine.assembleUVVolume(elementWrapper,_lfsu_s,_lfsv_s);
        engine.assembleVVolume(elementWrapper,_lfsv_s);

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
                                                                                ids > idn);

                        // bind lfsv_n
                        _lfsv_n.bind(*(iit->outside()));
                        engine.onBindLFSVOutside(skeletonIntersectionWrapper,
                                                 _lfsv_s,
                                                 _lfsv_n);

                        if (require_lfsu_n)
                          {
                            // bind lfsu_n
                            _lfsu_n.bind(*(iit->outside()));
                            engine.onBindLFSUVOutside(skeletonIntersectionWrapper,
                                                      _lfsu_s,_lfsv_s,
                                                      _lfsu_n,_lfsv_n);
                            engine.loadCoefficientsLFSUOutside(_lfsu_n);
                          }

                        engine.assembleUVSkeleton(skeletonIntersectionWrapper,_lfsu_s,_lfsv_s,_lfsu_n,_lfsv_n);
                        engine.assembleVSkeleton(skeletonIntersectionWrapper,_lfsv_s,_lfsv_n);

                        if (require_lfsv_coupling)
                          {
                            // bind lfsv_c
                            _lfsv_c.bind(*iit);
                            engine.onBindLFSVCoupling(skeletonIntersectionWrapper,
                                                      _lfsv_s,
                                                      _lfsv_n,
                                                      _lfsv_c);

                            if (require_uv_enriched_coupling)
                              {
                                // bind_lfsu_c
                                _lfsu_c.bind(*iit);
                                engine.onBindLFSUVCoupling(skeletonIntersectionWrapper,
                                                           _lfsu_s,_lfsv_s,
                                                           _lfsu_n,_lfsv_n,
                                                           _lfsu_c,_lfsv_c);
                                engine.loadCoefficientsLFSUCoupling(_lfsu_c);
                              }

                            engine.assembleUVEnrichedCoupling(skeletonIntersectionWrapper,_lfsu_s,_lfsv_s,_lfsu_n,_lfsv_n,_lfsu_c,_lfsv_c);
                            engine.assembleVEnrichedCoupling(skeletonIntersectionWrapper,_lfsv_s,_lfsv_n,_lfsv_c);

                            engine.onUnbindLFSUVCoupling(skeletonIntersectionWrapper,
                                                         _lfsu_s,_lfsv_s,
                                                         _lfsu_n,_lfsv_n,
                                                         _lfsu_c,_lfsv_c);
                            engine.onUnbindLFSVCoupling(skeletonIntersectionWrapper,
                                                        _lfsv_s,
                                                        _lfsv_n,
                                                        _lfsv_c);
                          }

                        engine.onUnbindLFSUVOutside(skeletonIntersectionWrapper,
                                                    _lfsu_s,_lfsv_s,
                                                    _lfsu_n,_lfsv_n);
                        engine.onUnbindLFSVOutside(skeletonIntersectionWrapper,
                                                   _lfsv_s,
                                                   _lfsv_n);
                      }
                    break;


                  case IntersectionType::boundary:

                    if (require_uv_boundary || require_v_boundary)
                      {
                        BoundaryIntersectionWrapper boundaryIntersectionWrapper(*iit,intersection_index,elementWrapper);
                        engine.assembleUVBoundary(boundaryIntersectionWrapper,_lfsu_s,_lfsv_s);
                        engine.assembleVBoundary(boundaryIntersectionWrapper,_lfsv_s);
                      }
                    break;


                  case IntersectionType::processor:

                    if (require_uv_processor || require_v_processor)
                      {
                        BoundaryIntersectionWrapper boundaryIntersectionWrapper(*iit,intersection_index,elementWrapper);
                        engine.assembleUVProcessor(boundaryIntersectionWrapper,_lfsu_s,_lfsv_s);
                        engine.assembleVProcessor(boundaryIntersectionWrapper,_lfsv_s);
                      }
                    break;

                  } // switch
              } // loop
          }  // if

        engine.assembleUVVolumePostSkeleton(elementWrapper,_lfsu_s,_lfsv_s);
        engine.assembleVVolumePostSkeleton(elementWrapper,_lfsv_s);

        engine.onUnbindLFSUV(elementWrapper,_lfsu_s,_lfsv_s);
        engine.onUnbindLFSV(elementWrapper,_lfsv_s);

      }

    engine.postAssembly();
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
