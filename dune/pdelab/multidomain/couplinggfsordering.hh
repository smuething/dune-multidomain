// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTIDOMAIN_COUPLINGGFSORDERING_HH
#define DUNE_PDELAB_MULTIDOMAIN_COUPLINGGFSORDERING_HH

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/localorderingbase.hh>
#include <dune/pdelab/ordering/orderingbase.hh>
#include <dune/pdelab/ordering/directleaflocalordering.hh>
#include <dune/pdelab/multidomain/dofmapper.hh>
#include <dune/pdelab/multidomain/istlhelpers.hh>

namespace Dune {
  namespace PDELab {
    namespace MultiDomain {


    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{


    template<typename FEM, typename GV, typename Predicate, typename DI, typename CI>
    class DirectIntersectionLeafLocalOrdering
      : public TypeTree::LeafNode
    {

      template<typename>
      friend class LeafIntersectionGridViewOrdering;

    public:

      typedef LocalOrderingTraits<GV,DI,CI> Traits;

      void map_local_index(const typename Traits::SizeType geometry_type_index,
                           const typename Traits::SizeType entity_index,
                           typename Traits::TreeIndexView mi,
                           typename Traits::ContainerIndex& ci) const
      {
        DUNE_THROW(NotImplemented,"not implemented");
      }

      template<typename ItIn, typename ItOut>
      void map_indices(const ItIn begin, const ItIn end, ItOut out) const
      {
        // don't do anything - this is handled by the specialized GridViewOrdering
      }

      typename Traits::SizeType size(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index) const
      {
        typedef typename Traits::SizeType size_type;
        if (_gt_used[geometry_type_index])
          {
            const size_type index = _gt_entity_offsets[geometry_type_index] + entity_index;
            return _entity_dof_offsets[index+1] - _entity_dof_offsets[index];
          }
        else
          {
            return 0;
          }
      }

      typename Traits::SizeType size(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index, const typename Traits::SizeType child_index) const
      {
        DUNE_THROW(NotImplemented,"not implemented");
      }

      typename Traits::SizeType offset(const typename Traits::SizeType geometry_type_index, const typename Traits::SizeType entity_index, const typename Traits::SizeType child_index) const
      {
        assert(child_index == 0);
        return 0;
      }

      DirectIntersectionLeafLocalOrdering(const shared_ptr<const FEM>& fem, const GV& gv, const Predicate& predicate)
        : _predicate(predicate)
        , _fem(fem)
        , _gv(gv)
        , _container_blocked(false)
      {}

      const typename Traits::GridView& gridView() const
      {
        return _gv;
      }

      const FEM& finiteElementMap() const
      {
        return *_fem;
      }

    private:

      typedef FiniteElementInterfaceSwitch<
      typename FEM::Traits::FiniteElement
      > FESwitch;

      typedef typename GV::Intersection Intersection;

      void pre_collect_used_geometry_types_from_intersection()
      {
        typedef typename Traits::SizeType size_type;
        const size_type dim = Traits::GridView::dimension;

        _codim_used.assign(dim + 1,0);
        _gt_used.assign(GlobalGeometryTypeIndex::size(dim),false);
        _max_local_size = 0;
      }


      void collect_used_geometry_types_from_intersection(const Intersection& intersection)
      {
        if (!_predicate(intersection))
          return;

        FESwitch::setStore(_fe_store,_fem->find(intersection));

        const typename FESwitch::Coefficients& coeffs =
          FESwitch::coefficients(*_fe_store);

        _max_local_size = std::max(_max_local_size,coeffs.size());

        const ReferenceElement<typename Traits::GridView::ctype, Traits::GridView::dimension-1>& ref_el = ReferenceElements<typename Traits::GridView::ctype,Traits::GridView::dimension-1>::general(intersection.type());

        for (std::size_t i = 0; i < coeffs.size(); ++i)
          {
            const LocalKey& key = coeffs.localKey(i);
            GeometryType gt = ref_el.type(key.subEntity(),key.codim());
            _gt_used[GlobalGeometryTypeIndex::index(gt)] = true;
            _codim_used[key.codim()] = true;
          }
      }


      template<typename It>
      void allocate_entity_offset_vector(It it, const It end)
      {
        _gt_entity_offsets.assign(GlobalGeometryTypeIndex::size(GV::dimension) + 1,0);
        for (; it != end; ++it)
          {
            if (_gt_used[GlobalGeometryTypeIndex::index(*it)])
              _gt_entity_offsets[GlobalGeometryTypeIndex::index(*it) + 1] = _gv.indexSet().size(*it);
          }
        std::partial_sum(_gt_entity_offsets.begin(),_gt_entity_offsets.end(),_gt_entity_offsets.begin());
        _entity_dof_offsets.assign(_gt_entity_offsets.back() + 1,0);
      }


      void extract_per_entity_sizes_from_intersection(const Intersection& intersection)
      {
        if (!_predicate(intersection))
          return;

        FESwitch::setStore(_fe_store,_fem->find(intersection));

        const typename FESwitch::Coefficients& coeffs =
          FESwitch::coefficients(*_fe_store);

        typedef typename Traits::SizeType size_type;

        const ReferenceElement<typename Traits::GridView::ctype,Traits::GridView::dimension-1>& ref_el =
          ReferenceElements<typename Traits::GridView::ctype,Traits::GridView::dimension-1>::general(intersection.type());

        DOFMapper<GV> dm(intersection);

        for (std::size_t i = 0; i < coeffs.size(); ++i)
          {
            const LocalKey& key = coeffs.localKey(i);
            GeometryType gt = ref_el.type(key.subEntity(),key.codim());
            const size_type geometry_type_index = GlobalGeometryTypeIndex::index(gt);

            const size_type entity_index = _gv.indexSet().subIndex(dm.element(),
                                                                   dm.mapSubIndex(key.subEntity(),key.codim()),
                                                                   key.codim() + 1
                                                                   );
            const size_type index = _gt_entity_offsets[geometry_type_index] + entity_index;

            _entity_dof_offsets[index+1] = std::max(_entity_dof_offsets[index+1],static_cast<size_type>(key.index() + 1));
          }

      }


      void finalize_non_fixed_size_update()
      {
        // convert per-entity sizes to offsets
        std::partial_sum(_entity_dof_offsets.begin(),_entity_dof_offsets.end(),_entity_dof_offsets.begin());
      }


      typename Traits::SizeType maxLocalSize() const
      {
        return _max_local_size;
      }

    private:

      const Predicate& _predicate;

    protected:

      shared_ptr<const FEM> _fem;
      typename FESwitch::Store _fe_store;

      GV _gv;
      typename Traits::SizeType _max_local_size;
      const bool _container_blocked;

      std::vector<bool> _codim_used;
      std::vector<bool> _gt_used;

      std::vector<typename Traits::SizeType> _gt_entity_offsets;
      std::vector<typename Traits::SizeType> _entity_dof_offsets;
    };


    template<typename LocalOrdering>
    class LeafIntersectionGridViewOrdering
      : public TypeTree::CompositeNode<LocalOrdering>
      , public VirtualOrderingBase<typename LocalOrdering::Traits::DOFIndex,
                                   typename LocalOrdering::Traits::ContainerIndex>
      , public OrderingBase<typename LocalOrdering::Traits::DOFIndex,
                            typename LocalOrdering::Traits::ContainerIndex>
    {
    public:
      typedef typename LocalOrdering::Traits Traits;

      static const bool has_dynamic_ordering_children = false;

      static const bool consume_tree_index = false;

    private:

      typedef typename Traits::GridView GV;

      typedef TypeTree::CompositeNode<LocalOrdering> NodeT;

      typedef OrderingBase<
        typename LocalOrdering::Traits::DOFIndex,
        typename LocalOrdering::Traits::ContainerIndex
        > BaseT;

    public:

      LocalOrdering& localOrdering()
      {
        return this->template child<0>();
      }

      const LocalOrdering& localOrdering() const
      {
        return this->template child<0>();
      }


      LeafIntersectionGridViewOrdering(const typename NodeT::NodeStorage& localOrdering, bool container_blocked)
        : NodeT(localOrdering)
        , BaseT(*this,container_blocked,this)
        , _gv(this->template child<0>().gridView())
      {}

      virtual void map_index_dynamic(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        mapIndex(di,ci);
      }

      typename Traits::ContainerIndex mapIndex(const typename Traits::DOFIndex& di) const
      {
        typename Traits::ContainerIndex ci;
        mapIndex(di.view(),ci);
        return ci;
      }

      void mapIndex(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {

        const typename Traits::SizeType geometry_type_index = Traits::DOFIndexAccessor::geometryType(di);
        const typename Traits::SizeType entity_index = Traits::DOFIndexAccessor::entityIndex(di);
        assert (di.treeIndex().size() == 1);
        ci.push_back(di.treeIndex().back());

        ci.back() += localOrdering()._entity_dof_offsets[localOrdering()._gt_entity_offsets[geometry_type_index] + entity_index];

      }


      template<typename ItIn, typename ItOut>
      void map_lfs_indices(const ItIn begin, const ItIn end, ItOut out) const
      {
        typedef typename Traits::SizeType size_type;

        for (ItIn in = begin; in != end; ++in, ++out)
          {
            assert(in->treeIndex().size() == 1);
            out->push_back(in->treeIndex().back());
            const size_type geometry_type_index = Traits::DOFIndexAccessor::geometryType(*in);
            const size_type entity_index = Traits::DOFIndexAccessor::entityIndex(*in);
            out->back() += localOrdering()._entity_dof_offsets[localOrdering()._gt_entity_offsets[geometry_type_index] + entity_index];
          }
      }


      void update()
      {
        LocalOrdering& lo = localOrdering();

        const std::size_t dim = GV::dimension;

        typedef typename Traits::SizeType size_type;
        typedef std::vector<GeometryType> GTVector;
        GTVector geom_types;

        for (size_type cc = 0; cc <= dim; ++cc)
          {
            const GTVector& per_codim_geom_types = _gv.indexSet().geomTypes(cc);
            std::copy(per_codim_geom_types.begin(),per_codim_geom_types.end(),std::back_inserter(geom_types));
          }

        lo.pre_collect_used_geometry_types_from_intersection();

        typedef typename GV::template Codim<0>::Iterator CellIterator;

        typedef typename GV::IntersectionIterator IntersectionIterator;

        const CellIterator end_it = _gv.template end<0>();
        for (CellIterator it = _gv.template begin<0>(); it != end_it; ++it)
          {
            for (IntersectionIterator iit = _gv.ibegin(*it),
                   end_iit = _gv.iend(*it);
                 iit != end_iit;
                 ++iit)
              {
                lo.collect_used_geometry_types_from_intersection(*iit);
              }
          }

        lo.allocate_entity_offset_vector(geom_types.begin(),geom_types.end());

        for (CellIterator it = _gv.template begin<0>(); it != end_it; ++it)
          {
            for (IntersectionIterator iit = _gv.ibegin(*it),
                   end_iit = _gv.iend(*it);
                 iit != end_iit;
                 ++iit)
              {
                lo.extract_per_entity_sizes_from_intersection(*iit);
              }
          }

        lo.finalize_non_fixed_size_update();

        _block_count = _size = lo._entity_dof_offsets.back();

        _fixed_size = false;
        _max_local_size = lo.maxLocalSize();

      }

      //! dofs are blocked per entity/intersection on the leafs
      bool blocked() const { return true; }

      //! \brief whether all entites of the same geometry type/all
      //!        intersections have the same number of dofs
      /**
       * On the leaf this is realized by iterating over the grid during update
       * an checking.
       *
       * \note Even if fixedSize()==true the number of dofs may still vary
       *       between entities od different geometry type or between entities
       *       and intersections.
       */
      bool fixedSize() const { return false; }

      //! \brief maximum number of dofs attached to any given element and all
      //!        of its subentities and intersections
      /**
       * This is generally not an exact maximum and may be bigger than the
       * actual maximum.  There is however one special case: it is guaranteed
       * to be the exact maximum for fixedSize()==true.
       */

    private:

      using BaseT::_max_local_size;
      using BaseT::_size;
      using BaseT::_block_count;
      using BaseT::_container_blocked;
      using BaseT::_fixed_size;

      typename Traits::GridView _gv;

    };


    template<typename GFS, typename Transformation>
    struct leaf_coupling_gfs_to_ordering_descriptor
    {

      static const bool recursive = false;

      typedef DirectIntersectionLeafLocalOrdering<
        typename GFS::Traits::FiniteElementMap,
        typename GFS::Traits::GridView,
        typename GFS::Predicate,
        typename Transformation::DOFIndex,
        typename Transformation::ContainerIndex
        > LocalOrdering;

      typedef LeafIntersectionGridViewOrdering<LocalOrdering> GridViewOrdering;

      typedef GridViewOrdering transformed_type;
      typedef shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const GFS& gfs, const Transformation& t)
      {
        return transformed_type(make_tuple(make_shared<LocalOrdering>(gfs.finiteElementMapStorage(),gfs.gridView(),gfs.predicate())),gfs.backend().blocked());
      }

      static transformed_storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t)
      {
        return make_shared<transformed_type>(make_tuple(make_shared<LocalOrdering>(gfs->finiteElementMapStorage(),gfs->gridView(),gfs->predicate())),gfs->backend().blocked());
      }

    };


    template<typename GridFunctionSpace, typename Params>
    leaf_coupling_gfs_to_ordering_descriptor<
      GridFunctionSpace,
      gfs_to_ordering<Params>
      >
    registerNodeTransformation(GridFunctionSpace* gfs, gfs_to_ordering<Params>* t, CouplingGridFunctionSpaceTag* tag);



   //! \} group GridFunctionSpace
    } // namespace MultiDomain
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_COUPLINGGFSORDERING_HH
