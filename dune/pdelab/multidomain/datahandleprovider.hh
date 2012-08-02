// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_MULTIDOMAIN_DATAHANDLEPROVIDER_HH
#define DUNE_PDELAB_MULTIDOMAIN_DATAHANDLEPROVIDER_HH

#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>

namespace Dune {
  namespace PDELab {
    namespace MultiDomain {

    namespace {

      template<typename Entity, typename Traits>
      struct get_size_for_entity
        : public TypeTree::DirectChildrenVisitor
        , public TypeTree::StaticTraversal
      {

        typedef typename Traits::DOFIndex::EntityIndex EntityIndex;

        template<typename GFS, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(const GFS& gfs, const Child& child, TreePath tp, ChildIndex childIndex)
        {
          accumulate_size(child,
                          gfs.ordering().template child<ChildIndex::value>(),
                          typename gfs_flavor_tag<Child>::type());
        }


        template<typename GFS, typename Ordering>
        void accumulate_size(const GFS& gfs, const Ordering& ordering, MultiDomainGFSTag tag)
        {
          Dune::PDELab::get_size_for_entity<EntityIndex> get_size(_entity_index);
          TypeTree::applyToTree(ordering,get_size);

          _size += get_size.size();
        }

        template<typename GFS, typename Ordering>
        void accumulate_size(const GFS& gfs, const Ordering& ordering, CouplingGFSTag tag)
        {
          Dune::PDELab::get_size_for_entity<EntityIndex> get_size(_entity_index);
          TypeTree::applyToTree(ordering,get_size);

          _size += get_size.size();
        }

        template<typename GFS, typename Ordering>
        void accumulate_size(const GFS& gfs, const Ordering& ordering, SubDomainGFSTag tag)
        {
          typedef typename GFS::Traits::GridViewType::template Codim<Entity::codimension>::EntityPointer SDEP;
          const SDEP ep = gfs.gridView().grid().subDomainEntityPointer(_entity);
          if (!gfs.gridView().indexSet().contains(*ep))
            return;

          EntityIndex ei;

          Traits::DOFIndexAccessor::GeometryIndex::store(
            ei,
            ep->type(),
            gfs.gridView().indexSet().index(*ep)
          );

          Dune::PDELab::get_size_for_entity<EntityIndex> get_size(ei);
          TypeTree::applyToTree(ordering,get_size);

          _size += get_size.size();
        }

        get_size_for_entity(const Entity& entity,
                            const EntityIndex& entity_index)
          : _size(0)
          , _entity(entity)
          , _entity_index(entity_index)
        {}

        std::size_t size() const
        {
          return _size;
        }

      private:

        std::size_t _size;
        const Entity& _entity;
        const EntityIndex& _entity_index;

      };


      template<typename Entity, typename Traits>
      struct indices_for_entity
        : public TypeTree::DirectChildrenVisitor
        , public TypeTree::StaticTraversal
      {

        typedef typename Traits::DOFIndex::EntityIndex EntityIndex;
        typedef typename Traits::ContainerIndex ContainerIndex;

        typedef std::size_t size_type;

        template<typename GFS, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(const GFS& gfs, const Child& child, TreePath tp, ChildIndex childIndex)
        {
          collect_indices(child,
                          gfs.ordering().template child<ChildIndex::value>(),
                          typename gfs_flavor_tag<Child>::type());
        }

        template<typename GFS, typename Ordering>
        void collect_indices(const GFS& gfs, const Ordering& ordering, MultiDomainGFSTag tag)
        {
          Dune::PDELab::indices_for_entity<
            EntityIndex,
            ContainerIndex,
            TypeTree::TreeInfo<Ordering>::depth
            > extract_indices(_entity_index,_indices,_it);
          TypeTree::applyToTree(ordering,extract_indices);
          _end = extract_indices.end();
        }

        template<typename GFS, typename Ordering>
        void collect_indices(const GFS& gfs, const Ordering& ordering, CouplingGFSTag tag)
        {
          Dune::PDELab::indices_for_entity<
            EntityIndex,
            ContainerIndex,
            TypeTree::TreeInfo<Ordering>::depth
            > extract_indices(_entity_index,_indices,_it);
          TypeTree::applyToTree(ordering,extract_indices);
          _end = extract_indices.end();
        }

        template<typename GFS, typename Ordering>
        void collect_indices(const GFS& gfs, const Ordering& ordering, SubDomainGFSTag tag)
        {
          typedef typename GFS::Traits::GridViewType::template Codim<Entity::codimension>::EntityPointer SDEP;
          const SDEP ep = gfs.gridView().grid().subDomainEntityPointer(_entity);
          if (!gfs.gridView().indexSet().contains(*ep))
            return;

          EntityIndex ei;

          Traits::DOFIndexAccessor::GeometryIndex::store(
            ei,
            ep->type(),
            gfs.gridView().indexSet().index(*ep)
          );

          Dune::PDELab::indices_for_entity<
            EntityIndex,
            ContainerIndex,
            TypeTree::TreeInfo<Ordering>::depth
            > extract_indices(ei,_indices,_it);
          TypeTree::applyToTree(ordering,extract_indices);
          _end = extract_indices.end();
        }

        template<typename Ordering, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(const Ordering& ordering, const Child& child, TreePath tp, ChildIndex childIndex)
        {
          ordering.containerIndices(_entity_index,
                                    childIndex,
                                    _it,
                                    _end);
          _it = _end;
        }


        indices_for_entity(const Entity& _entity,
                           const EntityIndex& entity_index,
                           std::vector<ContainerIndex>& indices)
          : _entity_index(entity_index)
          , _indices(indices)
          , _it(indices.begin())
          , _end(indices.begin())
        {}

      private:

        typedef typename std::vector<ContainerIndex>::iterator Iterator;

        const Entity& _entity;
        const EntityIndex& _entity_index;
        std::vector<ContainerIndex>& _indices;
        Iterator _it;
        Iterator _end;

      };

    } // anonymous namespace


    template<typename GFS>
    class DataHandleProvider
    {

    public:

      typedef std::size_t size_type;

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        return gfs().ordering().contains(codim);
      }

      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        return gfs().ordering().fixedSize(codim);
      }

      /*! how many objects of type DataType have to be sent for a given entity

        Note: Only the sender side needs to know this size.
      */
      template<class Entity>
      size_type dataHandleSize (const Entity& e) const
      {
        typedef typename GFS::Ordering Ordering;

        typedef typename Ordering::Traits::DOFIndex::EntityIndex EntityIndex;
        EntityIndex ei;

        Ordering::Traits::DOFIndexAccessor::GeometryIndex::store(
          ei,
          e.type(),
          gfs().gridView().indexSet().index(e)
        );

        get_size_for_entity<Entity,typename Ordering::Traits> get_size(e,ei);
        TypeTree::applyToTree(gfs(),get_size);

        return get_size.size();
      }

      //! return vector of global indices associated with the given entity
      template<typename Entity, typename ContainerIndex>
      size_type dataHandleContainerIndices (const Entity& e,
                                            std::vector<ContainerIndex>& indices) const
      {
        typedef typename GFS::Ordering Ordering;

        dune_static_assert((is_same<ContainerIndex,typename Ordering::Traits::ContainerIndex>::value),
                           "dataHandleContainerIndices() called with invalid ContainerIndex type.");

        // Clear index state
        for (typename std::vector<ContainerIndex>::iterator it = indices.begin(),
               endit = indices.end();
             it != endit;
             ++it)
          it->clear();

        typedef typename Ordering::Traits::DOFIndex::EntityIndex EntityIndex;
        EntityIndex ei;

        Ordering::Traits::DOFIndexAccessor::GeometryIndex::store(
          ei,
          e.type(),
          gfs().gridView().indexSet().index(e)
        );

        get_size_for_entity<Entity,typename Ordering::Traits> get_size(e,ei);
        TypeTree::applyToTree(gfs(),get_size);

        indices.resize(get_size.size());

        indices_for_entity<
          Entity,
          typename Ordering::Traits
          > extract_indices(e,ei,indices);
        TypeTree::applyToTree(gfs(),extract_indices);

        return get_size.size();
      }

    private:

      const GFS& gfs() const
      {
        return static_cast<const GFS&>(*this);
      }

    };

    } // namespace MultiDomain
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_DATAHANDLEPROVIDER
