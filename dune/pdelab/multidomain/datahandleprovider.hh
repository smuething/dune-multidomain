// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_MULTIDOMAIN_DATAHANDLEPROVIDER_HH
#define DUNE_PDELAB_MULTIDOMAIN_DATAHANDLEPROVIDER_HH

#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>
#include <dune/pdelab/multidomain/multidomainlocalfunctionspace.hh>

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


      template<typename Entity, typename Traits, bool map_dof_indices>
      struct indices_for_entity
        : public TypeTree::DirectChildrenVisitor
        , public TypeTree::StaticTraversal
      {

        typedef typename Traits::DOFIndex DOFIndex;
        typedef typename DOFIndex::EntityIndex EntityIndex;
        typedef typename Traits::ContainerIndex ContainerIndex;

        typedef std::size_t size_type;

        typedef typename std::vector<ContainerIndex>::iterator CIIterator;
        typedef typename std::conditional<
          map_dof_indices,
          typename std::vector<DOFIndex>::iterator,
          DummyDOFIndexIterator
          >::type DIIterator;

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
            DOFIndex,
            ContainerIndex,
            TypeTree::TreeInfo<Ordering>::depth,
            map_dof_indices
            > extract_indices(_entity_index,_ci_it,_di_it);
          TypeTree::applyToTree(ordering,extract_indices);
          _ci_end = extract_indices.ci_end();
          _di_end = extract_indices.di_end();
        }

        template<typename GFS, typename Ordering>
        void collect_indices(const GFS& gfs, const Ordering& ordering, CouplingGFSTag tag)
        {
          Dune::PDELab::indices_for_entity<
            DOFIndex,
            ContainerIndex,
            TypeTree::TreeInfo<Ordering>::depth,
            map_dof_indices
            > extract_indices(_entity_index,_ci_it,_di_it);
          TypeTree::applyToTree(ordering,extract_indices);
          _ci_end = extract_indices.ci_end();
          _di_end = extract_indices.di_end();
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
            DOFIndex,
            ContainerIndex,
            TypeTree::TreeInfo<Ordering>::depth,
            map_dof_indices
            > extract_indices(ei,_ci_it,_di_it);
          TypeTree::applyToTree(ordering,extract_indices);
          _ci_end = extract_indices.ci_end();
          _di_end = extract_indices.di_end();
        }

        template<typename GFS, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(const GFS& gfs, const Child& child, TreePath tp, ChildIndex childIndex)
        {
          gfs.ordering().extract_entity_indices(_entity_index,
                                                childIndex,
                                                _ci_it,
                                                _ci_end);

          if (GFS::Ordering::consume_tree_index)
            for (DIIterator it = _di_it;
                 it != _di_end;
                 ++it)
              it->treeIndex().push_back(childIndex);

          _ci_it = _ci_end;
          _di_it = _di_end;

        }


        indices_for_entity(const Entity& entity,
                           const EntityIndex& entity_index,
                           CIIterator ci_begin,
                           DIIterator di_begin = DIIterator())
          : _entity(entity)
          , _entity_index(entity_index)
          , _ci_it(ci_begin)
          , _ci_end(ci_begin)
          , _di_it(di_begin)
          , _di_end(di_begin)
        {}

      private:

        const Entity& _entity;
        const EntityIndex& _entity_index;
        CIIterator _ci_it;
        CIIterator _ci_end;
        DIIterator _di_it;
        DIIterator _di_end;

      };

    } // anonymous namespace


    template<typename GFS>
    class DataHandleProvider
      : public Dune::PDELab::DataHandleProvider<GFS>
    {

    public:

      typedef typename Dune::PDELab::DataHandleProvider<GFS>::size_type size_type;

    private:

      using Dune::PDELab::DataHandleProvider<GFS>::gfs;

    public:

      //------------------------------
      // generic data handle interface
      //------------------------------

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
      template<typename Entity, typename ContainerIndex, typename DOFIndex, bool map_dof_indices>
      size_type dataHandleIndices (const Entity& e,
                                   std::vector<ContainerIndex>& container_indices,
                                   std::vector<DOFIndex>& dof_indices,
                                   std::integral_constant<bool,map_dof_indices> map_dof_indices_value) const
      {
        typedef typename GFS::Ordering Ordering;

        dune_static_assert((is_same<ContainerIndex,typename Ordering::Traits::ContainerIndex>::value),
                           "dataHandleContainerIndices() called with invalid ContainerIndex type.");

        // Clear index state
        for (typename std::vector<ContainerIndex>::iterator it = container_indices.begin(),
               endit = container_indices.end();
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

        container_indices.resize(get_size.size());

        this->setup_dof_indices(dof_indices,get_size.size(),ei,map_dof_indices_value);

        indices_for_entity<
          Entity,
          typename Ordering::Traits,
          map_dof_indices
          > extract_indices(e,ei,container_indices.begin(),this->dof_indices_begin(dof_indices,map_dof_indices_value));
        TypeTree::applyToTree(gfs(),extract_indices);

        return get_size.size();
      }

    };

    } // namespace MultiDomain
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTIDOMAIN_DATAHANDLEPROVIDER
