#ifndef DUNE_MULTIDOMAIN_MULTIDOMAINGRIDFUNCTIONSPACE_HH
#define DUNE_MULTIDOMAIN_MULTIDOMAINGRIDFUNCTIONSPACE_HH

#include <tuple>
#include <type_traits>
#include <dune/pdelab/common/multitypetree.hh>
#include <dune/pdelab/common/countingptr.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/grid/multidomaingrid.hh>
#include <utility>

namespace Dune {

namespace PDELab {

namespace MultiDomain {


/**
 * TMP for deriving storage types in VariadicCompositeNode
 */
template<typename... OArgs>
struct replace;

template<typename OHead, typename... OTail>
struct replace<OHead,OTail...> {
  template<template<typename T> class Transform, typename... SArgs>
  struct with {
    typedef typename replace<OTail...>::template with<Transform,SArgs...,typename Transform<OHead>::type>::type type;
  };
};

/**
 * End of recursion - export the transformed argument list as a tuple
 */
template<>
struct replace<> {
  template<template<typename T> class Transform, typename... SArgs>
  struct with {
    typedef typename Transform<void>::template container<SArgs...>::type type;
  };
};

/**
 * wrapper to simplify calling the TMP
 */
template<template<typename T> class Transform, typename... Args>
struct transform {
  typedef typename replace<Args...>::template with<Transform>::type type;
};


/**
 * VariadicCompositeNode
 *
 */
template<typename P, typename... Children>
class VariadicCompositeNode
{

  template<typename T>
  struct forwarder {
    typedef typename P::template Storage<T>::Type type;

    template<typename... Args>
    struct container {
      typedef std::tuple<Args...> type;
    };
  };

  typedef std::tuple<Children...> OT;
  typedef typename transform<forwarder,Children...>::type ST;

public:

  enum { isLeaf = false };
  enum { isPower = false /**< */ };
  enum { isComposite = true /**< */ };
  enum { CHILDREN = sizeof...(Children) };

  VariadicCompositeNode (Children&&... c_) : c(P::convert(c_)...) {}

  template<int i>
  struct Child
  {
    typedef typename std::tuple_element<i,OT>::type Type;
  };

  template<int i>
  typename Child<i>::Type& getChild ()
  {
    return P::get(std::get<i>(c));
  }

  template<int i>
  const typename Child<i>::Type& getChild () const
  {
    return P::get(std::get<i>(c));
  }

private:
  ST c;
};



template<typename T, int n, int i>
struct MultiDomainGridFunctionSpaceVisitChildMetaProgram // visit child of inner node
{

  typedef MultiDomainGridFunctionSpaceVisitChildMetaProgram<T,n,i+1> NextChild;

  template<typename Int>
  static void setup (T& t, Int childGlobalSize[], Int childLocalSize[])
  {
    childGlobalSize[i] = t.template getChild<i>().globalSize();
    childLocalSize[i] = t.template getChild<i>().maxLocalSize();
    NextChild::setup(t,childGlobalSize,childLocalSize);
  }
  static void update (T& t)
  {
    t.template getChild<i>().update();
    NextChild::update(t);
  }
  static bool dataHandleContains (const T& t, int dim, int codim)
  {
    return t.template getChild<i>().dataHandleContains(dim,codim) ||
      NextChild::dataHandleContains(t,dim,codim);
  }
  static bool dataHandleFixedSize (const T& t, int dim, int codim)
  {
    return t.template getChild<i>().dataHandleFixedSize(dim,codim) &&
      NextChild::dataHandleFixedSize(t,dim,codim);
  }
  template<class EntityType>
  static size_t dataHandleSize (const T& t, const EntityType& e)
  {
    return t.template getChild<i>().dataHandleSize(e) +
      NextChild::dataHandleSize(t,e);
  }
  template<class EntityType, class C>
  static void dataHandleGlobalIndices (const T& t, const EntityType& e, C& global, size_t ng, C& childglobal)
  {
    size_t nc=t.template getChild<i>().dataHandleSize(e);
    childglobal.resize(nc);
    t.template getChild<i>().dataHandleGlobalIndices(e,childglobal);
    for (size_t j=0; j<childglobal.size(); j++)
      global[ng+j] = t.template subMap<i>(childglobal[j]);
    NextChild::dataHandleGlobalIndices(t,e,global,ng+nc,childglobal);
  }
  static void verifyChild(const T& t)
  {
    typedef typename T::Traits::GridViewType::Grid MultiDomainGrid;
    typedef typename MultiDomainGrid::SubDomainGrid SubDomainGrid;
    typedef typename T::template Child<i>::Type::Traits::GridViewType::Grid ChildGrid;
    dune_static_assert((std::is_same<MultiDomainGrid,ChildGrid>::value || std::is_same<SubDomainGrid,ChildGrid>::value)
                       ,"MultiDomainGridFunctionSpace only works with a MultiDomainGrid and its associated SubDomainGrids.");
    doVerify(t.grid(),t.template getChild<i>().gridview().grid());
    NextChild::verifyChild(t);
  }

  static void doVerify(const typename T::Traits::GridViewType::Grid& g, const typename T::Traits::GridViewType::Grid& cg)
  {
    assert(&g == &cg);
  }

  static void doVerify(const typename T::Traits::GridViewType::Grid& g, const typename T::Traits::GridViewType::Grid::SubDomainGrid& cg)
  {
    assert(&g == &cg.multiDomainGrid());
  }

  template<typename ST>
  static void doVerify(const typename T::Traits::GridViewType::Grid& g, const ST& cg)
  {
    // this is only here to keep the compiler from complaining about a missing function
    // if ST is not a Multi-/SubDomainGrid
  }

};

template<typename T, int n>
struct MultiDomainGridFunctionSpaceVisitChildMetaProgram<T,n,n> // end of child recursion
{
  template<typename Int>
  static void setup (T& t, Int childGlobalSize[], Int childLocalSize[])
  {
  }
  static void update (T& t)
  {
  }
  static bool dataHandleContains (const T& t, int dim, int codim)
  {
    return false;
  }
  static bool dataHandleFixedSize (const T& t, int dim, int codim)
  {
    return true;
  }
  template<class EntityType>
  static size_t dataHandleSize (const T& t, const EntityType& e)
  {
    return 0;
  }
  template<class EntityType, class C>
  static void dataHandleGlobalIndices (const T& t, const EntityType& e, C& global, size_t n_, C& childglobal)
  {
  }
  static void verifyChild(const T& t)
  {
  }
};

template<typename G, typename... Children>
class MultiDomainGridFunctionSpace : public Countable, public VariadicCompositeNode<CopyStoragePolicy,Children...>
{

  typedef VariadicCompositeNode<CopyStoragePolicy,Children...> BaseT;
  typedef MultiDomainGridFunctionSpaceVisitChildMetaProgram<MultiDomainGridFunctionSpace,sizeof...(Children),0> VisitChildTMP;

public:
  //! export traits class
  typedef PowerCompositeGridFunctionSpaceTraits<typename BaseT::template Child<0>::Type::Traits::GridViewType,
                                                typename BaseT::template Child<0>::Type::Traits::BackendType,
                                                CopyStoragePolicy,
                                                sizeof...(Children)>
  Traits;

  //! extract type of container storing Es
  template<typename E>
  struct VectorContainer
  {
    //! \brief define Type as the Type of a container of E's
    typedef typename Traits::BackendType::template VectorContainer<MultiDomainGridFunctionSpace,E> Type;
  private:
    VectorContainer () {}
  };

  //! extract type for storing constraints
  template<typename E>
  struct ConstraintsContainer
  {
    //! \brief define Type as the Type of a container of E's
    typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;
  private:
    ConstraintsContainer () {}
  };

  // define local function space parametrized by self
  //typedef Dune::PDELab::CompositeLocalFunctionSpace<CompositeGridFunctionSpace> LocalFunctionSpace;


  MultiDomainGridFunctionSpace (G& g, Children&... children) : BaseT(children...), _g(g)
  {
    dune_static_assert(Dune::mdgrid::GridType<G>::v == Dune::mdgrid::multiDomainGrid,
                       "MultiDomainGridFunctionSpace only works on a MultiDomainGrid");
    VisitChildTMP::verifyChild(*this);
  }

  // get grid view
  const typename Traits::GridViewType& gridview () const
  {
    return this->template getChild<0>().gridview();
  }

  //! get dimension of root finite element space
  typename Traits::SizeType globalSize () const
  {
    return offset[BaseT::CHILDREN];
  }

  //! get dimension of this finite element space
  typename Traits::SizeType size () const
  {
    return offset[BaseT::CHILDREN];
  }

  // get max dimension of shape function space
  typename Traits::SizeType maxLocalSize () const
  {
    // this is bullshit !
    return maxlocalsize;
  }

  //! map index from our index set [0,size()-1] to root index set
  typename Traits::SizeType upMap (typename Traits::SizeType i) const
  {
    return i;
  }

  //! map index from child i's index set into our index set
  template<int i>
  typename Traits::SizeType subMap (typename Traits::SizeType j) const
  {
    return offset[i]+j;
  }


  //------------------------------
  // generic data handle interface
  //------------------------------

  //! returns true if data for this codim should be communicated
  bool dataHandleContains (int dim, int codim) const
  {
    return VisitChildTMP::dataHandleContains(*this,dim,codim);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool dataHandleFixedSize (int dim, int codim) const
  {
    return VisitChildTMP::dataHandleFixedSize(*this,dim,codim);
  }

  /*! how many objects of type DataType have to be sent for a given entity

    Note: Only the sender side needs to know this size.
*/
  template<class EntityType>
  size_t dataHandleSize (const EntityType& e) const
  {
    return VisitChildTMP::dataHandleSize(*this,e);
  }

  //! return vector of global indices associated with the given entity
  template<class EntityType>
  void dataHandleGlobalIndices (const EntityType& e,
                                std::vector<typename Traits::SizeType>& global) const
  {
    size_t n=dataHandleSize(e);
    global.resize(n);
    VisitChildTMP::dataHandleGlobalIndices(*this,e,global,0,childglobal);
  }

  //------------------------------



  // recalculate sizes
  void update ()
  {
    VisitChildTMP::update(*this);
    setup();
  }

  G& grid() {
    return _g;
  }

  const G& grid() const {
    return _g;
  }

private:

  const G& _g;

  void setup ()
  {
    Dune::dinfo << "multi domain grid function space:"
                << std::endl;

    VisitChildTMP::setup(*this,childGlobalSize,childLocalSize);

    Dune::dinfo << "( ";
    offset[0] = 0;
    maxlocalsize = 0;
    for (int i=0; i<BaseT::CHILDREN; i++)
      {
        Dune::dinfo << childGlobalSize[i] << " ";
        offset[i+1] = offset[i]+childGlobalSize[i];
        maxlocalsize += childLocalSize[i];
      }
    Dune::dinfo << ") total size = " << offset[BaseT::CHILDREN]
                << " max local size = " << maxlocalsize
                << std::endl;
    childglobal.resize(maxlocalsize);
  }

  typename Traits::SizeType childGlobalSize[BaseT::CHILDREN];
  typename Traits::SizeType childLocalSize[BaseT::CHILDREN];
  typename Traits::SizeType offset[BaseT::CHILDREN+1];
  typename Traits::SizeType maxlocalsize;
  mutable std::vector<typename Traits::SizeType> childglobal;


};


} // namespace MultiDomain

} // namespace PDELab

} // namespace Dune

#endif // DUNE_MULTIDOMAIN_MULTIDOMAINGRIDFUNCTIONSPACE_HH
