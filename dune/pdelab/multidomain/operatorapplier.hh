#ifndef DUNE_MULTIDOMAIN_OPERATORAPPLIER_HH
#define DUNE_MULTIDOMAIN_OPERATORAPPLIER_HH

namespace Dune {

namespace PDELab {

namespace MultiDomain {


template<typename Applier, typename Operator, typename Extractor, std::size_t i, std::size_t n>
struct apply_operator_helper
{
  static void apply(Applier& applier, Operator& op)
  {
    op(applier,Extractor::template get<i>(applier.gos()));
    apply_operator_helper<Applier,Operator,Extractor,i+1,n>::apply(applier,op);
  }
};

// end of recursion
template<typename Applier, typename Operator, typename Extractor, std::size_t n>
struct apply_operator_helper<Applier, Operator, Extractor, n, n>
{
  static void apply(Applier& applier, Operator& op)
  {
  }
};


template<typename Applier, typename Operator, typename Extractor, std::size_t i, bool do_apply>
struct conditional_apply;

template<typename Applier, typename Operator, typename Extractor, std::size_t i>
struct conditional_apply<Applier,Operator,Extractor,i,true>
{
  static void apply(Applier& applier, Operator& op)
  {
    op(applier,Extractor::template get<i>(applier.gos()));
  }
};

template<typename Applier, typename Operator, typename Extractor, std::size_t i>
struct conditional_apply<Applier,Operator,Extractor,i,false>
{
  static void apply(Applier& applier, Operator& op)
  {
  }
};

template<typename Applier, typename Condition, typename Operator, typename Extractor, std::size_t i, std::size_t n>
struct conditional_apply_operator_helper
{
  static void apply(Applier& applier, Operator& op)
  {
    conditional_apply<Applier,Operator,Extractor,i,Condition::template test<typename Extractor::template GetType<i>::Type>::value>::apply(applier,op);
    conditional_apply_operator_helper<Applier,Condition,Operator,Extractor,i+1,n>::apply(applier,op);
  }
};

// end of recursion
template<typename Applier, typename Condition, typename Operator, typename Extractor, std::size_t n>
struct conditional_apply_operator_helper<Applier, Condition, Operator, Extractor, n, n>
{
  static void apply(Applier& applier, Operator& op)
  {
  }
};





// data namespace ************************************************


namespace data {

template<typename MDGOS>
class LocalTrialFunctionSpace
{

public:

  typedef typename MDGOS::Traits::TrialGridFunctionSpace::LocalFunctionSpace LFSU;

  LocalTrialFunctionSpace() :
    _plfsu(NULL)
  {}

  void setlfsu(const LFSU& lfsu)
  {
    _plfsu = &lfsu;
  }

  const LFSU& lfsu() const
  {
    assert(_plfsu != NULL);
    return *_plfsu;
  }

private:
  const LFSU* _plfsu;
};

template<typename MDGOS>
class LocalTestFunctionSpace
{

public:

  typedef typename MDGOS::Traits::TestGridFunctionSpace::LocalFunctionSpace LFSV;

  LocalTestFunctionSpace() :
    _plfsv(NULL)
  {}

  void setlfsv(const LFSV& lfsv)
  {
    _plfsv = &lfsv;
  }

  const LFSV& lfsv() const
  {
    assert(_plfsv != NULL);
    return *_plfsv;
  }

private:
  const LFSV* _plfsv;
};

template<typename MDGOS>
class LocalFunctionSpaces :
    public LocalTrialFunctionSpace<MDGOS>,
    public LocalTestFunctionSpace<MDGOS>
{
};


template<typename MDGOS>
class NeighborTrialFunctionSpace
{

public:

  typedef typename MDGOS::Traits::TrialGridFunctionSpace::LocalFunctionSpace LFSU_;

  NeighborTrialFunctionSpace() :
    _plfsun(NULL)
  {}

  void setlfsun(const LFSU_& lfsun)
  {
    _plfsun = &lfsun;
  }

  const LFSU_& lfsun() const
  {
    assert(_plfsun != NULL);
    return *_plfsun;
  }

private:
  const LFSU_* _plfsun;
};

template<typename MDGOS>
class NeighborTestFunctionSpace
{

public:

  typedef typename MDGOS::Traits::TestGridFunctionSpace::LocalFunctionSpace LFSV_;

  NeighborTestFunctionSpace() :
    _plfsvn(NULL)
  {}

  void setlfsvn(const LFSV_& lfsvn)
  {
    _plfsvn = &lfsvn;
  }

  const LFSV_& lfsvn() const
  {
    assert(_plfsvn != NULL);
    return *_plfsvn;
  }

private:
  const LFSV_* _plfsvn;
};

template<typename MDGOS>
class NeighborFunctionSpaces :
    public NeighborTrialFunctionSpace<MDGOS>,
    public NeighborTestFunctionSpace<MDGOS>
{
};


template<typename MDGOS>
class CouplingTrialFunctionSpace
{

public:

  typedef typename MDGOS::Traits::TrialGridFunctionSpace::CouplingLocalFunctionSpace CouplingLFSU;

  CouplingTrialFunctionSpace() :
    _pcouplinglfsu(NULL)
  {}

  void setcouplinglfsu(const CouplingLFSU& couplinglfsu)
  {
    _pcouplinglfsu = &couplinglfsu;
  }

  const CouplingLFSU& couplinglfsu() const
  {
    assert(_pcouplinglfsu != NULL);
    return *_pcouplinglfsu;
  }

private:
  const CouplingLFSU* _pcouplinglfsu;
};


template<typename MDGOS>
class CouplingTestFunctionSpace
{

public:

  typedef typename MDGOS::Traits::TestGridFunctionSpace::CouplingLocalFunctionSpace CouplingLFSV;

  CouplingTestFunctionSpace() :
    _pcouplinglfsv(NULL)
  {}

  void setcouplinglfsv(const CouplingLFSV& couplinglfsv)
  {
    _pcouplinglfsv = &couplinglfsv;
  }

  const CouplingLFSV& couplinglfsv() const
  {
    assert(_pcouplinglfsv != NULL);
    return *_pcouplinglfsv;
  }

private:
  const CouplingLFSV* _pcouplinglfsv;
};


template<typename MDGOS>
class CouplingFunctionSpaces :
    public CouplingTrialFunctionSpace<MDGOS>,
    public CouplingTestFunctionSpace<MDGOS>
{
};


template<typename MDGOS>
class ElementReference
{

public:

  typedef typename MDGOS::Traits::TrialGridFunctionSpace::Traits::GridType Grid;
  typedef typename Grid::template Codim<0>::Entity Element;

  ElementReference() :
    _element(NULL)
  {}

  void setElement(const Element& element) {
    _element = &element;
  }

  const Element& element() const {
    assert(_element != NULL);
    return *_element;
  }

private:
  const Element* _element;

};

template<typename MDGOS>
class IntersectionReference
{

public:

  typedef typename MDGOS::Traits::TrialGridFunctionSpace::Traits::GridType Grid;
  typedef typename Grid::LeafGridView::Intersection Intersection;

  IntersectionReference() :
    _intersection(NULL)
  {}

  void setIntersection(const Intersection& intersection, int intersectionIndex) {
    _intersection = &intersection;
    _intersectionIndex = intersectionIndex;
  }

  const Intersection& intersection() const {
    assert(_intersection != NULL);
    return *_intersection;
  }

  int intersectionIndex() const {
    assert(_intersection != NULL);
    return _intersectionIndex;
  }


private:
  const Intersection* _intersection;
  int _intersectionIndex;
};

template<typename MDGOS>
class ElementSubDomains
{

public:

  typedef typename MDGOS::Traits::TrialGridFunctionSpace::Traits::GridType Grid;
  typedef typename Grid::Traits::LeafIndexSet::SubDomainSet ElementSubDomainSet;

  const ElementSubDomainSet& elementSubDomains() const {
    return _elementSet;
  }

  void setElementSubDomains(const ElementSubDomainSet& elementSet) {
    _elementSet = elementSet;
  }

private:
  ElementSubDomainSet _elementSet;

};

template<typename MDGOS>
class NeighborSubDomains
{

public:

  typedef typename MDGOS::Traits::TrialGridFunctionSpace::Traits::GridType Grid;
  typedef typename Grid::Traits::LeafIndexSet::SubDomainSet ElementSubDomainSet;

  const ElementSubDomainSet& neighborSubDomains() const {
    return _neighborElementSet;
  }

  void setNeighborSubDomains(const ElementSubDomainSet& neighborElementSet) {
    _neighborElementSet = neighborElementSet;
  }

private:
  ElementSubDomainSet _neighborElementSet;

};

template<typename MDGOS>
class SkeletonInvocationTracker
{
public:
  SkeletonInvocationTracker() :
    _alphaSkeletonInvoked(false)
  {}

  void setAlphaSkeletonInvoked() const {
    _alphaSkeletonInvoked = true;
  }

  void clearAlphaSkeletonInvoked() {
    _alphaSkeletonInvoked = false;
  }

  bool alphaSkeletonInvoked() const {
    return _alphaSkeletonInvoked;
  }

private:
  mutable bool _alphaSkeletonInvoked;

};


template<typename MDGOS>
class EnrichedCouplingInvocationTracker
{
public:
  EnrichedCouplingInvocationTracker() :
    _alphaEnrichedCouplingInvoked(false)
  {}

  void setAlphaEnrichedCouplingInvoked() const {
    _alphaEnrichedCouplingInvoked = true;
  }

  void clearAlphaEnrichedCouplingInvoked() {
    _alphaEnrichedCouplingInvoked = false;
  }

  bool alphaEnrichedCouplingInvoked() const {
    return _alphaEnrichedCouplingInvoked;
  }

private:
  mutable bool _alphaEnrichedCouplingInvoked;

};

template<typename MDGOS>
class ElementData :
  public LocalFunctionSpaces<MDGOS>,
  public ElementReference<MDGOS>,
  public ElementSubDomains<MDGOS>
{};

template<typename MDGOS>
class NeighborData :
  public NeighborFunctionSpaces<MDGOS>,
  public NeighborSubDomains<MDGOS>
{};


}

template<typename MDGOS_, template<typename> class... Data>
class operator_applier : public Data<MDGOS_>...
{

public:


  typedef MDGOS_ MDGOS;
  typedef MDGOS_ GOS;

  operator_applier(MDGOS& mdgos) :
    _mdgos(mdgos)
  {}

  template<typename Extractor, typename Operator>
  void apply(Operator&& op)
  {
    apply_operator_helper<operator_applier,Operator,Extractor,0,Extractor::N>::apply(*this,op);
  }

  template<typename Extractor, typename Condition, typename Operator>
  void conditional(Operator&& op)
  {
    conditional_apply_operator_helper<operator_applier,Condition,Operator,Extractor,0,Extractor::N>::apply(*this,op);
  }

  MDGOS& gos() {
    return _mdgos;
  }

  const MDGOS& gos() const {
    return _mdgos;
  }

private:

  MDGOS& _mdgos;

};


} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_OPERATORAPPLIER_HH
