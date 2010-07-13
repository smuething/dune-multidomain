#ifndef DUNE_MULTIDOMAIN_OPERATORAPPLIER_HH
#define DUNE_MULTIDOMAIN_OPERATORAPPLIER_HH

namespace Dune {

namespace PDELab {

namespace MultiDomain {

template<typename MDGOS, typename Condition, template<bool,bool> class BooleanOp, bool start_value, std::size_t i, std::size_t n>
struct child_condition
{
  static const bool value = BooleanOp<Condition::template test<typename MDGOS::template Child<i>::Type>::value,
                                      child_condition<MDGOS,Condition,BooleanOp,start_value,i+1,n>::value
                                     >::value;
};

template<typename MDGOS, typename Condition, template<bool,bool> class BooleanOp, bool start_value, std::size_t n>
struct child_condition<MDGOS,Condition,BooleanOp,start_value,n,n>
{
  static const bool value = start_value;
};

template<bool a, bool b>
struct and_
{
  static const bool value = a && b;
};

template<bool a, bool b>
struct or_
{
  static const bool value = a || b;
};

template<typename MDGOS, typename Condition>
struct all_childs : public child_condition<MDGOS,Condition,and_,true,0,MDGOS::CHILDREN> {};

template<typename MDGOS, typename Condition>
struct any_child : public child_condition<MDGOS,Condition,or_,false,0,MDGOS::CHILDREN> {};

template<typename Applier, typename Operator, std::size_t i, std::size_t n>
struct apply_operator_helper
{
  static void apply(Applier& applier, Operator& op)
  {
    op(applier,applier.gos().template getChild<i>());
    apply_operator_helper<Applier,Operator,i+1,n>::apply(applier,op);
  }
};

// end of recursion
template<typename Applier, typename Operator, std::size_t n>
struct apply_operator_helper<Applier, Operator, n,n>
{
  static void apply(Applier& applier, Operator& op)
  {
  }
};


template<typename Applier, typename Operator, std::size_t i, bool do_apply>
struct conditional_apply;

template<typename Applier, typename Operator, std::size_t i>
struct conditional_apply<Applier,Operator,i,true>
{
  static void apply(Applier& applier, Operator& op)
  {
    op(applier,applier.gos().template getChild<i>());
  }
};

template<typename Applier, typename Operator, std::size_t i>
struct conditional_apply<Applier,Operator,i,false>
{
  static void apply(Applier& applier, Operator& op)
  {
  }
};

template<typename Applier, typename Condition, typename Operator, std::size_t i, std::size_t n>
struct conditional_apply_operator_helper
{
  static void apply(Applier& applier, Operator& op)
  {
    conditional_apply<Applier,Operator,i,Condition::template test<typename Applier::MDGOS::template Child<i>::Type>::value>::apply(applier,op);
    conditional_apply_operator_helper<Applier,Condition,Operator,i+1,n>::apply(applier,op);
  }
};

// end of recursion
template<typename Applier, typename Condition, typename Operator, std::size_t n>
struct conditional_apply_operator_helper<Applier, Condition, Operator, n,n>
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

  typedef typename MDGOS::Traits::TrialGridFunctionSpace::LocalFunctionSpace LFSU;

  NeighborTrialFunctionSpace() :
    _plfsun(NULL)
  {}

  void setlfsun(const LFSU& lfsun)
  {
    _plfsun = &lfsun;
  }

  const LFSU& lfsun() const
  {
    assert(_plfsun != NULL);
    return *_plfsun;
  }

private:
  const LFSU* _plfsun;
};

template<typename MDGOS>
class NeighborTestFunctionSpace
{

public:

  typedef typename MDGOS::Traits::TestGridFunctionSpace::LocalFunctionSpace LFSV;

  NeighborTestFunctionSpace() :
    _plfsvn(NULL)
  {}

  void setlfsvn(const LFSV& lfsvn)
  {
    _plfsvn = &lfsvn;
  }

  const LFSV& lfsvn() const
  {
    assert(_plfsvn != NULL);
    return *_plfsvn;
  }

private:
  const LFSV* _plfsvn;
};

template<typename MDGOS>
class NeighborFunctionSpaces :
    public NeighborTrialFunctionSpace<MDGOS>,
    public NeighborTestFunctionSpace<MDGOS>
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

  template<typename Operator>
  void operator()(Operator&& op)
  {
    apply_operator_helper<operator_applier,Operator,0,MDGOS::CHILDREN>::apply(*this,op);
  }

  template<typename Condition, typename Operator>
  void conditional(Operator&& op)
  {
    conditional_apply_operator_helper<operator_applier,Condition,Operator,0,MDGOS::CHILDREN>::apply(*this,op);
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


struct SpatialOperator
{

  template<typename SubProblem>
  struct ExtractType
  {
    typedef typename SubProblem::Traits::LocalOperator Type;
  };

  template<typename SubProblem>
  static typename SubProblem::Traits::LocalOperator& extract(SubProblem& subProblem) {
    return subProblem.localOperator();
  }

  template<typename SubProblem>
  static const typename SubProblem::Traits::LocalOperator& extract(const SubProblem& subProblem) {
    return subProblem.localOperator();
  }
};

struct TemporalOperator
{

  template<typename SubProblem>
  struct ExtractType
  {
    typedef typename SubProblem::Traits::TemporalOperator Type;
  };

  template<typename SubProblem>
  static typename SubProblem::Traits::TemporalOperator& extract(SubProblem& subProblem) {
    return subProblem.temporalOperator();
  }

  template<typename SubProblem>
  static const typename SubProblem::Traits::TemporalOperator& extract(const SubProblem& subProblem) {
    return subProblem.temporalOperator();
  }
};


template<typename Operator = SpatialOperator>
struct do_pattern_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doPatternSkeleton;
  };
};

template<typename Operator = SpatialOperator>
struct do_pattern_volume
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doPatternVolume;
  };
};

template<typename Operator = SpatialOperator>
struct do_alpha_volume
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doAlphaVolume;
  };
};

template<typename Operator = SpatialOperator>
struct do_alpha_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doAlphaSkeleton;
  };
};

template<typename Operator = SpatialOperator>
struct do_alpha_boundary
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doAlphaBoundary;
  };
};

template<typename Operator = SpatialOperator>
struct do_alpha_skeleton_or_boundary
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doAlphaSkeleton || Operator::template ExtractType<T>::Type::doAlphaBoundary;
  };
};

template<typename Operator = SpatialOperator>
struct do_lambda_volume
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doLambdaVolume;
  };
};

template<typename Operator = SpatialOperator>
struct do_lambda_boundary
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doLambdaBoundary;
  };
};

template<typename Operator = SpatialOperator>
struct do_alpha_volume_post_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doAlphaVolumePostSkeleton;
  };
};

template<typename Operator = SpatialOperator>
struct do_lambda_volume_post_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doLambdaVolumePostSkeleton;
  };
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_OPERATORAPPLIER_HH
