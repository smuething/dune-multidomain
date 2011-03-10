#ifndef DUNE_MULTIDOMAIN_OPERATORFLAGTESTS_HH
#define DUNE_MULTIDOMAIN_OPERATORFLAGTESTS_HH

namespace Dune {

namespace PDELab {

namespace MultiDomain {

namespace {

  template<typename Predicate>
  struct predicate_wrapper
  {

    typedef bool result_type;

    template<typename Node, typename TreePath>
    struct doVisit
    {
      static const bool value = Node::isLeaf;
    };

    template<typename Node, typename TreePath>
    struct visit
    {
      static const bool result = Predicate::template test<Node>::value;
    };

  };

} // anonymous namespace

template<typename Tree, typename Predicate>
struct any_child
{
  static const bool value = Dune::PDELab::TypeTree::AccumulateValue<Tree,
                                                                    predicate_wrapper<Predicate>,
                                                                    Dune::PDELab::TypeTree::or_<bool>,
                                                                    false>::result;
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

struct CouplingOperator
{

  template<typename Coupling>
  struct ExtractType
  {
    typedef typename Coupling::Traits::CouplingOperator Type;
  };

  template<typename Coupling>
  static typename Coupling::Traits::CouplingOperator& extract(Coupling& coupling) {
    return coupling.couplingOperator();
  }

  template<typename Coupling>
  static const typename Coupling::Traits::CouplingOperator& extract(const Coupling& coupling) {
    return coupling.couplingOperator();
  }
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
struct do_pattern_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doPatternSkeleton;
  };
};

template<typename Operator = SpatialOperator>
struct do_pattern_boundary
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doPatternBoundary;
  };
};

template<typename Operator = SpatialOperator>
struct do_pattern_skeleton_or_boundary
{
  template<typename T>
  struct test {
    static const bool value =
      Operator::template ExtractType<T>::Type::doPatternSkeleton ||
      Operator::template ExtractType<T>::Type::doPatternBoundary;
  };
};

template<typename Operator = CouplingOperator>
struct do_pattern_coupling
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doPatternCoupling;
  };
};

template<typename Operator = CouplingOperator>
struct do_pattern_enriched_coupling
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doPatternEnrichedCoupling;
  };
};

template<typename Operator = SpatialOperator>
struct do_pattern_volume_post_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doPatternVolumePostSkeleton;
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
    static const bool value =
      Operator::template ExtractType<T>::Type::doAlphaSkeleton ||
      Operator::template ExtractType<T>::Type::doAlphaBoundary;
  };
};

template<typename Operator = CouplingOperator>
struct do_alpha_coupling
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doAlphaCoupling;
  };
};

template<typename Operator = CouplingOperator>
struct do_alpha_enriched_coupling
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doAlphaEnrichedCoupling;
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
struct do_lambda_volume
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doLambdaVolume;
  };
};

template<typename Operator = SpatialOperator>
struct do_lambda_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doLambdaSkeleton;
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
struct do_lambda_skeleton_or_boundary
{
  template<typename T>
  struct test {
    static const bool value =
      Operator::template ExtractType<T>::Type::doLambdaSkeleton ||
      Operator::template ExtractType<T>::Type::doLambdaBoundary;
  };
};

template<typename Operator = CouplingOperator>
struct do_lambda_coupling
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doLambdaCoupling;
  };
};

template<typename Operator = CouplingOperator>
struct do_lambda_enriched_coupling
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doLambdaEnrichedCoupling;
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

template<typename Operator = SpatialOperator>
struct do_skeleton_two_sided
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doSkeletonTwoSided;
  };
};

// ********************************************************************************
// tests for constraints operators
// ********************************************************************************

struct do_constraints_volume
{
  template<typename T>
  struct test {
    static const bool value = T::doVolume;
  };
};

struct do_constraints_skeleton
{
  template<typename T>
  struct test {
    static const bool value = T::doSkeleton;
  };
};

struct do_constraints_boundary
{
  template<typename T>
  struct test {
    static const bool value = T::doBoundary;
  };
};

struct do_constraints_processor
{
  template<typename T>
  struct test {
    static const bool value = T::doProcessor;
  };
};



} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_OPERATORFLAGTESTS_HH
