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
  static const bool value = TypeTree::AccumulateValue<
    Tree,
    predicate_wrapper<Predicate>,
    TypeTree::or_<bool>,
    false
    >::result;
};




struct DefaultOperator
{

  template<typename Participant>
  struct ExtractType
  {
    typedef typename Participant::Traits::LocalOperator Type;
  };

  template<typename Participant>
  static typename Participant::Traits::LocalOperator& extract(Participant& participant) {
    return participant.localOperator();
  }

  template<typename Participant>
  static const typename Participant::Traits::LocalOperator& extract(const Participant& participant) {
    return participant.localOperator();
  }
};


struct IdentityExtractor
{

  template<typename Descriptor>
  struct ExtractType
  {
    typedef Descriptor Type;
  };

  template<typename Descriptor>
  static Descriptor& extract(Descriptor& descriptor) {
    return descriptor;
  }

  template<typename Descriptor>
  static const Descriptor& extract(const Descriptor& descriptor) {
    return descriptor;
  }
};

struct LFSConstraintsExtractor
{

  template<typename LFS>
  struct ExtractType
  {
    typedef typename LFS::Traits::ConstraintsType Type;
  };

  template<typename LFS>
  static typename LFS::Traits::ConstraintsType& extract(LFS& lfs) {
    return lfs.constraints();
  }

  template<typename LFS>
  static const typename LFS::Traits::ConstraintsType& extract(const LFS& lfs) {
    return lfs.constraints();
  }
};


template<typename Operator = DefaultOperator>
struct do_pattern_volume
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doPatternVolume;
  };
};

template<typename Operator = DefaultOperator>
struct do_pattern_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doPatternSkeleton;
  };
};

template<typename Operator = DefaultOperator>
struct do_pattern_boundary
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doPatternBoundary;
  };
};

template<typename Operator = DefaultOperator>
struct do_pattern_skeleton_or_boundary
{
  template<typename T>
  struct test {
    static const bool value =
      Operator::template ExtractType<T>::Type::doPatternSkeleton ||
      Operator::template ExtractType<T>::Type::doPatternBoundary;
  };
};

template<typename Operator = DefaultOperator>
struct do_pattern_coupling
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doPatternCoupling;
  };
};

template<typename Operator = DefaultOperator>
struct do_pattern_enriched_coupling
{
  template<typename T>
  struct test {
    static const bool value =
      Operator::template ExtractType<T>::Type::doPatternEnrichedCouplingToSubProblems ||
      Operator::template ExtractType<T>::Type::doPatternEnrichedCoupling;
  };
};

template<typename Operator = DefaultOperator>
struct do_pattern_volume_post_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doPatternVolumePostSkeleton;
  };
};

template<typename Operator = DefaultOperator>
struct do_alpha_volume
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doAlphaVolume;
  };
};

template<typename Operator = DefaultOperator>
struct do_alpha_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doAlphaSkeleton;
  };
};

template<typename Operator = DefaultOperator>
struct do_alpha_boundary
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doAlphaBoundary;
  };
};

template<typename Operator = DefaultOperator>
struct do_alpha_skeleton_or_boundary
{
  template<typename T>
  struct test {
    static const bool value =
      Operator::template ExtractType<T>::Type::doAlphaSkeleton ||
      Operator::template ExtractType<T>::Type::doAlphaBoundary;
  };
};

template<typename Operator = DefaultOperator>
struct do_alpha_coupling
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doAlphaCoupling;
  };
};

template<typename Operator = DefaultOperator>
struct do_alpha_enriched_coupling
{
  template<typename T>
  struct test {
    static const bool value =
      Operator::template ExtractType<T>::Type::doAlphaEnrichedCouplingToSubProblems ||
      Operator::template ExtractType<T>::Type::doAlphaEnrichedCoupling;
  };
};

template<typename Operator = DefaultOperator>
struct do_alpha_volume_post_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doAlphaVolumePostSkeleton;
  };
};



template<typename Operator = DefaultOperator>
struct do_lambda_volume
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doLambdaVolume;
  };
};

template<typename Operator = DefaultOperator>
struct do_lambda_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doLambdaSkeleton;
  };
};

template<typename Operator = DefaultOperator>
struct do_lambda_boundary
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doLambdaBoundary;
  };
};

template<typename Operator = DefaultOperator>
struct do_lambda_skeleton_or_boundary
{
  template<typename T>
  struct test {
    static const bool value =
      Operator::template ExtractType<T>::Type::doLambdaSkeleton ||
      Operator::template ExtractType<T>::Type::doLambdaBoundary;
  };
};

template<typename Operator = DefaultOperator>
struct do_lambda_coupling
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doLambdaCoupling;
  };
};

template<typename Operator = DefaultOperator>
struct do_lambda_enriched_coupling
{
  template<typename T>
  struct test {
    static const bool value =
      Operator::template ExtractType<T>::Type::doLambdaEnrichedCouplingToSubProblems ||
      Operator::template ExtractType<T>::Type::doLambdaEnrichedCoupling;
  };
};

template<typename Operator = DefaultOperator>
struct do_lambda_volume_post_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Operator::template ExtractType<T>::Type::doLambdaVolumePostSkeleton;
  };
};

template<typename Operator = DefaultOperator>
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

template<typename Extractor = IdentityExtractor>
struct do_constraints_volume
{
  template<typename T>
  struct test {
    static const bool value = Extractor::template ExtractType<T>::Type::doVolume;
  };
};

template<typename Extractor = IdentityExtractor>
struct do_constraints_skeleton
{
  template<typename T>
  struct test {
    static const bool value = Extractor::template ExtractType<T>::Type::doSkeleton;
  };
};

template<typename Extractor = IdentityExtractor>
struct do_constraints_boundary
{
  template<typename T>
  struct test {
    static const bool value = Extractor::template ExtractType<T>::Type::doBoundary;
  };
};

template<typename Extractor = IdentityExtractor>
struct do_constraints_skeleton_or_boundary
{
  template<typename T>
  struct test {
    static const bool value =
      Extractor::template ExtractType<T>::Type::doSkeleton ||
      Extractor::template ExtractType<T>::Type::doBoundary;
  };
};


template<typename Extractor = IdentityExtractor>
struct do_constraints_processor
{
  template<typename T>
  struct test {
    static const bool value = Extractor::template ExtractType<T>::Type::doProcessor;
  };
};



} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_OPERATORFLAGTESTS_HH
