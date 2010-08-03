#ifndef SIMPLETIMEOPERATOR_HH
#define SIMPLETIMEOPERATOR_HH

class SimpleTimeOperator
  : public Dune::PDELab::NumericalJacobianApplyVolume<SimpleTimeOperator>,
    public Dune::PDELab::NumericalJacobianVolume<SimpleTimeOperator>,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  // constructor remembers parameters
  SimpleTimeOperator (unsigned int intorder_=2)
    : intorder(intorder_) {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {

    // domain and range field type (assume both components have same RF)
    typedef typename LFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianType;
    typedef typename LFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename LFSU::Traits::SizeType size_type;

    // dimensions
    const int dim = EG::Geometry::dimension;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate basis functions on reference element
        std::vector<RangeType> phi(lfsu.size());
        lfsu.localFiniteElement().localBasis().evaluateFunction(it->position(),phi);

        // compute u_0, u_1 at integration point
        RF u=0.0;
        for (size_type i=0; i<lfsu.size(); i++) u += x[lfsu.localIndex(i)]*phi[i];

        // integration
        RF factor = it->weight() * eg.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsu.size(); i++)
          r[lfsu.localIndex(i)] += u*phi[i]*factor;
      }
  }
private:
  unsigned int intorder;
};

#endif // SIMPLETIMEOPERATOR_HH
