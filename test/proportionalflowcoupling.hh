#ifndef PROPORTIONALFLOWCOUPLING_HH
#define PROPORTIONALFLOWCOUPLING_HH


#include <dune/pdelab/multidomain/couplingutilities.hh>
#include <dune/pdelab/localoperator/pattern.hh>

class ProportionalFlowCoupling :
  public Dune::PDELab::MultiDomain::NumericalJacobianCoupling<ProportionalFlowCoupling>,
  public Dune::PDELab::MultiDomain::NumericalJacobianApplyCoupling<ProportionalFlowCoupling>,
  public Dune::PDELab::MultiDomain::FullCouplingPattern
{

public:

  ProportionalFlowCoupling(double intensity)
    : _intensity(intensity)
  {}

  static const bool doAlphaCoupling = true;
  static const bool doPatternCoupling = true;

  template<typename IG, typename LFSU1, typename LFSU2, typename X, typename LFSV1, typename LFSV2,
           typename R>
  void alpha_coupling
  ( const IG& ig,
    const LFSU1& lfsu_s, const X& x_s, const LFSV1& lfsv_s,
    const LFSU2& lfsu_n, const X& x_n, const LFSV2& lfsv_n,
    R& r_s, R& r_n) const
  {
    // domain and range field type
    typedef typename LFSU1::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU1::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU1::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;

    const int intorder = 4;

    typedef typename LFSU1::Traits::SizeType size_type;

    // dimensions
    const int dim = IG::dimension;

    // select quadrature rule
    Dune::GeometryType gtface = ig.geometryInInside().type();
    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

    // loop over quadrature points and integrate normal flux
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // position of quadrature point in local coordinates of element
        Dune::FieldVector<DF,dim> local1 = ig.geometryInInside().global(it->position());
        Dune::FieldVector<DF,dim> local2 = ig.geometryInOutside().global(it->position());

        // evaluate ansatz shape functions (assume Galerkin for now)
        std::vector<RangeType> phi1(lfsv_s.size());
        lfsv_s.localFiniteElement().localBasis().evaluateFunction(local1,phi1);

        std::vector<RangeType> phi2(lfsv_n.size());
        lfsv_n.localFiniteElement().localBasis().evaluateFunction(local2,phi2);

        RF u_s(0.0);
        for (size_t i=0; i<lfsu_s.size(); i++)
          u_s += x_s[lfsu_s.localIndex(i)] * phi1[i];

        RF u_n(0.0);
        for (size_t i=0; i<lfsu_n.size(); i++)
          u_n += x_n[lfsu_n.localIndex(i)] * phi2[i];

        RF u_diff = _intensity*(u_s - u_n);

        // integrate J
        RF factor = it->weight()*ig.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsv_s.size(); i++)
          r_s[lfsu_s.localIndex(i)] += u_diff*phi1[i]*factor;
        for (size_type i=0; i<lfsv_n.size(); i++)
          r_n[lfsu_n.localIndex(i)] -= u_diff*phi2[i]*factor;

      }
  }

private:
  const double _intensity;

};

#endif // PROPORTIONALFLOWCOUPLING_HH
