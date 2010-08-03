// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef INSTATIONARY_POISSON_OPERATOR_HH
#define INSTATIONARY_POISSON_OPERATOR_HH

#include<vector>

#include <dune/pdelab/localoperator/poisson.hh>

template<typename F, typename B, typename J, int qorder=2>
class InstationaryPoisson : public Dune::PDELab::NumericalJacobianApplyVolume<InstationaryPoisson<F,B,J,qorder> >,
                            public Dune::PDELab::NumericalJacobianVolume<InstationaryPoisson<F,B,J,qorder> >,
                            public Dune::PDELab::FullVolumePattern,
                            public Dune::PDELab::LocalOperatorDefaultFlags,
                            public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doLambdaVolume = true };
  enum { doLambdaBoundary = true };

  InstationaryPoisson (F& f_, B& b_, J& j_)
    : f(f_), b(b_), j(j_)
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // domain and range field type
    typedef typename LFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianType;

    f.setTime(this->getTime());
    j.setTime(this->getTime());

    // dimensions
    const int dim = EG::Geometry::dimension;
    const int dimw = EG::Geometry::dimensionworld;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
        std::vector<JacobianType> js(lfsu.size());
        lfsu.localFiniteElement().localBasis().evaluateJacobian(it->position(),js);

        // transform gradient to real element
        const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
        for (size_t i=0; i<lfsu.size(); i++)
          {
            gradphi[i] = 0.0;
            jac.umv(js[i][0],gradphi[i]);
          }

        // compute gradient of u
        Dune::FieldVector<RF,dim> gradu(0.0);
        for (size_t i=0; i<lfsu.size(); i++)
          gradu.axpy(x[i],gradphi[i]);

        // integrate grad u * grad phi_i
        RF factor = it->weight() * eg.geometry().integrationElement(it->position());
        for (size_t i=0; i<lfsu.size(); i++)
          r[i] += (gradu*gradphi[i])*factor;
      }
  }

  // volume integral depending only on test functions
  template<typename EG, typename LFSV, typename R>
  void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
  {
    // domain and range field type
    typedef typename LFSV::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSV::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSV::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;

    f.setTime(this->getTime());
    j.setTime(this->getTime());

    // dimensions
    const int dim = EG::Geometry::dimension;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate shape functions
        std::vector<RangeType> phi(lfsv.size());
        lfsv.localFiniteElement().localBasis().evaluateFunction(it->position(),phi);

        // evaluate right hand side parameter function
        typename F::Traits::RangeType y;
        f.evaluate(eg.entity(),it->position(),y);

        // integrate f
        RF factor = it->weight() * eg.geometry().integrationElement(it->position());
        for (size_t i=0; i<lfsv.size(); i++)
          r[i] -= y*phi[i]*factor;
      }
  }

  // boundary integral independen of ansatz functions
  template<typename IG, typename LFSV, typename R>
  void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
  {
    // domain and range field type
    typedef typename LFSV::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSV::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSV::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;

    f.setTime(this->getTime());
    j.setTime(this->getTime());

    // dimensions
    const int dim = IG::dimension;

    // select quadrature rule
    Dune::GeometryType gtface = ig.geometryInInside().type();
    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

    // loop over quadrature points and integrate normal flux
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate boundary condition type
        typename B::Traits::RangeType bctype;
        b.evaluate(ig,it->position(),bctype);

        // skip rest if we are on Dirichlet boundary
        if (isDirichlet(bctype)) continue;

        // position of quadrature point in local coordinates of element
        Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

        // evaluate test shape functions
        std::vector<RangeType> phi(lfsv.size());
        lfsv.localFiniteElement().localBasis().evaluateFunction(local,phi);

        // evaluate flux boundary condition
        typename J::Traits::RangeType y;
        j.evaluate(*(ig.inside()),local,y);

        // integrate J
        RF factor = it->weight()*ig.geometry().integrationElement(it->position());
        for (size_t i=0; i<lfsv.size(); i++)
          r[i] += y*phi[i]*factor;
      }
  }

protected:
  F& f;
  B& b;
  J& j;
};

#endif
