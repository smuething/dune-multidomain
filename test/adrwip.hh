// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PM_ADRWIP_HH
#define DUNE_PM_ADRWIP_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

#include"adrparam.hh"

namespace Dune {
  namespace PM {

    struct ADRWIPMethod
    {
      enum Type { NIPG, SIPG, OBB };
    };

    struct ADRWIPWeights
    {
      enum Type { weightsOn, weightsOff };
    };

    /** \brief Weighted interior penalty discontinous Galerkin scheme

        Param : parameter class, see above
        cylindrical : if true (default) use use cylindrical coordinates (assumes dim==2)
    */
    template<typename Param>
    class ADRWIP
      : public Dune::PDELab::NumericalJacobianVolume<ADRWIP<Param> >,
        public Dune::PDELab::NumericalJacobianApplyVolume<ADRWIP<Param> >,
        public Dune::PDELab::NumericalJacobianSkeleton<ADRWIP<Param> >,
        public Dune::PDELab::NumericalJacobianApplySkeleton<ADRWIP<Param> >,
        public Dune::PDELab::NumericalJacobianBoundary<ADRWIP<Param> >,
        public Dune::PDELab::NumericalJacobianApplyBoundary<ADRWIP<Param> >,
        public Dune::PDELab::FullSkeletonPattern,
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::LocalOperatorDefaultFlags,
        public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename Param::Traits::RangeFieldType>
    {
      enum { dim = Param::Traits::GridViewType::dimension };

      typedef typename Param::Traits::RangeFieldType Real;
      typedef typename ADRBoundaryConditions::Type BCType;

      // Real reg (Real x) const
      // {
      //   const Real eps=1e-4;
      //   if (x>=1.0) return x;
      //   // return x*x*( (2*eps-1)*x + 2-3*eps) + eps;
      //   return eps + (1.0-eps)*x;
      // }

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume  = true };

      //! constructor: pass parameter object
      ADRWIP (Param& param_,
              ADRWIPMethod::Type method_=ADRWIPMethod::OBB,
              ADRWIPWeights::Type weights_=ADRWIPWeights::weightsOff,
              Real alpha_=0.0,
              Real beta_=0.0,
              int quadrature_factor_=2)
        : Dune::PDELab::NumericalJacobianVolume<ADRWIP<Param> >(1e-7),
          Dune::PDELab::NumericalJacobianApplyVolume<ADRWIP<Param> >(1e-7),
          Dune::PDELab::NumericalJacobianSkeleton<ADRWIP<Param> >(1e-7),
          Dune::PDELab::NumericalJacobianApplySkeleton<ADRWIP<Param> >(1e-7),
          Dune::PDELab::NumericalJacobianBoundary<ADRWIP<Param> >(1e-7),
          Dune::PDELab::NumericalJacobianApplyBoundary<ADRWIP<Param> >(1e-7),
          param(param_), method(method_), weights(weights_),
          alpha(alpha_), beta(beta_), quadrature_factor(quadrature_factor_)
      {
        epsilon = 1.0;
        if (method==ADRWIPMethod::SIPG) epsilon = -1.0;
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
		// domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;
        const int intorder = quadrature_factor*lfsu.finiteElement().localBasis().order();

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        typename Param::Traits::PermTensorType tensor;
        Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
        tensor = param.K(eg.entity(),localcenter);

        // transformation
        Dune::FieldMatrix<DF,dimw,dim> jac;

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> js(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // transform gradients of shape functions to real element
            jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            for (size_type i=0; i<lfsu.size(); i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (size_type i=0; i<lfsu.size(); i++)
              gradu.axpy(x(lfsu,i),gradphi[i]);

            // compute K * gradient of u
            Dune::FieldVector<RF,dim> Kgradu(0.0);
            tensor.umv(gradu,Kgradu);

            // evaluate velocity field
            typename Param::Traits::RangeType b = param.b(eg.entity(),it->position());

            // evaluate reaction term
            Real a = param.a(eg.entity(),it->position());

            // evaluate nonlinear diffusion term
            Real w = param.w(eg.entity(),it->position(),u);

            // integrate (K grad u - bu)*grad phi_i + a*u*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsu.size(); i++)
              r.accumulate(lfsu,i,( w*(Kgradu*gradphi[i]) - b*gradphi[i] + a*u*phi[i] )*factor);
          }
      }

      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
		// domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSV::Traits::SizeType size_type;

        // dimensions
        const int dim = IG::dimension;
        const int intorder = quadrature_factor*std::max(lfsu_s.finiteElement().localBasis().order(),
                                        lfsu_n.finiteElement().localBasis().order());

        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>&
          outside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);
        typename Param::Traits::PermTensorType K_s, K_n;
        K_s = param.K(*(ig.inside()),inside_local);
        K_n = param.K(*(ig.outside()),outside_local);

        // face diameter
        RF h_F = 0.0;
        if (dim==1) h_F = 1.0; else
          {
            Dune::FieldVector<DF,dim> x0 = ig.geometry().corner(0);
            for (int i=1; i<ig.geometry().corners(); i++)
              {
                Dune::FieldVector<DF,dim> x = ig.geometry().corner(i);
                x -= x0;
                h_F = std::max(h_F,x.two_norm());
              }
          }

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // compute weights
            const Dune::FieldVector<DF,dim> n_F = ig.unitOuterNormal(it->position());
            Dune::FieldVector<RF,dim> Kn_F_s;
            K_s.mv(n_F,Kn_F_s);
            RF delta_s = (Kn_F_s*n_F);
            Dune::FieldVector<RF,dim> Kn_F_n;
            K_n.mv(n_F,Kn_F_n);
            RF delta_n = (Kn_F_n*n_F);
            RF omega_s;
            RF omega_n;
            if (weights==ADRWIPWeights::weightsOn)
              {
                omega_s = delta_n/(delta_s+delta_n);
                omega_n = delta_s/(delta_s+delta_n);
              }
            else
              omega_s = omega_n = 0.5;

            // position of quadrature point in local coordinates of elements
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());
            Dune::FieldVector<DF,dim> iplocal_n = ig.geometryInOutside().global(it->position());

            // evaluate basis functions
            std::vector<RangeType> phi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);
            std::vector<RangeType> phi_n(lfsu_n.size());
            lfsu_n.finiteElement().localBasis().evaluateFunction(iplocal_n,phi_n);

            // evaluate u
            RF u_s=0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              u_s += x_s(lfsu_s,i)*phi_s[i];
            RF u_n=0.0;
            for (size_type i=0; i<lfsu_n.size(); i++)
              u_n += x_n(lfsu_n,i)*phi_n[i];

            // evaluate nonlinear diffusion term
            Real w_s = param.w(*(ig.inside()),iplocal_s,u_s);
            Real w_n = param.w(*(ig.outside()),iplocal_n,u_n);

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
            std::vector<JacobianType> gradphi_n(lfsu_n.size());
            lfsu_n.finiteElement().localBasis().evaluateJacobian(iplocal_n,gradphi_n);

            // transform gradients of shape functions to real element
            jac = ig.inside()->geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
            jac = ig.outside()->geometry().jacobianInverseTransposed(iplocal_n);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_n(lfsu_n.size());
            for (size_type i=0; i<lfsu_n.size(); i++) jac.mv(gradphi_n[i][0],tgradphi_n[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu_s(0.0);
            for (size_type i=0; i<lfsu_s.size(); i++)
              gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);
            Dune::FieldVector<RF,dim> gradu_n(0.0);
            for (size_type i=0; i<lfsu_n.size(); i++)
              gradu_n.axpy(x_n(lfsu_n,i),tgradphi_n[i]);

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            typename Param::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
            RF normalflux = b*n_F;
            RF omegaup_s, omegaup_n;
            if (normalflux>=0.0)
              {
                omegaup_s = 1.0;
                omegaup_n = 0.0;
              }
            else
              {
                omegaup_s = 0.0;
                omegaup_n = 1.0;
              }

            // integration factor
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());

            // convection term
            RF term1 = /*(omegaup_s*u_s + omegaup_n*u_n) * */ normalflux *factor;
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s.accumulate(lfsu_s,i,term1 * phi_s[i]);
            for (size_type i=0; i<lfsu_n.size(); i++)
              r_n.accumulate(lfsu_n,i,-term1 * phi_n[i]);

            // diffusion term
            RF diff_flux = -(omega_s*(Kn_F_s*gradu_s) + omega_n*(Kn_F_n*gradu_n));
            RF w_upwind;
            if (diff_flux>=0) w_upwind = w_s; else w_upwind = w_n;
            RF term2 =  w_upwind * diff_flux * factor;
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s.accumulate(lfsu_s,i,term2 * phi_s[i]);
            for (size_type i=0; i<lfsu_n.size(); i++)
              r_n.accumulate(lfsu_n,i,-term2 * phi_n[i]);

            // (non-)symmetric IP term
            RF term3 = (u_s-u_n) * factor;
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s.accumulate(lfsu_s,i,term3 * epsilon * omega_s * w_upwind * (Kn_F_s*tgradphi_s[i]));
            for (size_type i=0; i<lfsu_n.size(); i++)
              r_n.accumulate(lfsu_n,i,term3 * epsilon * omega_n * w_upwind * (Kn_F_n*tgradphi_n[i]));

            // standard IP term
            RF term4;
            if (weights==ADRWIPWeights::weightsOn)
              {
                RF gamma = alpha*w_upwind*(delta_s*delta_n/(delta_s+delta_n))/h_F;
                term4 = ((u_s-u_n) * gamma + beta*0.5*std::abs(normalflux)) * factor;
              }
            else
              {
                RF gamma = alpha/h_F;
                term4 = (u_s-u_n) * gamma * factor;
              }
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s.accumulate(lfsu_s,i,term4 * phi_s[i]);
            for (size_type i=0; i<lfsu_n.size(); i++)
              r_n.accumulate(lfsu_n,i,-term4 * phi_n[i]);
          }
      }

      // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
		// domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSV::Traits::SizeType size_type;

        // dimensions
        const int dim = IG::dimension;
        const int intorder = quadrature_factor*lfsu_s.finiteElement().localBasis().order();

        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        typename Param::Traits::PermTensorType K_s;
        K_s = param.K(*(ig.inside()),inside_local);

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // face diameter
        RF h_F = 0.0;
        if (dim==1) h_F = 1.0; else
          {
            Dune::FieldVector<DF,dim> x0 = ig.geometry().corner(0);
            for (int i=1; i<ig.geometry().corners(); i++)
              {
                Dune::FieldVector<DF,dim> x = ig.geometry().corner(i);
                x -= x0;
                h_F = std::max(h_F,x.two_norm());
              }
          }

        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

        // evaluate boundary condition
        const Dune::FieldVector<DF,dim-1>
          face_local = Dune::GenericReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        BCType bctype = param.bc(ig.intersection(),face_local);

        if (bctype == ADRBoundaryConditions::None)
          return;

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // compute weights
            const Dune::FieldVector<DF,dim> n_F = ig.unitOuterNormal(it->position());
            Dune::FieldVector<RF,dim> Kn_F_s;
            K_s.mv(n_F,Kn_F_s);
            RF delta_s = (Kn_F_s*n_F);

            // position of quadrature point in local coordinates of elements
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());

            // evaluate basis functions
            std::vector<RangeType> phi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);

            // integration factor
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());

            if (bctype == ADRBoundaryConditions::Flux)
              {
                // evaluate flux boundary condition
                RF j = param.j(ig.intersection(),it->position());

                // integrate
                for (size_type i=0; i<lfsv_s.size(); i++)
                  r_s.accumulate(lfsv_s,i,j * phi_s[i] * factor);

                continue;
              }

            // evaluate u
            RF u_s=0.0; for (size_type i=0; i<lfsu_s.size(); i++) u_s += x_s(lfsu_s,i)*phi_s[i];

            // evaluate nonlinear diffusion term
            Real w_s = param.w(*(ig.inside()),iplocal_s,u_s);

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);

            // transform gradients of shape functions to real element
            jac = ig.inside()->geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu_s(0.0);
            for (size_type i=0; i<lfsu_s.size(); i++) gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            typename Param::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
            RF normalflux = b*n_F;

            // convection term
            RF term1 = normalflux * factor;
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s.accumulate(lfsu_s,i,term1 * phi_s[i]);

            if (bctype == ADRBoundaryConditions::Outflow) continue;

            // evaluate Dirichlet boundary condition
            RF g = param.g(ig.intersection(),it->position());

            // interior penalty parameter
            RF gamma;
            if (weights==ADRWIPWeights::weightsOn)
              gamma = alpha*w_s*delta_s/h_F + beta*0.5*std::abs(normalflux)/(u_s-g);
            else
              gamma = alpha/h_F;

            // diffusion term
            RF term2 =  w_s*(Kn_F_s*gradu_s) * factor;
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s.accumulate(lfsu_s,i,-term2 * phi_s[i]);

            // (non-)symmetric IP term
            RF term3 = (u_s-g) * factor;
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s.accumulate(lfsu_s,i,term3 * epsilon * w_s * (Kn_F_s*tgradphi_s[i]));

            // standard IP term
            RF term4 = (u_s-g) * gamma * factor;
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s.accumulate(lfsu_s,i,term4 * phi_s[i]);
          }
      }

 	  // volume integral depending only on test functions
	  template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
		// domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int intorder = 2*lfsv.finiteElement().localBasis().order();

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate shape functions
            std::vector<RangeType> phi(lfsv.size());
            lfsv.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate right hand side parameter function
            Real f;
            f = param.f(eg.entity(),it->position());

            // integrate f
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsv.size(); i++)
              r.accumulate(lfsv,i,-f*phi[i]*factor);
          }
      }

      //! set time for subsequent evaluation
      void setTime (typename Param::Traits::RangeFieldType t)
      {
        time = t;
        param.setTime(t);
      }

    private:
      Param& param;  // two phase parameter class
      ADRWIPMethod::Type method;
      ADRWIPWeights::Type weights;
      Real alpha, beta;
      int quadrature_factor;
      Real epsilon;
      mutable Real time;
    };


    /** a local operator for the mass operator (L_2 integral)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     */
	class ADRWIPTemporal : public Dune::PDELab::NumericalJacobianApplyVolume<ADRWIPTemporal>,
                           public Dune::PDELab::FullVolumePattern,
                           public Dune::PDELab::LocalOperatorDefaultFlags,
                           public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };

	  // residual assembly flags
      enum { doAlphaVolume = true };

      ADRWIPTemporal ()
        : Dune::PDELab::NumericalJacobianApplyVolume<ADRWIPTemporal>(1e-7)
      {}

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
		// domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;
        //const int dimw = EG::Geometry::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const int intorder = 2*lfsu.finiteElement().localBasis().order();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // u*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsu.size(); i++)
              r.accumulate(lfsu,i,u*phi[i]*factor);
          }
	  }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
	  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
		// domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;
        //const int dimw = EG::Geometry::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const int intorder = 2*lfsu.finiteElement().localBasis().order();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // integrate phi_j*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat.accumulate(lfsu,i,lfsu,j,phi[j]*phi[i]*factor);
          }
      }
	};


  }
}

#endif
