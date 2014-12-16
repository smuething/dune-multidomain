// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef MULTIDOMAIN_TEST_NEUMANDIRICHLETCOUPLING_HH
#define MULTIDOMAIN_TEST_NEUMANDIRICHLETCOUPLING_HH

#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>

/**
   \todo update quadrature order to work with lfsv != lfsu
   \todo update alpha_* to work with lfsv != lfsu (./)
   \todo update jacobian_* to work with lfsv != lfsu
   \todo update caches to work with lfsv != lfsu
 */

namespace Dune {
  namespace PDELab {

    enum class CouplingMode
    {
      Neumann,
      Dirichlet,
      Both,
    };

    /** a local operator for solving the convection-diffusion equation with discontinuous Galerkin
     *
     * \f{align*}{
     *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \\
     *                                              u &=& g \mbox{ on } \partial\Omega_D \\
     *                (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
     *                        -(A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_O
     * \f}
     * Note:
     *  - This formulation is valid for velocity fields which are non-divergence free.
     *  - Outflow boundary conditions should only be set on the outflow boundary
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T>
    class ConvectionDiffusionDGNeumannDirichletCoupling
      : public Dune::PDELab::MultiDomain::CouplingOperatorDefaultFlags,
        public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>
    {
      enum { dim = T::Traits::GridViewType::dimension };

      typedef typename T::Traits::RangeFieldType Real;
      typedef typename ConvectionDiffusionBoundaryConditions::Type BCType;

    public:
      // pattern assembly flags
      // there is no additional pattern relative to the subproblems

      // residual assembly flags
      enum { doAlphaCoupling  = true };

      //! constructor: pass parameter object
      ConvectionDiffusionDGNeumannDirichletCoupling(
        CouplingMode coupling_mode,
        T& param_,
        ConvectionDiffusionDGMethod::Type method_=ConvectionDiffusionDGMethod::NIPG,
        ConvectionDiffusionDGWeights::Type weights_=ConvectionDiffusionDGWeights::weightsOff,
        Real alpha_=0.0,
        int intorderadd_=0)
        : _coupling_mode(coupling_mode)
        , param(param_)
        , method(method_)
        , weights(weights_)
        , alpha(alpha_)
        , intorderadd(intorderadd_)
        , quadrature_factor(2)
      {
        theta = 1.0;
        if (method==ConvectionDiffusionDGMethod::SIPG) theta = -1.0;
      }

      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<
        typename IG,
        typename LFSU, typename LFSV,
        typename RemoteLFSU, typename RemoteLFSV,
        typename X, typename R
        >
      void alpha_coupling(
        const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const RemoteLFSU& lfsu_r, const X& x_r, const RemoteLFSV& lfsv_r,
        R& r_s, const R& r_r) const
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
        const int order = std::max(
            std::max(lfsu_s.finiteElement().localBasis().order(),
                lfsu_r.finiteElement().localBasis().order()),
            std::max(lfsv_s.finiteElement().localBasis().order(),
                lfsv_r.finiteElement().localBasis().order())
            );
        const int intorder = intorderadd+quadrature_factor*order;

        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>&
          outside_local = Dune::ReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);
        typename T::Traits::PermTensorType A_s, A_n;
        A_s = param.A(*(ig.inside()),inside_local);
        A_n = param.A(*(ig.outside()),outside_local);

        // face diameter; this should be revised for anisotropic meshes?
        DF h_s, h_n;
        DF hmax_s = 0.;
        DF hmax_n = 0.;
        element_size(ig.inside()->geometry(),h_s,hmax_s);
        element_size(ig.outside()->geometry(),h_n,hmax_n);
        RF h_F = std::min(h_s,h_n);
        h_F = std::min(ig.inside()->geometry().volume(),ig.outside()->geometry().volume())/ig.geometry().volume(); // Houston!

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // transformation
        typename IG::Entity::Geometry::JacobianInverseTransposed jac;

        // tensor times normal
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        Dune::FieldVector<RF,dim> An_F_n;
        A_n.mv(n_F,An_F_n);

        // compute weights
        RF omega_s;
        RF omega_n;
        RF harmonic_average(0.0);
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          {
            RF delta_s = (An_F_s*n_F);
            RF delta_n = (An_F_n*n_F);
            omega_s = delta_n/(delta_s+delta_n+1e-20);
            omega_n = delta_s/(delta_s+delta_n+1e-20);
            harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
          }
        else
          {
            omega_s = omega_n = 0.5;
            harmonic_average = 1.0;
          }

        // get polynomial degree
        const int order_s = lfsu_s.finiteElement().localBasis().order();
        const int order_n = lfsu_r.finiteElement().localBasis().order();
        int degree = std::max( order_s, order_n );

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // exact normal
            const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(it->position());

            // position of quadrature point in local coordinates of elements
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());
            Dune::FieldVector<DF,dim> iplocal_n = ig.geometryInOutside().global(it->position());


            // evaluate basis functions
            std::vector<RangeType> psi_s(lfsv_s.size());
            lfsv_s.finiteElement().localBasis().evaluateFunction(iplocal_s,psi_s);

            // integration factor
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());

            // evaluate velocity field and upwinding, assume H(div) velocity field => may choose any side
            typename T::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
            RF normalflux = b*n_F_local;

            std::vector<RangeType> phi_r(lfsu_r.size());
            lfsu_r.finiteElement().localBasis().evaluateFunction(iplocal_n,phi_r);
            RF u_r=0.0;
            for (size_type i=0; i<lfsu_r.size(); i++)
              u_r += x_r(lfsu_r,i)*phi_r[i];


            if (_coupling_mode == CouplingMode::Neumann || _coupling_mode == CouplingMode::Both)
              {
                // evaluate flux on neighbor

                std::vector<JacobianType> gradphi_r(lfsu_r.size());
                lfsu_r.finiteElement().localBasis().evaluateJacobian(iplocal_n,gradphi_r);

                jac = ig.outside()->geometry().jacobianInverseTransposed(iplocal_n);
                std::vector<Dune::FieldVector<RF,dim> > tgradphi_r(lfsu_r.size());
                for (size_type i=0; i<lfsu_r.size(); i++) jac.mv(gradphi_r[i][0],tgradphi_r[i]);

                Dune::FieldVector<RF,dim> gradu_r(0.0);
                for (size_type i=0; i<lfsu_r.size(); i++)
                  gradu_r.axpy(x_r(lfsu_r,i),tgradphi_r[i]);

                // flux in normal direction
                RF j = normalflux * u_r - gradu_r * n_F_local;
                //j = 0;
                // integrate
                for (size_type i=0; i<lfsv_s.size(); i++)
                  r_s.accumulate(lfsv_s,i, j * psi_s[i] * factor);

                if (_coupling_mode == CouplingMode::Neumann)
                  continue;
              }

            std::vector<RangeType> phi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);

            // evaluate u
            RF u_s=0.0;
            for (size_type i=0; i<lfsu_s.size(); i++)
              u_s += x_s(lfsu_s,i)*phi_s[i];

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
            std::vector<JacobianType> gradpsi_s(lfsv_s.size());
            lfsv_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradpsi_s);

            // transform gradients of shape functions to real element
            jac = ig.inside()->geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
            std::vector<Dune::FieldVector<RF,dim> > tgradpsi_s(lfsv_s.size());
            for (size_type i=0; i<lfsv_s.size(); i++) jac.mv(gradpsi_s[i][0],tgradpsi_s[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu_s(0.0);
            for (size_type i=0; i<lfsu_s.size(); i++)
              gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);

            RF omegaup_s, omegaup_r;
            if (normalflux>=0.0)
              {
                omegaup_s = 1.0;
                omegaup_r = 0.0;
              }
            else
              {
                omegaup_s = 0.0;
                omegaup_r = 1.0;
              }

            // convection term
            RF term1 = (omegaup_s*u_s + omegaup_r*u_r) * normalflux *factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsu_s,i,term1 * psi_s[i]);

            // diffusion term
            RF term2 =  -factor * (An_F_s*gradu_s);
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsv_s,i,term2 * psi_s[i]);

            // (non-)symmetric IP term
            RF term3 = (u_s-u_r) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsv_s,i,term3 * theta * (An_F_s*tgradpsi_s[i]));

            // standard IP term integral
            RF term4 = penalty_factor * (u_s-u_r) * factor;
            for (size_type i=0; i<lfsv_s.size(); i++)
              r_s.accumulate(lfsv_s,i,term4 * psi_s[i]);
          }
      }

      template<
        typename IG,
        typename LFSU, typename LFSV,
        typename RemoteLFSU, typename RemoteLFSV,
        typename X, typename M
        >
      void jacobian_coupling(
        const IG& ig,
        const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
        const RemoteLFSU& lfsu_r, const X& x_r, const RemoteLFSV& lfsv_r,
        M& mat_ss, const M& mat_sr,
        const M& mat_rs, const M& mat_rr) const
      {

        if (_coupling_mode == CouplingMode::Neumann)
          return;

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
        const int order = std::max(
            std::max(lfsu_s.finiteElement().localBasis().order(),
                lfsu_r.finiteElement().localBasis().order()),
            std::max(lfsv_s.finiteElement().localBasis().order(),
                lfsv_r.finiteElement().localBasis().order())
            );
        const int intorder = intorderadd+quadrature_factor*order;

        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>&
          outside_local = Dune::ReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);
        typename T::Traits::PermTensorType A_s, A_n;
        A_s = param.A(*(ig.inside()),inside_local);
        A_n = param.A(*(ig.outside()),outside_local);

        // face diameter; this should be revised for anisotropic meshes?
        DF h_s, h_n;
        DF hmax_s = 0., hmax_n = 0.;
        element_size(ig.inside()->geometry(),h_s,hmax_s);
        element_size(ig.outside()->geometry(),h_n,hmax_n);
        RF h_F = std::min(h_s,h_n);
        h_F = std::min(ig.inside()->geometry().volume(),ig.outside()->geometry().volume())/ig.geometry().volume(); // Houston!

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // transformation
        typename IG::Entity::Geometry::JacobianInverseTransposed jac;

        // tensor times normal
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        Dune::FieldVector<RF,dim> An_F_n;
        A_n.mv(n_F,An_F_n);

        // compute weights
        RF omega_s;
        RF omega_n;
        RF harmonic_average(0.0);
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          {
            RF delta_s = (An_F_s*n_F);
            RF delta_n = (An_F_n*n_F);
            omega_s = delta_n/(delta_s+delta_n+1e-20);
            omega_n = delta_s/(delta_s+delta_n+1e-20);
            harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
          }
        else
          {
            omega_s = omega_n = 0.5;
            harmonic_average = 1.0;
          }

        // get polynomial degree
        const int order_s = lfsu_s.finiteElement().localBasis().order();
        const int order_n = lfsu_r.finiteElement().localBasis().order();
        int degree = std::max( order_s, order_n );

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // exact normal
            const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(it->position());

            // position of quadrature point in local coordinates of elements
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());

            // evaluate basis functions
            std::vector<RangeType> phi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateFunction(iplocal_s,phi_s);

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> gradphi_s(lfsu_s.size());
            lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);

            // transform gradients of shape functions to real element
            jac = ig.inside()->geometry().jacobianInverseTransposed(iplocal_s);
            std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);

            // evaluate velocity field and upwinding, assume H(div) velocity field => may choose any side
            typename T::Traits::RangeType b = param.b(*(ig.inside()),iplocal_s);
            RF normalflux = b*n_F_local;
            RF omegaup_s, omegaup_r;
            if (normalflux>=0.0)
              {
                omegaup_s = 1.0;
                omegaup_r = 0.0;
              }
            else
              {
                omegaup_s = 0.0;
                omegaup_r = 1.0;
              }

            // integration factor
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());
            RF ipfactor = penalty_factor * factor;

            // do all terms in the order: I convection, II diffusion, III consistency, IV ip
            for (size_type j=0; j<lfsu_s.size(); j++) {
              RF temp1 = -(An_F_s*tgradphi_s[j])*factor;
              for (size_type i=0; i<lfsu_s.size(); i++) {
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,omegaup_s * phi_s[j] * normalflux *factor * phi_s[i]);
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,temp1 * phi_s[i]);
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * factor * theta * (An_F_s*tgradphi_s[i]));
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * ipfactor * phi_s[i]);
              }
            }
          }
      }


      //! set time in parameter class
      void setTime (double t)
      {
        param.setTime(t);
      }

    private:
      CouplingMode _coupling_mode;
      T& param;  // two phase parameter class
      ConvectionDiffusionDGMethod::Type method;
      ConvectionDiffusionDGWeights::Type weights;
      Real alpha, beta;
      int intorderadd;
      int quadrature_factor;
      Real theta;

      template<class GEO>
      void element_size (const GEO& geo, typename GEO::ctype& hmin, typename GEO::ctype hmax) const
      {
        typedef typename GEO::ctype DF;
        hmin = 1.0E100;
        hmax = -1.0E00;
        const int dim = GEO::coorddimension;
        if (dim==1)
          {
            Dune::FieldVector<DF,dim> x = geo.corner(0);
            x -= geo.corner(1);
            hmin = hmax = x.two_norm();
            return;
          }
        else
          {
            Dune::GeometryType gt = geo.type();
            for (int i=0; i<Dune::ReferenceElements<DF,dim>::general(gt).size(dim-1); i++)
              {
                Dune::FieldVector<DF,dim> x = geo.corner(Dune::ReferenceElements<DF,dim>::general(gt).subEntity(i,dim-1,0,dim));
                x -= geo.corner(Dune::ReferenceElements<DF,dim>::general(gt).subEntity(i,dim-1,1,dim));
                hmin = std::min(hmin,x.two_norm());
                hmax = std::max(hmax,x.two_norm());
              }
            return;
          }
      }
    };
  }
}

#endif // MULTIDOMAIN_TEST_NEUMANDIRICHLETCOUPLING_HH
