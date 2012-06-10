#ifndef STOKESDARCYCOUPLINGOPERATOR_HH
#define STOKESDARCYCOUPLINGOPERATOR_HH

#include <dune/pdelab/multidomain/couplingutilities.hh>
#include <dune/common/parametertree.hh>
#include <dune/pdelab/localoperator/stokesparameter.hh>

template<typename GV, typename RF>
struct StokesDarcyCouplingParameterTraits
  : public Dune::PDELab::NavierStokesParameterTraits<GV,RF>
{

  typedef Dune::FieldMatrix<RF,GV::dimension,GV::dimension> PermeabilityTensor;

};

/*
template<typename GV, typename RF>
class CouplingParameters
{

public:

  typedef StokesDarcyCouplingParameterTraits<GV, RF> Traits;

  double viscosity() const
  {
    return viscosity_;
  }

  double density() const
  {
    return density_;
  }

  double gravity() const
  {
    return gravity_;
  }

  double alpha() const
  {
    return alpha_;
  }

  double porosity() const
  {
    return porosity_;
  }

  double gamma() const
  {
    return gamma_;
  }

  typename Traits::PermeabilityTensor kabs() const
  {
    return _kabs;
  }

  double epsilon() const
  {
    return epsilon_;
  }

  CouplingParameters(const Dune::ParameterTree& params)
    : alpha_(params.get("coupling.alpha",1.0))
    , gamma_(params.get("coupling.gamma",1.0))
    , porosity_(params.get<double>("soil.porosity"))
    , gravity_(params.get("gravity",9.81))
    , density_(params.get<double>("fluid.density"))
    , viscosity_(params.get<double>("fluid.viscosity"))
    , epsilon_(params.get("epsilon",1e-8))
  {
    for (int i = 0; i < GV::dimension; ++i)
      for (int j = 0; j < GV::dimension; ++j)
        _kabs[i][j] = (i == j ? params.get<double>("soil.permeability") : 0.0);
  }

private:
  PermeabilityTensor _kabs;
  const double alpha_;
  const double gamma_;
  const double porosity_;
  const double gravity_;
  const double density_;
  const double viscosity_;
  const double epsilon_;

};
*/

#define NUMERIC_DIFF

template<typename Parameters>
class StokesDarcyCouplingOperator
  : public Dune::PDELab::MultiDomain::CouplingOperatorDefaultFlags
  , public Dune::PDELab::MultiDomain::FullCouplingPattern
#ifdef NUMERIC_DIFF
  , public Dune::PDELab::MultiDomain::NumericalJacobianCoupling<StokesDarcyCouplingOperator<Parameters> >
  , public Dune::PDELab::MultiDomain::NumericalJacobianApplyCoupling<StokesDarcyCouplingOperator<Parameters> >
#endif
{

public:

  static const bool doPatternCoupling = true;
  static const bool doAlphaCoupling = true;

  StokesDarcyCouplingOperator(const Parameters& params)
    :
#ifdef NUMERIC_DIFF
    Dune::PDELab::MultiDomain::NumericalJacobianCoupling<StokesDarcyCouplingOperator<Parameters> >(params.epsilon()),
    Dune::PDELab::MultiDomain::NumericalJacobianApplyCoupling<StokesDarcyCouplingOperator<Parameters> >(params.epsilon()),
#endif
    parameters(params)
  {}

  template<typename IG,
           typename StokesLFSU,typename StokesLFSV,
           typename DarcyLFSU, typename DarcyLFSV,
           typename X, typename R>
  void alpha_coupling(const IG& ig,
                      const StokesLFSU& stokeslfsu, const X& stokesx, const StokesLFSV& stokeslsfv,
                      const DarcyLFSU& darcylfsu, const X& darcyx, const DarcyLFSV& darcylfsv,
                      R& stokesr, R& darcyr) const
  {
    // dimensions
    const int dim = IG::dimension;

    // extract local function spaces
    typedef typename StokesLFSU::template Child<0>::Type LFSU_V_PFS;
    const LFSU_V_PFS& lfsu_v_pfs = stokeslfsu.template child<0>();

    typedef typename LFSU_V_PFS::template Child<0>::Type LFSU_V;
    const unsigned int vsize = lfsu_v_pfs.child(0).size();

    // domain and range field type
    typedef Dune::FiniteElementInterfaceSwitch<typename LFSU_V::Traits::FiniteElementType > FESwitch_V;
    typedef Dune::BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
    typedef typename BasisSwitch_V::DomainField DF;
    typedef typename BasisSwitch_V::RangeField RF;
    typedef typename BasisSwitch_V::Range RT_V;
    typedef typename LFSU_V::Traits::SizeType size_type;

    typedef typename StokesLFSU::template Child<1>::Type LFSU_P;
    //const LFSU_P& lfsu_p = stokeslfsu.template getChild<1>();
    //const unsigned int psize = lfsu_p.size();

    typedef Dune::FiniteElementInterfaceSwitch<typename LFSU_P::Traits::FiniteElementType > FESwitch_P;
    typedef Dune::BasisInterfaceSwitch<typename FESwitch_P::Basis > BasisSwitch_P;
    typedef typename BasisSwitch_P::Range RT_P;

    typedef Dune::FiniteElementInterfaceSwitch<typename DarcyLFSU::Traits::FiniteElementType > FESwitch_D;
    typedef Dune::BasisInterfaceSwitch<typename FESwitch_D::Basis > BasisSwitch_D;
    typedef typename BasisSwitch_D::Range RT_D;
    // typedef typename DarcyLFSU::Traits::FiniteElementType::
    //   Traits::LocalBasisType::Traits::JacobianType JacobianType_D;
    const unsigned int dsize = darcylfsu.size();

    typedef typename IG::Geometry::LocalCoordinate LC;
    typedef typename IG::Geometry::GlobalCoordinate GC;

    // select quadrature rule
    Dune::GeometryType gt = ig.geometry().type();
    const int qorder = 2 *
        std::max(
            FESwitch_V::basis(lfsu_v_pfs.template child(0).finiteElement()).order(),
            FESwitch_D::basis(darcylfsu.finiteElement()).order());

    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gt,qorder);

    // const typename IG::Element& darcyCell = ig.outsideElement();

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {

        const GC pos = ig.geometry().global(it->position());
        const GC stokesPos = ig.geometryInInside().global(it->position());
        const GC darcyPos = ig.geometryInOutside().global(it->position());

        const RF g = 0; //parameters.gravity();
        const RF alpha = parameters.alpha(ig,it->position());
        const RF nu = parameters.viscosity(ig,it->position());
        // const RF porosity = parameters.porosity(ig,it->position());
        // const RF gamma = parameters.gamma(ig,it->position());
        const RF rho = parameters.density(ig,it->position());
        const typename Parameters::Traits::PermeabilityTensor kabs = parameters.K(ig,it->position());
        RF tracePi = 0.0;
        for (int i = 0; i < dim; ++i)
          tracePi += kabs[i][i];
        tracePi *= nu/g/rho;

        // integration weight
        const RF factor = it->weight() * ig.geometry().integrationElement(it->position());

        std::vector<RT_V> v(vsize);
        FESwitch_V::basis(lfsu_v_pfs.child(0).finiteElement()).evaluateFunction(stokesPos,v);

        std::vector<RT_D> psi(dsize);
        FESwitch_D::basis(darcylfsu.finiteElement()).evaluateFunction(darcyPos,psi);

        // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
        std::vector<Dune::FieldMatrix<RF,1,dim> > gradpsi(dsize);
        BasisSwitch_D::gradient(FESwitch_D::basis(darcylfsu.finiteElement()),
            ig.outside()->geometry(), darcyPos, gradpsi);

        // calculate phi and grad phi
        RT_D phi = 0.0;
        GC gradphi(0.0);
        for (size_type i = 0; i < darcylfsu.size(); ++i)
          {
            phi += darcyx(darcylfsu,i) * psi[i];
            gradphi.axpy(darcyx(darcylfsu,i),gradpsi[i][0]);
          }

        Dune::FieldVector<RF,dim> u(0.0);
        const GC n = ig.unitOuterNormal(it->position());

        // calculate u
        for (int d = 0; d < dim; ++d)
          {
            const LFSU_V& lfsu_v = lfsu_v_pfs.child(d);
            // calculate d-th component of u
            for (size_type i = 0; i < lfsu_v.size(); ++i)
              u[d] += stokesx(lfsu_v,i) * v[i];
          }

        GC tangentialFlow(0.0);
        //kabs.mv(gradphi,tangentialFlow);
        //tangentialFlow /= porosity;
        tangentialFlow += u;
        // project into tangential plane
        GC scaledNormal = n;
        scaledNormal *= (tangentialFlow * n);
        tangentialFlow -= scaledNormal;

        const RF h1 = parameters.h1(ig,it->position(),n);
        const GC h2 = parameters.h2(ig,it->position(),n);
        const RF h3 = parameters.h3(ig,it->position(),n);

        GC normalStress = n;
        normalStress *= h1 + phi;

        tangentialFlow *= alpha / sqrt(1);

        normalStress += tangentialFlow;
        normalStress += h2;

        for (size_type i = 0; i < darcylfsu.size(); ++i)
          darcyr.accumulate(darcylfsv,i, ((u * n) + h3) * psi[i] * factor);

        for (int d = 0; d < dim; ++d)
          {
            const LFSU_V& lfsu_v = lfsu_v_pfs.child(d);
            for (size_type i = 0; i < lfsu_v.size(); ++i)
              {
                // stokesr.accumulate(lfsu_v,i, - rho * g * (phi - pos[dim-1]) * v[i] * n[d] * factor);
                stokesr.accumulate(lfsu_v,i, normalStress[d] * v[i] * factor);
              }
          }

      }
  }

  template<typename IG,
           typename StokesLFSU,typename StokesLFSV,
           typename DarcyLFSU, typename DarcyLFSV,
           typename X, typename M>
  void jacobian_coupling(const IG& ig,
                         const StokesLFSU& stokeslfsu, const X& stokesx, const StokesLFSV& stokeslfsv,
                         const DarcyLFSU& darcylfsu, const X& darcyx, const DarcyLFSV& darcylfsv,
                         M& jac_stokes_stokes, M& jac_stokes_darcy,
                         M& jac_darcy_stokes, M& jac_darcy_darcy) const
  {
    // dimensions
    const int dim = IG::dimension;

    // extract local function spaces
    typedef typename StokesLFSU::template Child<0>::Type LFSU_V_PFS;
    const LFSU_V_PFS& lfsu_v_pfs = stokeslfsu.template child<0>();

    typedef typename LFSU_V_PFS::template Child<0>::Type LFSU_V;
    const unsigned int vsize = lfsu_v_pfs.child(0).size();

    typedef typename StokesLFSV::template Child<0>::Type LFSV_V_PFS;
    const LFSV_V_PFS& lfsv_v_pfs = stokeslfsv.template child<0>();

    typedef typename LFSV_V_PFS::template Child<0>::Type LFSV_V;

    // domain and range field type
    typedef Dune::FiniteElementInterfaceSwitch<typename LFSU_V::Traits::FiniteElementType > FESwitch_V;
    typedef Dune::BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
    typedef typename BasisSwitch_V::DomainField DF;
    typedef typename BasisSwitch_V::RangeField RF;
    typedef typename BasisSwitch_V::Range RT_V;
    typedef typename LFSU_V::Traits::SizeType size_type;

    typedef typename StokesLFSU::template Child<1>::Type LFSU_P;

    typedef Dune::FiniteElementInterfaceSwitch<typename LFSU_P::Traits::FiniteElementType > FESwitch_P;
    typedef Dune::BasisInterfaceSwitch<typename FESwitch_P::Basis > BasisSwitch_P;
    typedef typename BasisSwitch_P::Range RT_P;

    typedef Dune::FiniteElementInterfaceSwitch<typename DarcyLFSU::Traits::FiniteElementType > FESwitch_D;
    typedef Dune::BasisInterfaceSwitch<typename FESwitch_D::Basis > BasisSwitch_D;
    typedef typename BasisSwitch_D::Range RT_D;
    // typedef typename DarcyLFSU::Traits::FiniteElementType::
    //   Traits::LocalBasisType::Traits::JacobianType JacobianType_D;
    const unsigned int dsize = darcylfsu.size();

    typedef typename IG::Geometry::LocalCoordinate LC;
    typedef typename IG::Geometry::GlobalCoordinate GC;

    // select quadrature rule
    Dune::GeometryType gt = ig.geometry().type();
    const int qorder = 2 *
        std::max(
            FESwitch_V::basis(lfsu_v_pfs.template child(0).finiteElement()).order(),
            FESwitch_D::basis(darcylfsu.finiteElement()).order());

    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gt,qorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {

        const GC pos = ig.geometry().global(it->position());
        const GC stokesPos = ig.geometryInInside().global(it->position());
        const GC darcyPos = ig.geometryInOutside().global(it->position());

        const RF g = 0; //parameters.gravity();
        const RF alpha = parameters.alpha(ig,it->position());
        const RF nu = parameters.viscosity(ig,it->position());
        // const RF porosity = parameters.porosity(ig,it->position());
        // const RF gamma = parameters.gamma(ig,it->position());
        const RF rho = parameters.density(ig,it->position());
        const typename Parameters::Traits::PermeabilityTensor kabs = parameters.K(ig,it->position());
        RF tracePi = 0.0;
        for (int i = 0; i < dim; ++i)
          tracePi += kabs[i][i];
        tracePi *= nu/g/rho;

        // integration weight
        const RF factor = it->weight() * ig.geometry().integrationElement(it->position());

        std::vector<RT_V> v(vsize);
        FESwitch_V::basis(lfsu_v_pfs.child(0).finiteElement()).evaluateFunction(stokesPos,v);

        std::vector<RT_D> psi(dsize);
        FESwitch_D::basis(darcylfsu.finiteElement()).evaluateFunction(darcyPos,psi);

        const GC n = ig.unitOuterNormal(it->position());

        for (int d = 0; d < dim; ++d)
          {
            const LFSU_V& lfsu_v = lfsu_v_pfs.child(d);
            for (size_type i = 0; i < darcylfsv.size(); ++i)
              for (size_type j = 0; j < lfsu_v.size(); ++j)
                jac_darcy_stokes.accumulate(darcylfsv,i,lfsu_v,j, v[j] * n[d] * psi[i] * factor);
          }

         for (int d = 0; d < dim; ++d)
          {
            const LFSV_V& lfsv_v = lfsv_v_pfs.child(d);

            for (size_type i = 0; i < lfsv_v.size(); ++i)
              for (size_type j = 0; j < darcylfsu.size(); ++j)
                {
                  jac_stokes_darcy.accumulate(lfsv_v,i,darcylfsu,j, psi[j] * v[i] * n[d] * factor);
                }


            for (size_type i = 0; i < lfsv_v.size(); ++i)
              for (int dd = 0; dd < dim; ++dd)
                {
                  const LFSU_V& lfsu_v = lfsu_v_pfs.child(dd);
                  const RF dim_factor = ((d == dd ? 1.0 : 0.0) - n[d] * n[dd]) * factor;
                  for (size_type j = 0; j < lfsu_v.size(); ++j)
                    {
                      jac_stokes_stokes.accumulate(lfsv_v,i,lfsu_v,j, alpha / sqrt(1) * v[j] * v[i] * dim_factor);
                    }
                }
          }

      }
  }


private:

  const Parameters& parameters;

};


#endif // STOKESDARCYCOUPLINGOPERATOR_HH
