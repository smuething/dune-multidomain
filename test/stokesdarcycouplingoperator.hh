#ifndef STOKESDARCYCOUPLINGOPERATOR_HH
#define STOKESDARCYCOUPLINGOPERATOR_HH

#include <dune/pdelab/multidomain/couplingutilities.hh>
#include <dune/common/parametertree.hh>

template<typename GV>
class CouplingParameters
{

public:

  typedef Dune::FieldMatrix<double,GV::dimension,GV::dimension> PermeabilityTensor;

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

  PermeabilityTensor kabs() const
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


template<typename Parameters>
class StokesDarcyCouplingOperator
  : public Dune::PDELab::MultiDomain::CouplingOperatorDefaultFlags
  , public Dune::PDELab::MultiDomain::FullCouplingPattern
  , public Dune::PDELab::MultiDomain::NumericalJacobianCoupling<StokesDarcyCouplingOperator<Parameters> >
  , public Dune::PDELab::MultiDomain::NumericalJacobianApplyCoupling<StokesDarcyCouplingOperator<Parameters> >
{

public:

  static const bool doPatternCoupling = true;
  static const bool doAlphaCoupling = true;

  StokesDarcyCouplingOperator(const Parameters& params)
    : Dune::PDELab::MultiDomain::NumericalJacobianCoupling<StokesDarcyCouplingOperator<Parameters> >(params.epsilon())
    , Dune::PDELab::MultiDomain::NumericalJacobianApplyCoupling<StokesDarcyCouplingOperator<Parameters> >(params.epsilon())
    , parameters(params)
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
    const int dimw = IG::dimensionworld;

    // extract local function spaces
    typedef typename StokesLFSU::template Child<0>::Type LFSU_V_PFS;
    const LFSU_V_PFS& lfsu_v_pfs = stokeslfsu.template child<0>();

    typedef typename LFSU_V_PFS::template Child<0>::Type LFSU_V;
    const unsigned int vsize = lfsu_v_pfs.child(0).size();

    // domain and range field type
    typedef typename LFSU_V::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU_V::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RT_V;
    typedef typename LFSU_V::Traits::SizeType size_type;

    typedef typename StokesLFSU::template Child<1>::Type LFSU_P;
    //const LFSU_P& lfsu_p = stokeslfsu.template getChild<1>();
    //const unsigned int psize = lfsu_p.size();

    typedef typename LFSU_P::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU_P::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RT_P;

    typedef typename DarcyLFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RT_D;
    typedef typename DarcyLFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianType_D;
    const unsigned int dsize = darcylfsu.size();

    typedef typename IG::Geometry::LocalCoordinate LC;
    typedef typename IG::Geometry::GlobalCoordinate GC;

    // select quadrature rule
    Dune::GeometryType gt = ig.geometry().type();
    const int qorder = 2 * std::max(lfsu_v_pfs.template child(0).finiteElement().localBasis().order(),
                                    darcylfsu.finiteElement().localBasis().order());

    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gt,qorder);

    const typename IG::EntityPointer darcyCell = ig.outside();

    const RF g = parameters.gravity();
    const RF alpha = parameters.alpha();
    const RF nu = parameters.viscosity();
    const RF porosity = parameters.porosity();
    const RF gamma = parameters.gamma();
    const RF rho = parameters.density();
    const typename Parameters::PermeabilityTensor kabs = parameters.kabs();
    RF tracePi = 0.0;
    for (int i = 0; i < dim; ++i)
      tracePi += kabs[i][i];
    tracePi *= nu/g/rho;

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {

        const GC pos = ig.geometry().global(it->position());
        const GC stokesPos = ig.geometryInInside().global(it->position());
        const GC darcyPos = ig.geometryInOutside().global(it->position());

        // integration weight
        const RF factor = it->weight() * ig.geometry().integrationElement(it->position());

        std::vector<RT_V> v(vsize);
        lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateFunction(stokesPos,v);

        std::vector<RT_D> psi(dsize);
        darcylfsu.finiteElement().localBasis().evaluateFunction(darcyPos,psi);

        // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
        std::vector<JacobianType_D> js(dsize);
        darcylfsu.finiteElement().localBasis().evaluateJacobian(darcyPos,js);

        // transform gradient to real element
        const Dune::FieldMatrix<DF,dimw,dim> jac = darcyCell->geometry().jacobianInverseTransposed(darcyPos);
        std::vector<Dune::FieldVector<RF,dim> > gradpsi(dsize);
        for (size_type i=0; i<vsize; i++)
            jac.mv(js[i][0],gradpsi[i]);

        // calculate phi and grad phi
        RT_D phi = 0.0;
        GC gradphi(0.0);
        for (size_type i = 0; i < darcylfsu.size(); ++i)
          {
            phi += darcyx(darcylfsu,i) * psi[i];
            gradphi.axpy(darcyx(darcylfsu,i),gradpsi[i]);
          }

        Dune::FieldVector<RF,dim> u(0.0);
        const GC n = ig.unitOuterNormal(it->position());

        // calculate u
        for (int d = 0; d < dim; ++d)
          {
            const LFSU_V& lfsu_v = lfsu_v_pfs.child(d);
            // calculate d-th component of u
            for (size_type i = 0; i < lfsu_v.size(); ++i)
              u[d] += stokesx[lfsu_v.localIndex(i)] * v[i];
          }


        for (size_type i = 0; i < darcylfsu.size(); ++i)
          darcyr[darcylfsu.localIndex(i)] -= gamma * porosity * (u * n) * psi[i] * factor;

        Dune::FieldVector<RF,dim> tangentialFlow(0.0);
        kabs.mv(gradphi,tangentialFlow);
        tangentialFlow /= porosity;
        tangentialFlow += u;
        // project into tangential plane
        GC scaledNormal = n;
        scaledNormal *= (tangentialFlow * n);
        tangentialFlow -= scaledNormal;

        for (int d = 0; d < dim; ++d)
          {
            const LFSU_V& lfsu_v = lfsu_v_pfs.child(d);
            for (size_type i = 0; i < lfsu_v.size(); ++i)
              {
                stokesr[lfsu_v.localIndex(i)] -= rho * g * (phi - pos[dim-1]) * v[i] * n[d] * factor;
              }

            for (size_type i = 0; i < lfsu_v.size(); ++i)
              {
                stokesr[lfsu_v.localIndex(i)] += alpha * sqrt(dim) / sqrt(tracePi) * tangentialFlow[d] * v[i] * factor;
              }
          }

      }
  }

private:

  const Parameters& parameters;

};


#endif // STOKESDARCYCOUPLINGOPERATOR_HH
