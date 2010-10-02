#ifndef STOKESDARCYCOUPLINGOPERATOR_HH
#define STOKESDARCYCOUPLINGOPERATOR_HH

#include <dune/pdelab/multidomain/couplingutilities.hh>

template<typename GV>
class CouplingParameters
{

public:

  typedef Dune::FieldMatrix<double,GV::dimension,GV::dimension> PermeabilityTensor;

  double viscosity() const
  {
    return 1.0;
  }

  double gravity() const
  {
    return 1.0;
  }

  double alpha() const
  {
    return 1.0;
  }

  double porosity() const
  {
    return 1.0;
  }

  PermeabilityTensor permeability() const
  {
    return _permeability;
  }

  CouplingParameters()
  {
    for (int i = 0; i < GV::dimension; ++i)
      for (int j = 0; j < GV::dimension; ++j)
        _permeability[i][j] = (i == j ? 1.0 : 0.0);
  }

private:
  PermeabilityTensor _permeability;

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
    : parameters(params)
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
    const LFSU_V_PFS& lfsu_v_pfs = stokeslfsu.template getChild<0>();

    typedef typename LFSU_V_PFS::template Child<0>::Type LFSU_V;
    const unsigned int vsize = lfsu_v_pfs.getChild(0).size();

    // domain and range field type
    typedef typename LFSU_V::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU_V::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RT_V;
    typedef typename LFSU_V::Traits::SizeType size_type;

    typedef typename StokesLFSU::template Child<1>::Type LFSU_P;
    //const LFSU_P& lfsu_p = stokeslfsu.template getChild<1>();
    //const unsigned int psize = lfsu_p.size();

    typedef typename LFSU_P::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU_P::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RT_P;

    typedef typename DarcyLFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RT_D;
    typedef typename DarcyLFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianType_D;
    const unsigned int dsize = darcylfsu.size();

    typedef typename IG::Geometry::LocalCoordinate LC;
    typedef typename IG::Geometry::GlobalCoordinate GC;

    // select quadrature rule
    Dune::GeometryType gt = ig.geometry().type();
    const int qorder = 2 * std::max(lfsu_v_pfs.template getChild(0).localFiniteElement().localBasis().order(),
                                    darcylfsu.localFiniteElement().localBasis().order());

    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gt,qorder);

    const typename IG::EntityPointer darcyCell = ig.outside();

    const RF g = parameters.gravity();
    const RF alpha = parameters.alpha();
    const RF nu = parameters.viscosity();
    const RF S = parameters.porosity();
    const typename Parameters::PermeabilityTensor Pi = parameters.permeability();
    RF tracePi = 0.0;
    for (int i = 0; i < dim; ++i)
      tracePi += Pi[i][i];

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {

        const GC pos = ig.geometry().global(it->position());
        const GC stokesPos = ig.geometryInInside().global(it->position());
        const GC darcyPos = ig.geometryInOutside().global(it->position());

        // integration weight
        const RF factor = it->weight() * ig.geometry().integrationElement(it->position());

        std::vector<RT_V> v(vsize);
        lfsu_v_pfs.getChild(0).localFiniteElement().localBasis().evaluateFunction(stokesPos,v);

        std::vector<RT_D> psi(dsize);
        darcylfsu.localFiniteElement().localBasis().evaluateFunction(darcyPos,psi);

        // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
        std::vector<JacobianType_D> js(dsize);
        darcylfsu.localFiniteElement().localBasis().evaluateJacobian(darcyPos,js);

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
            phi += darcyx[darcylfsu.localIndex(i)] * psi[i];
            gradphi.axpy(darcyx[darcylfsu.localIndex(i)],gradpsi[i]);
          }

        Dune::FieldVector<RF,dim> u(0.0);
        const GC n = ig.unitOuterNormal(it->position());

        // calculate u
        for (int d = 0; d < dim; ++d)
          {
            const LFSU_V& lfsu_v = lfsu_v_pfs.getChild(d);
            // calculate d-th component of u
            for (size_type i = 0; i < lfsu_v.size(); ++i)
              u[d] += stokesx[lfsu_v.localIndex(i)] * v[i];
          }

        for (size_type i = 0; i < darcylfsu.size(); ++i)
          darcyr[darcylfsu.localIndex(i)] -= 1.0/S * (u * n) * psi[i] * factor;

        Dune::FieldVector<RF,dim> tangentialFlow(0.0);
        Pi.mv(gradphi,tangentialFlow);
        tangentialFlow *= g/nu;
        tangentialFlow += u;
        // project into tangential plane
        GC scaledNormal = n;
        scaledNormal *= (tangentialFlow * n);
        tangentialFlow -= scaledNormal;

        for (int d = 0; d < dim; ++d)
          {
            const LFSU_V& lfsu_v = lfsu_v_pfs.getChild(d);
            for (size_type i = 0; i < lfsu_v.size(); ++i)
                stokesr[lfsu_v.localIndex(i)] += g * (phi - pos[dim]) * v[i] * n[d] * factor;

            for (size_type i = 0; i < lfsu_v.size(); ++i)
                stokesr[lfsu_v.localIndex(i)] += nu * alpha * sqrt(dim) / tracePi * tangentialFlow[d] * v[i] * factor;
          }

      }
  }

private:

  const Parameters& parameters;

};


#endif // STOKESDARCYCOUPLINGOPERATOR_HH
