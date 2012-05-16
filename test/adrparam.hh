// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PM_ADVECTIONDIFFUSIONREACTIONPARAMETERS_HH
#define DUNE_PM_ADVECTIONDIFFUSIONREACTIONPARAMETERS_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/geometry/type.hh>

namespace Dune {
  namespace PM {
    //! \addtogroup LocalOperator
    //! \ingroup PM
    //! \{

    /** Parameter class for the general advection-diffusion-reaction equation in the form
     *
     * \f{align*}{
     *   \nabla\cdot\{ b u - K(x) w(x,u) \nabla u \} + a u &=& f \mbox{ in } \Omega,  \\
     *                                                   u &=& g \mbox{ on } \Gamma_D \\
     *                  (b u - K(x)w(x,u)\nabla u) \cdot n &=& j \mbox{ on } \Gamma_F \\
     *                       (K(x)(w(x,u)\nabla u) \cdot n &=& 0 \mbox{ on } \Gamma_O
      * \f}
     */

	//! traits class for two phase parameter class
	template<typename GV, typename RF>
	struct ADRTraits
	{
	  //! \brief the grid view
	  typedef GV GridViewType;

	  //! \brief Enum for domain dimension
	  enum {
		//! \brief dimension of the domain
		dimDomain = GV::dimension
	  };

	  //! \brief Export type for domain field
	  typedef typename GV::Grid::ctype DomainFieldType;

	  //! \brief domain type
	  typedef Dune::FieldVector<DomainFieldType,dimDomain> DomainType;

	  //! \brief domain type
	  typedef Dune::FieldVector<DomainFieldType,dimDomain-1> IntersectionDomainType;

	  //! \brief Export type for range field
	  typedef RF RangeFieldType;

	  //! \brief range type
	  typedef Dune::FieldVector<RF,GV::dimensionworld> RangeType;

      //! \brief permeability tensor type
      typedef Dune::FieldMatrix<RangeFieldType,dimDomain,dimDomain> PermTensorType;

	  //! grid types
	  typedef typename GV::Traits::template Codim<0>::Entity ElementType;
	  typedef typename GV::Intersection IntersectionType;
	};

    struct ADRBoundaryConditions
    {
      enum Type { Dirichlet, Flux, Outflow, None };
    };

    //! Definition of parameters, T is the traits class
	template<typename GV, typename RF>
	class ADRIF
	{
	public:
	  typedef ADRTraits<GV,RF> Traits;

      typedef ADRBoundaryConditions::Type BCType;

	  //! velocity field
	  typename Traits::RangeType
	  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return typename Traits::RangeType();
	  }

	  //! tensor diffusion coefficient
	  typename Traits::PermTensorType
	  K (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return typename Traits::PermTensorType();
	  }

	  //! nonlinear scalar diffusion term
	  typename Traits::RangeFieldType
	  w (const typename Traits::ElementType& e, const typename Traits::DomainType& x, typename Traits::RangeFieldType u) const
	  {
		return 1.0;
	  }

	  //! reaction term
	  typename Traits::RangeFieldType
	  a (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return 0.0;
	  }

	  //! source term
	  typename Traits::RangeFieldType
	  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return 0.0;
	  }

	  //! boundary condition type function
	  BCType
	  bc (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
	  {
		return ADRBoundaryConditions::Dirichlet;
	  }

      //! Dirichlet boundary condition value
      typename Traits::RangeFieldType
      g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        return 0;
      }

      //! Neumann boundary condition
      typename Traits::RangeFieldType
      j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        return 0.0;
      }

      //! set time for subsequent evaluation
      void setTime (typename Traits::RangeFieldType t)
      {
      }

	};

   //! \} group PM
  } // namespace PM
} // namespace Dune

#endif
