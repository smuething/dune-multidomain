#ifndef DUNE_MULTIDOMAIN_COUPLINGUTILITIES_HH
#define DUNE_MULTIDOMAIN_COUPLINGUTILITIES_HH

namespace Dune {

namespace PDELab {

namespace MultiDomain {


class FullCouplingPattern
{

public:

  template<typename LFSU1, typename LFSV1, typename LFSU2, typename LFSV2>
  void pattern_coupling (const LFSU1& lfsu_s, const LFSV1& lfsv_s, const LFSU2& lfsu_n, const LFSV2& lfsv_n,
                         Dune::PDELab::LocalSparsityPattern& pattern_sn,
                         Dune::PDELab::LocalSparsityPattern& pattern_ns) const
  {
    for (unsigned int i=0; i<lfsv_s.size(); ++i)
      for (unsigned int j=0; j<lfsu_n.size(); ++j)
          pattern_sn.push_back(Dune::PDELab::SparsityLink(i,j));
    for (unsigned int i=0; i<lfsv_n.size(); ++i)
      for (unsigned int j=0; j<lfsu_s.size(); ++j)
          pattern_ns.push_back(Dune::PDELab::SparsityLink(i,j));
  }

};


//! Implement jacobian_skeleton() based on alpha_skeleton()
/**
 * Derive from this class to add numerical jacobian for skeleton.  The
 * derived class needs to implement alpha_skeleton().
 *
 * \tparam Imp Type of the derived class (CRTP-trick).
 */
template<typename Imp>
class NumericalJacobianCoupling
{
public:
  NumericalJacobianCoupling ()
    : epsilon(1e-11)
  {}

  NumericalJacobianCoupling (double epsilon_)
    : epsilon(epsilon_)
  {}

  //! compute local jacobian of the skeleton term
  template<typename IG, typename LFSU_S, typename LFSU_N,
           typename X, typename LFSV_S, typename LFSV_N,
           typename R>
  void jacobian_coupling
  ( const IG& ig,
    const LFSU_S& lfsu_s, const X& x_s, const LFSV_S& lfsv_s,
    const LFSU_N& lfsu_n, const X& x_n, const LFSV_N& lfsv_n,
    LocalMatrix<R>& mat_ss, LocalMatrix<R>& mat_sn,
    LocalMatrix<R>& mat_ns, LocalMatrix<R>& mat_nn) const
  {
    const int m_s=lfsv_s.size();
    const int m_n=lfsv_n.size();
    const int n_s=lfsu_s.size();
    const int n_n=lfsu_n.size();

    X u_s(x_s);
    X u_n(x_n);
    std::vector<R> down_s(m_s,0.0),up_s(m_s);
    std::vector<R> down_n(m_n,0.0),up_n(m_n);

    // base line
    asImp().alpha_coupling(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,down_s,down_n);

    // jiggle in self
    for (int j=0; j<n_s; j++)
      {
        for (int k=0; k<m_s; k++) up_s[k]=0.0;
        for (int k=0; k<m_n; k++) up_n[k]=0.0;
        R delta = epsilon*(1.0+std::abs(u_s[j]));
        u_s[j] += delta;
        asImp().alpha_coupling(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,up_s,up_n);
        for (int i=0; i<m_s; i++)
          mat_ss(i,j) += (up_s[i]-down_s[i])/delta;
        for (int i=0; i<m_n; i++)
          mat_ns(i,j) += (up_n[i]-down_n[i])/delta;
        u_s[j] = x_s[j];
      }

    // jiggle in neighbor
    for (int j=0; j<n_n; j++)
      {
        for (int k=0; k<m_s; k++) up_s[k]=0.0;
        for (int k=0; k<m_n; k++) up_n[k]=0.0;
        R delta = epsilon*(1.0+std::abs(u_n[j]));
        u_n[j] += delta;
        asImp().alpha_coupling(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,up_s,up_n);
        for (int i=0; i<m_s; i++)
          mat_sn(i,j) += (up_s[i]-down_s[i])/delta;
        for (int i=0; i<m_n; i++)
          mat_nn(i,j) += (up_n[i]-down_n[i])/delta;
        u_n[j] = x_n[j];
      }
  }

private:
  const double epsilon; // problem: this depends on data type R!
  Imp& asImp () { return static_cast<Imp &> (*this); }
  const Imp& asImp () const { return static_cast<const Imp &>(*this); }
};


//! Implement jacobian_apply_skeleton() based on alpha_skeleton()
/**
 * Derive from this class to add numerical jacobian application for
 * skeleton.  The derived class needs to implement alpha_skeleton().
 *
 * \tparam Imp Type of the derived class (CRTP-trick).
 */
template<typename Imp>
class NumericalJacobianApplyCoupling
{
public:
  NumericalJacobianApplyCoupling ()
    : epsilon(1e-11)
  {}

  NumericalJacobianApplyCoupling (double epsilon_)
    : epsilon(epsilon_)
  {}

  //! apply local jacobian of the skeleton term
  template<typename IG, typename LFSU_S, typename LFSU_N,
           typename X, typename LFSV_S, typename LFSV_N,
           typename Y>
  void jacobian_apply_coupling
  ( const IG& ig,
    const LFSU_S& lfsu_s, const X& x_s, const LFSV_S& lfsv_s,
    const LFSU_S& lfsu_n, const X& x_n, const LFSV_S& lfsv_n,
    Y& y_s, Y& y_n) const
  {
    typedef typename X::value_type R;
    const int m_s=lfsv_s.size();
    const int m_n=lfsv_n.size();
    const int n_s=lfsu_s.size();
    const int n_n=lfsu_n.size();

    X u_s(x_s);
    X u_n(x_n);
    std::vector<R> down_s(m_s,0.0),up_s(m_s);
    std::vector<R> down_n(m_n,0.0),up_n(m_n);

    // base line
    asImp().alpha_coupling(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,down_s,down_n);

    // jiggle in self
    for (int j=0; j<n_s; j++)
      {
        for (int k=0; k<m_s; k++) up_s[k]=0.0;
        for (int k=0; k<m_n; k++) up_n[k]=0.0;
        R delta = epsilon*(1.0+std::abs(u_s[j]));
        u_s[j] += delta;
        asImp().alpha_coupling(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,up_s,up_n);
        for (int i=0; i<m_s; i++)
          y_s[i] += ((up_s[i]-down_s[i])/delta)*x_s[j];
        for (int i=0; i<m_n; i++)
          y_n[i] += ((up_n[i]-down_n[i])/delta)*x_s[j];
        u_s[j] = x_s[j];
      }

    // jiggle in neighbor
    for (int j=0; j<n_n; j++)
      {
        for (int k=0; k<m_s; k++) up_s[k]=0.0;
        for (int k=0; k<m_n; k++) up_n[k]=0.0;
        R delta = epsilon*(1.0+std::abs(u_n[j]));
        u_n[j] += delta;
        asImp().alpha_coupling(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,up_s,up_n);
        for (int i=0; i<m_s; i++)
          y_s[i] += ((up_s[i]-down_s[i])/delta)*x_n[j];
        for (int i=0; i<m_n; i++)
          y_n[i] += ((up_n[i]-down_n[i])/delta)*x_n[j];
        u_n[j] = x_n[j];
      }
  }

private:
  const double epsilon; // problem: this depends on data type R!
  Imp& asImp () { return static_cast<Imp &> (*this); }
  const Imp& asImp () const { return static_cast<const Imp &>(*this); }
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_COUPLINGUTILITIES_HH
