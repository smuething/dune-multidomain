#ifndef DUNE_MULTIDOMAIN_COUPLINGUTILITIES_HH
#define DUNE_MULTIDOMAIN_COUPLINGUTILITIES_HH

namespace Dune {

namespace PDELab {

namespace MultiDomain {


class CouplingOperatorDefaultFlags
{
public:
  static const bool doSkeletonTwoSided = false;
  static const bool doPatternCoupling = false;
  static const bool doPatternEnrichedCouplingToSubProblems = false;
  static const bool doPatternEnrichedCoupling = false;
  static const bool doAlphaCoupling = false;
  static const bool doAlphaEnrichedCouplingToSubProblems = false;
  static const bool doAlphaEnrichedCoupling = false;
  static const bool doLambdaCoupling = false;
  static const bool doLambdaEnrichedCouplingToSubProblems = false;
  static const bool doLambdaEnrichedCoupling = false;
};


class FullCouplingPattern
{

public:

  template<typename LFSU1, typename LFSV1, typename LFSU2, typename LFSV2>
  void pattern_coupling (const LFSU1& lfsu_s, const LFSV1& lfsv_s, const LFSU2& lfsu_n, const LFSV2& lfsv_n,
                         Dune::PDELab::LocalSparsityPattern& pattern_ss,
                         Dune::PDELab::LocalSparsityPattern& pattern_sn,
                         Dune::PDELab::LocalSparsityPattern& pattern_ns,
                         Dune::PDELab::LocalSparsityPattern& pattern_nn) const
  {
    for (unsigned int i=0; i<lfsv_s.size(); ++i)
      for (unsigned int j=0; j<lfsu_n.size(); ++j)
        pattern_sn.push_back(Dune::PDELab::SparsityLink(lfsv_s.localIndex(i),lfsu_n.localIndex(j)));
    for (unsigned int i=0; i<lfsv_n.size(); ++i)
      for (unsigned int j=0; j<lfsu_s.size(); ++j)
        pattern_ns.push_back(Dune::PDELab::SparsityLink(lfsv_n.localIndex(i),lfsu_s.localIndex(j)));
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
    : epsilon(1e-8)
  {}

  NumericalJacobianCoupling (double epsilon_)
    : epsilon(epsilon_)
  {}

  //! compute local jacobian of the skeleton term
  template<typename IG, typename LFSU_S, typename LFSU_N,
           typename X, typename LFSV_S, typename LFSV_N,
           typename Jacobian>
  void jacobian_coupling
  ( const IG& ig,
    const LFSU_S& lfsu_s, const X& x_s, const LFSV_S& lfsv_s,
    const LFSU_N& lfsu_n, const X& x_n, const LFSV_N& lfsv_n,
    Jacobian& mat_ss, Jacobian& mat_sn,
    Jacobian& mat_ns, Jacobian& mat_nn) const
  {
    typedef typename X::value_type D;
    typedef typename Jacobian::value_type R;
    typedef LocalVector<R,TestSpaceTag,typename Jacobian::weight_type> ResidualVector;
    typedef typename ResidualVector::WeightedAccumulationView ResidualView;

    const int m_s=lfsv_s.size();
    const int m_n=lfsv_n.size();
    const int n_s=lfsu_s.size();
    const int n_n=lfsu_n.size();

    X u_s(x_s);
    X u_n(x_n);

    ResidualVector down_s(mat_ss.nrows()),up_s(mat_ss.nrows());
    ResidualView downview_s = down_s.weightedAccumulationView(1.0);
    ResidualView upview_s = up_s.weightedAccumulationView(1.0);

    ResidualVector down_n(mat_nn.nrows()),up_n(mat_nn.nrows());
    ResidualView downview_n = down_n.weightedAccumulationView(1.0);
    ResidualView upview_n = up_n.weightedAccumulationView(1.0);

    // base line
    asImp().alpha_coupling(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,downview_s,downview_n);

    // jiggle in self
    for (int j=0; j<n_s; j++)
      {
        up_s = 0.0;
        up_n = 0.0;
        D delta = epsilon*(1.0+std::abs(u_s(lfsu_s,j)));
        u_s(lfsu_s,j) += delta;
        asImp().alpha_coupling(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,upview_s,upview_n);
        for (int i=0; i<m_s; i++)
          mat_ss.accumulate(lfsv_s,i,lfsu_s,j,(up_s(lfsv_s,i)-down_s(lfsv_s,i))/delta);
        for (int i=0; i<m_n; i++)
          mat_ns.accumulate(lfsv_n,i,lfsu_s,j,(up_n(lfsv_n,i)-down_n(lfsv_n,i))/delta);
        u_s(lfsu_s,j) = x_s(lfsu_s,j);
      }

    // jiggle in neighbor
    for (int j=0; j<n_n; j++)
      {
        up_s = 0.0;
        up_n = 0.0;
        D delta = epsilon*(1.0+std::abs(u_n(lfsu_n,j)));
        u_n(lfsu_n,j) += delta;
        asImp().alpha_coupling(ig,lfsu_s,u_s,lfsv_s,lfsu_n,u_n,lfsv_n,upview_s,upview_n);
        for (int i=0; i<m_s; i++)
          mat_sn.accumulate(lfsv_s,i,lfsu_n,j,(up_s(lfsv_s,i)-down_s(lfsv_s,i))/delta);
        for (int i=0; i<m_n; i++)
          mat_nn.accumulate(lfsv_n,i,lfsu_n,j,(up_n(lfsv_n,i)-down_n(lfsv_n,i))/delta);
        u_n(lfsu_n,j) = x_n(lfsu_n,j);
      }
  }

private:
  const double epsilon; // problem: this depends on data type R!
  Imp& asImp () { return static_cast<Imp &> (*this); }
  const Imp& asImp () const { return static_cast<const Imp &>(*this); }
};


//! FIXME: This still needs to be ported!!!
template<typename Imp>
class NumericalJacobianApplyCoupling
{
public:
  NumericalJacobianApplyCoupling ()
    : epsilon(1e-8)
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


class FullEnrichedCouplingFirstSubProblemPattern
{

public:

  template<typename LFSU1, typename LFSV1,
           typename CouplingLFSU, typename CouplingLFSV
           >
  void pattern_enriched_coupling_first (const LFSU1& lfsu_s, const LFSV1& lfsv_s,
                                        const CouplingLFSU& coupling_lfsu, const CouplingLFSV& coupling_lfsv,
                                        Dune::PDELab::LocalSparsityPattern& pattern_ss,
                                        Dune::PDELab::LocalSparsityPattern& pattern_sc,
                                        Dune::PDELab::LocalSparsityPattern& pattern_cs,
                                        Dune::PDELab::LocalSparsityPattern& pattern_cc) const
  {
    for (unsigned int i=0; i<lfsv_s.size(); ++i)
      for (unsigned int j=0; j<coupling_lfsu.size(); ++j)
        pattern_sc.push_back(Dune::PDELab::SparsityLink(lfsv_s.localIndex(i),coupling_lfsu.localIndex(j)));
    for (unsigned int i=0; i<coupling_lfsv.size(); ++i)
      for (unsigned int j=0; j<lfsu_s.size(); ++j)
        pattern_cs.push_back(Dune::PDELab::SparsityLink(coupling_lfsv.localIndex(i),lfsu_s.localIndex(j)));
  }

};


class FullEnrichedCouplingSecondSubProblemPattern
{

public:

  template<typename LFSU2, typename LFSV2,
           typename CouplingLFSU, typename CouplingLFSV
           >
  void pattern_enriched_coupling_second (const LFSU2& lfsu_n, const LFSV2& lfsv_n,
                                         const CouplingLFSU& coupling_lfsu, const CouplingLFSV& coupling_lfsv,
                                         Dune::PDELab::LocalSparsityPattern& pattern_nn,
                                         Dune::PDELab::LocalSparsityPattern& pattern_nc,
                                         Dune::PDELab::LocalSparsityPattern& pattern_cn,
                                         Dune::PDELab::LocalSparsityPattern& pattern_cc) const
  {
    for (unsigned int i=0; i<lfsv_n.size(); ++i)
      for (unsigned int j=0; j<coupling_lfsu.size(); ++j)
        pattern_nc.push_back(Dune::PDELab::SparsityLink(lfsv_n.localIndex(i),coupling_lfsu.localIndex(j)));
    for (unsigned int i=0; i<coupling_lfsv.size(); ++i)
      for (unsigned int j=0; j<lfsu_n.size(); ++j)
        pattern_cn.push_back(Dune::PDELab::SparsityLink(coupling_lfsv.localIndex(i),lfsu_n.localIndex(j)));
  }

};


class FullEnrichedCouplingPattern
{

public:

  template<
           typename CouplingLFSU, typename CouplingLFSV
           >
  void pattern_enriched_coupling (const CouplingLFSU& coupling_lfsu, const CouplingLFSV& coupling_lfsv,
                                  Dune::PDELab::LocalSparsityPattern& pattern_cc
                                  ) const
  {
    for (unsigned int i=0; i<coupling_lfsv.size(); ++i)
      for (unsigned int j=0; j<coupling_lfsu.size(); ++j)
        pattern_cc.push_back(Dune::PDELab::SparsityLink(coupling_lfsv.localIndex(i),coupling_lfsu.localIndex(j)));
  }

};


//! Implement jacobian_enriched_coupling_first() based on numerical differentiation of the residuals.
template<typename Imp>
class NumericalJacobianEnrichedCouplingToFirstSubProblem
{
public:
  NumericalJacobianEnrichedCouplingToFirstSubProblem ()
    : epsilon(1e-8)
  {}

  NumericalJacobianEnrichedCouplingToFirstSubProblem (double epsilon_)
    : epsilon(epsilon_)
  {}

  //! compute local jacobian of the coupling term between mortar space and first subproblem
  template<typename IG, typename X,
           typename LFSU_S, typename LFSV_S,
           typename LFSU_C, typename LFSV_C,
           typename Jacobian>
  void jacobian_enriched_coupling_first
  ( const IG& ig,
    const LFSU_S& lfsu_s, const X& x_s, const LFSV_S& lfsv_s,
    const LFSU_C& lfsu_c, const X& x_c, const LFSV_C& lfsv_c,
    Jacobian& mat_ss, Jacobian& mat_sc,
    Jacobian& mat_cs, Jacobian& mat_cc) const
  {
    typedef typename X::value_type D;
    typedef typename Jacobian::value_type R;
    typedef LocalVector<R,TestSpaceTag,typename Jacobian::weight_type> ResidualVector;
    typedef typename ResidualVector::WeightedAccumulationView ResidualView;

    const int m_s=lfsv_s.size();
    const int m_c=lfsv_c.size();
    const int n_s=lfsu_s.size();
    const int n_c=lfsu_c.size();

    X u_s(x_s);
    X u_c(x_c);

    ResidualVector down_s(mat_ss.nrows()),up_s(mat_ss.nrows());
    ResidualView downview_s = down_s.weightedAccumulationView(1.0);
    ResidualView upview_s = up_s.weightedAccumulationView(1.0);

    ResidualVector down_c(mat_cc.nrows()),up_c(mat_cc.nrows());
    ResidualView downview_c = down_c.weightedAccumulationView(1.0);
    ResidualView upview_c = up_c.weightedAccumulationView(1.0);

    // base line
    asImp().alpha_enriched_coupling_first(ig,
                                          lfsu_s,u_s,lfsv_s,
                                          lfsu_c,u_c,lfsv_c,
                                          downview_s,downview_c);

    // jiggle in self
    for (int j=0; j<n_s; ++j)
      {
        up_s = 0.0;
        up_c = 0.0;

        D delta = epsilon*(1.0+std::abs(u_s(lfsu_s,j)));
        u_s(lfsu_s,j) += delta;
        asImp().alpha_enriched_coupling_first(ig,
                                              lfsu_s,u_s,lfsv_s,
                                              lfsu_c,u_c,lfsv_c,
                                              upview_s,upview_c);
        for (int i=0; i<m_s; ++i)
          mat_ss.accumulate(lfsv_s,i,lfsu_s,j,(up_s(lfsv_s,i)-down_s(lfsv_s,i))/delta);
        for (int i=0; i<m_c; ++i)
          mat_cs.accumulate(lfsv_c,i,lfsu_s,j,(up_c(lfsv_c,i)-down_c(lfsv_c,i))/delta);
        u_s(lfsu_s,j) = x_s(lfsu_s,j);
      }

    // jiggle in coupling
    for (int j=0; j<n_c; j++)
      {
        up_s = 0.0;
        up_c = 0.0;

        D delta = epsilon*(1.0+std::abs(u_c(lfsu_c,j)));
        u_c(lfsu_c,j) += delta;
        asImp().alpha_enriched_coupling_first(ig,
                                              lfsu_s,u_s,lfsv_s,
                                              lfsu_c,u_c,lfsv_c,
                                              upview_s,upview_c);
        for (int i=0; i<m_s; ++i)
          mat_sc.accumulate(lfsv_s,i,lfsu_c,j,(up_s(lfsv_s,i)-down_s(lfsv_s,i))/delta);
        for (int i=0; i<m_c; ++i)
          mat_cc.accumulate(lfsv_c,i,lfsu_c,j,(up_c(lfsv_c,i)-down_c(lfsv_c,i))/delta);
        u_c(lfsu_c,j) = x_c(lfsu_c,j);
      }
  }

private:
  const double epsilon; // problem: this depends on data type R!
  Imp& asImp () { return static_cast<Imp &> (*this); }
  const Imp& asImp () const { return static_cast<const Imp &>(*this); }
};


//! Implement jacobian_enriched_coupling_second() based on numerical differentiation of the residuals.
template<typename Imp>
class NumericalJacobianEnrichedCouplingToSecondSubProblem
{
public:
  NumericalJacobianEnrichedCouplingToSecondSubProblem ()
    : epsilon(1e-8)
  {}

  NumericalJacobianEnrichedCouplingToSecondSubProblem (double epsilon_)
    : epsilon(epsilon_)
  {}

  //! compute local jacobian of the coupling term between the mortar space and the second subproblem
  template<typename IG, typename X,
           typename LFSU_N, typename LFSV_N,
           typename LFSU_C, typename LFSV_C,
           typename Jacobian>
  void jacobian_enriched_coupling_second
  ( const IG& ig,
    const LFSU_N& lfsu_n, const X& x_n, const LFSV_N& lfsv_n,
    const LFSU_C& lfsu_c, const X& x_c, const LFSV_C& lfsv_c,
    Jacobian& mat_nn, Jacobian& mat_nc,
    Jacobian& mat_cn, Jacobian& mat_cc) const
  {
    typedef typename X::value_type D;
    typedef typename Jacobian::value_type R;
    typedef LocalVector<R,TestSpaceTag,typename Jacobian::weight_type> ResidualVector;
    typedef typename ResidualVector::WeightedAccumulationView ResidualView;

    const int m_n=lfsv_n.size();
    const int m_c=lfsv_c.size();
    const int n_n=lfsu_n.size();
    const int n_c=lfsu_c.size();

    X u_n(x_n);
    X u_c(x_c);

    ResidualVector down_n(mat_nn.nrows()),up_n(mat_nn.nrows());
    ResidualView downview_n = down_n.weightedAccumulationView(1.0);
    ResidualView upview_n = up_n.weightedAccumulationView(1.0);

    ResidualVector down_c(mat_cc.nrows()),up_c(mat_cc.nrows());
    ResidualView downview_c = down_c.weightedAccumulationView(1.0);
    ResidualView upview_c = up_c.weightedAccumulationView(1.0);

    // base line
    asImp().alpha_enriched_coupling_second(ig,
                                           lfsu_n,u_n,lfsv_n,
                                           lfsu_c,u_c,lfsv_c,
                                           downview_n,downview_c);

    // jiggle in neighbor
    for (int j=0; j<n_n; ++j)
      {
        up_n = 0.0;
        up_c = 0.0;

        D delta = epsilon*(1.0+std::abs(u_n(lfsu_n,j)));
        u_n(lfsu_n,j) += delta;
        asImp().alpha_enriched_coupling_second(ig,
                                               lfsu_n,u_n,lfsv_n,
                                               lfsu_c,u_c,lfsv_c,
                                               upview_n,upview_c);
        for (int i=0; i<m_n; ++i)
          mat_nn.accumulate(lfsv_n,i,lfsu_n,j,(up_n(lfsv_n,i)-down_n(lfsv_n,i))/delta);
        for (int i=0; i<m_c; ++i)
          mat_cn.accumulate(lfsv_c,i,lfsu_n,j,(up_c(lfsv_c,i)-down_c(lfsv_c,i))/delta);
        u_n(lfsu_n,j) = x_n(lfsu_n,j);
      }

    // jiggle in coupling
    for (int j=0; j<n_c; j++)
      {
        up_n = 0.0;
        up_c = 0.0;

        D delta = epsilon*(1.0+std::abs(u_c(lfsu_c,j)));
        u_c(lfsu_c,j) += delta;
        asImp().alpha_enriched_coupling_second(ig,
                                               lfsu_n,u_n,lfsv_n,
                                               lfsu_c,u_c,lfsv_c,
                                               upview_n,upview_c);
        for (int i=0; i<m_n; ++i)
          mat_nc.accumulate(lfsv_n,i,lfsu_c,j,(up_n(lfsv_n,i)-down_n(lfsv_n,i))/delta);
        for (int i=0; i<m_c; ++i)
          mat_cc.accumulate(lfsv_c,i,lfsu_c,j,(up_c(lfsv_c,i)-down_c(lfsv_c,i))/delta);
        u_c(lfsu_c,j) = x_c(lfsu_c,j);
      }
  }

private:
  const double epsilon; // problem: this depends on data type R!
  Imp& asImp () { return static_cast<Imp &> (*this); }
  const Imp& asImp () const { return static_cast<const Imp &>(*this); }
};



template<typename Imp>
class NumericalJacobianEnrichedCoupling
{
public:
  NumericalJacobianEnrichedCoupling ()
    : epsilon(1e-8)
  {}

  NumericalJacobianEnrichedCoupling (double epsilon_)
    : epsilon(epsilon_)
  {}

  //! compute local jacobian of the enriched coupling term
  template<typename IG, typename X,
           typename LFSU_C, typename LFSV_C,
           typename Jacobian>
  void jacobian_enriched_coupling
  ( const IG& ig,
    const LFSU_C& lfsu_c, const X& x_c, const LFSV_C& lfsv_c,
    Jacobian& mat_cc) const
  {
    typedef typename X::value_type D;
    typedef typename Jacobian::value_type R;
    typedef LocalVector<R,TestSpaceTag,typename Jacobian::weight_type> ResidualVector;
    typedef typename ResidualVector::WeightedAccumulationView ResidualView;

    const int m_c=lfsv_c.size();
    const int n_c=lfsu_c.size();

    X u_c(x_c);

    ResidualVector down_c(mat_cc.nrows()),up_c(mat_cc.nrows());
    ResidualView downview_c = down_c.weightedAccumulationView(1.0);
    ResidualView upview_c = up_c.weightedAccumulationView(1.0);

    // base line
    asImp().alpha_enriched_coupling(ig,
                                    lfsu_c,u_c,lfsv_c,
                                    downview_c);

    // jiggle in coupling
    for (int j=0; j<n_c; j++)
      {
        up_c = 0.0;

        D delta = epsilon*(1.0+std::abs(u_c(lfsu_c,j)));
        u_c(lfsu_c,j) += delta;
        asImp().alpha_enriched_coupling(ig,
                                        lfsu_c,u_c,lfsv_c,
                                        upview_c);
        for (int i=0; i<m_c; ++i)
          mat_cc.accumulate(lfsv_c,i,lfsu_c,j,(up_c(lfsv_c,i)-down_c(lfsv_c,i))/delta);
        u_c(lfsu_c,j) = x_c(lfsu_c,j);
      }
  }

private:
  const double epsilon; // problem: this depends on data type R!
  Imp& asImp () { return static_cast<Imp &> (*this); }
  const Imp& asImp () const { return static_cast<const Imp &>(*this); }

};


//! FIXME: This still needs to be ported!!!
template<typename Imp>
class NumericalJacobianApplyEnrichedCoupling
{
public:
  NumericalJacobianApplyEnrichedCoupling ()
    : epsilon(1e-11)
  {}

  NumericalJacobianApplyEnrichedCoupling (double epsilon_)
    : epsilon(epsilon_)
  {}

  //! apply local jacobian of the skeleton term
  template<typename IG, typename LFSU_S, typename LFSU_N,
           typename X, typename LFSV_S, typename LFSV_N,
           typename LFSU_C, typename LFSV_C,
           typename Y>
  void jacobian_apply_coupling
  ( const IG& ig,
    const LFSU_S& lfsu_s, const X& x_s, const LFSV_S& lfsv_s,
    const LFSU_S& lfsu_n, const X& x_n, const LFSV_S& lfsv_n,
    const LFSU_C& lfsu_c, const X& x_c, const LFSV_C& lfsv_c,
    Y& y_s, Y& y_n, Y& y_c) const
  {
    typedef typename X::value_type R;
    const int m_s=lfsv_s.size();
    const int m_n=lfsv_n.size();
    const int m_c=lfsv_c.size();
    const int n_s=lfsu_s.size();
    const int n_n=lfsu_n.size();
    const int n_c=lfsu_c.size();

    X u_s(x_s);
    X u_n(x_n);
    X u_c(x_c);
    std::vector<R> down_s(m_s,0.0),up_s(m_s);
    std::vector<R> down_n(m_n,0.0),up_n(m_n);
    std::vector<R> down_c(m_c,0.0),up_c(m_c);

    // base line
    asImp().alpha_enriched_coupling(ig,
                                    lfsu_s,u_s,lfsv_s,
                                    lfsu_n,u_n,lfsv_n,
                                    lfsu_c,u_c,lfsv_c,
                                    down_s,down_n,down_c);

    // jiggle in self
    for (int j=0; j<n_s; j++)
      {
        std::fill(up_s.begin(),up_s.end(),0.0);
        std::fill(up_n.begin(),up_n.end(),0.0);
        std::fill(up_c.begin(),up_c.end(),0.0);

        R delta = epsilon*(1.0+std::abs(u_s[j]));
        u_s[j] += delta;
        asImp().alpha_enriched_coupling(ig,
                                        lfsu_s,u_s,lfsv_s,
                                        lfsu_n,u_n,lfsv_n,
                                        lfsu_c,u_c,lfsv_c,
                                        up_s,up_n,up_c);
        for (int i=0; i<m_s; i++)
          y_s[i] += ((up_s[i]-down_s[i])/delta)*x_s[j];
        for (int i=0; i<m_n; i++)
          y_n[i] += ((up_n[i]-down_n[i])/delta)*x_s[j];
        for (int i=0; i<m_c; i++)
          y_c[i] += ((up_c[i]-down_c[i])/delta)*x_s[j];
        u_s[j] = x_s[j];
      }

    // jiggle in neighbor
    for (int j=0; j<n_n; j++)
      {
        std::fill(up_s.begin(),up_s.end(),0.0);
        std::fill(up_n.begin(),up_n.end(),0.0);
        std::fill(up_c.begin(),up_c.end(),0.0);

        R delta = epsilon*(1.0+std::abs(u_n[j]));
        u_n[j] += delta;
        asImp().alpha_enriched_coupling(ig,
                                        lfsu_s,u_s,lfsv_s,
                                        lfsu_n,u_n,lfsv_n,
                                        lfsu_c,u_c,lfsv_c,
                                        up_s,up_n,up_c);
        for (int i=0; i<m_s; i++)
          y_s[i] += ((up_s[i]-down_s[i])/delta)*x_n[j];
        for (int i=0; i<m_n; i++)
          y_n[i] += ((up_n[i]-down_n[i])/delta)*x_n[j];
        for (int i=0; i<m_c; i++)
          y_c[i] += ((up_c[i]-down_c[i])/delta)*x_n[j];
        u_n[j] = x_n[j];
      }

    // jiggle in coupling
    for (int j=0; j<n_c; j++)
      {
        std::fill(up_s.begin(),up_s.end(),0.0);
        std::fill(up_n.begin(),up_n.end(),0.0);
        std::fill(up_c.begin(),up_c.end(),0.0);

        R delta = epsilon*(1.0+std::abs(u_c[j]));
        u_c[j] += delta;
        asImp().alpha_enriched_coupling(ig,
                                        lfsu_s,u_s,lfsv_s,
                                        lfsu_n,u_n,lfsv_n,
                                        lfsu_c,u_c,lfsv_c,
                                        up_s,up_n,up_c);
        for (int i=0; i<m_s; i++)
          y_s[i] += ((up_s[i]-down_s[i])/delta)*x_c[j];
        for (int i=0; i<m_n; i++)
          y_n[i] += ((up_n[i]-down_n[i])/delta)*x_c[j];
        for (int i=0; i<m_c; i++)
          y_c[i] += ((up_c[i]-down_c[i])/delta)*x_c[j];
        u_c[j] = x_c[j];
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
