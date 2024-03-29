#ifndef REVTPSSC_EPS_H
#define REVTPSSC_EPS_H

#include "functional.h"
#include "constants.h"
#include "pbec_eps.h"

namespace revtpssc_eps
{
  using namespace pbec_eps;

  template<class num, class T>
  static num revtpssA(const num &eps, const T &u3, const num &beta_tpss)
  {
    using xc_constants::param_gamma;
    num beta_gamma = beta_tpss/param_gamma;
    return beta_gamma/expm1(-eps/(param_gamma*u3));
  }

  template<class num, class T>
  static num revtpssH(const num &d2, const num &eps, const T &u3, const num &beta_tpss)
  {
    num d2A = d2*revtpssA(eps,u3,beta_tpss);
    using xc_constants::param_gamma;
    num beta_gamma = beta_tpss/param_gamma;
    return param_gamma*u3*
      log(1+beta_gamma*d2*(1 + d2A)/(1+d2A*(1+d2A)));
  }

  template<class num>
  static num revtpss_beta(const num &dens)
  {
    num r_s = cbrt(3/(4*PI*dens));
    using xc_constants::param_beta_pbe_paper;
    return param_beta_pbe_paper*(1+0.1*r_s)/(1+0.1778*r_s);
  }

  template<class num>
  static num revtpss_pbec_eps(const densvars<num> &d)
  {
    num beta_tpss = revtpss_beta(d.n);
    num eps = pw92eps::pw92eps(d);
    num u = phi(d);
    // Avoiding the square root of d.gnn here
    num d2 = pow(1.0/12*pow(3,5.0/6.0)/pow(M_PI,-1.0/6),2)*
      d.gnn/(u*u*pow(d.n,7.0/3.0));
    return (eps + revtpssH(d2,eps,pow3(u),beta_tpss));
  }

  template<class num>
  static num revtpss_pbec_eps_polarized(const num &a, const num &gaa)
  {
    num eps = pw92eps::pw92eps_polarized(a);
    parameter u = pow(2.0,-1.0/3.0); //phi(d) for alpha or beta density =0
    num beta_tpss = revtpss_beta(a);
    // Avoiding the square root of d.gnn here
    num d2 = pow(1.0/12*pow(3,5.0/6.0)/pow(M_PI,-1.0/6),2)*
      gaa/(u*u*pow(a,7.0/3.0));
    return (eps + revtpssH(d2,eps,pow3(u),beta_tpss));
  }


  template<class num>
  static num C(const densvars<num> &d)
  {
    num gzeta2 = ( pow2(d.n)*d.gss - 2*d.n*d.s*d.gns + pow2(d.s)*d.gnn )/pow(d.n,4); // (grad zeta)^2
    num xi2 = gzeta2 / (4*pow(3*pow2(M_PI)*d.n,2.0/3.0));  // xi^2
    num C0 = 0.59 + 0.9269*pow2(d.zeta) + 0.6225*pow(d.zeta,4) + 2.1540*pow(d.zeta,6); // C(zeta,0)
    return C0/pow( 1 + 0.5*xi2*(pow(1+d.zeta,-4.0/3.0) + pow(1-d.zeta,-4.0/3.0)),4);
  }

  template<class num>
  static num epsc_summax(const densvars<num> &d) /* sum_sigma n_sigma*epsc_max/n is the sum in eq [12] of the reference */
  {
    num epsc_pbe = revtpss_pbec_eps(d);
    num epsc_pbe_a = revtpss_pbec_eps_polarized(d.a,d.gaa);
    num epsc_pbe_b = revtpss_pbec_eps_polarized(d.b,d.gbb);
    return ( d.a*max(epsc_pbe,epsc_pbe_a) + d.b*max(epsc_pbe,epsc_pbe_b) )/d.n;
  }

  template<class num>
  static num epsc_revpkzb(const densvars<num> &d) /* eps_c^{rev PKZB} it is the eq [12] from the reference */
  {
    num tauwtau2 = pow2(d.gnn/(8.0*d.n*d.tau));
    num epsc_sum = epsc_summax(d);
    num epsc_pbe = revtpss_pbec_eps(d);
    num C_zeta_xi = C(d);
    return epsc_pbe*(1 + C_zeta_xi*tauwtau2) - (1 + C_zeta_xi)*tauwtau2*epsc_sum;
  }

  template<class num>
  static num revtpssc_eps(const densvars<num> &d)
  {
    num eps_pkzb = epsc_revpkzb(d);
    num tauwtau2 = pow2(d.gnn/(8.0*d.n*d.tau));
    parameter dd = 2.8;
    return eps_pkzb*(1 + dd*eps_pkzb*tauwtau2);
  }


}

#endif
