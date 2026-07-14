#ifndef __END_EcoMugFlux__
#define __END_EcoMugFlux__

// Sea-level cosmic-muon flux and Metropolis-Hastings sampler (ECoMUG,
// Pagano et al. NIM A 1014 (2021) 165732). Pure C++, no Geant4 — so the
// physics can be unit-tested and the integrated rate recomputed against the
// reference Fortran (ecomug.f90) without a RAT build.
//
// Units: momentum p [GeV/c], x = cos(zenith) in [0,1], flux in
// muons m^-2 sr^-1 s^-1 (GeV/c)^-1.

#include <algorithm>
#include <cmath>

namespace END {
namespace ecomug {

constexpr double PI          = 3.141592653589793;
constexpr double MU_MASS_GEV = 0.1056584;

// Generation window; the MH chain samples this truncated support, so the
// integrated rate must use the same limits.
constexpr double PMU_LO = 1.0, PMU_HI = 2000.0;   // GeV/c
constexpr double CTH_LO = 0.0, CTH_HI = 1.0;

// Metropolis proposal widths and thinning (see ecomug.f90).
constexpr double DLGP_STEP = 0.15;    // sigma of log10(p) proposal
constexpr double DCT_STEP  = 0.10;    // sigma of cos(theta) proposal
constexpr int    NTHIN     = 10;      // MH steps per recorded sample

// ECoMUG sea-level differential flux at the surface (depth 0).
inline double dalt_flux(double p, double x) {
  double npo = std::max(0.1, 2.856 - 0.655 * std::log(p));
  return 1600.0 * std::pow(p + 2.68, -3.175) * std::pow(p, 0.279) * std::pow(x, npo);
}

// Horizontal-surface PDF kernel: flux projected onto a flat plane (flux * cos).
// This is also the integrand of the plane-crossing rate.
inline double horpdf(double p, double x) {
  if (p < PMU_LO || p > PMU_HI || x < CTH_LO || x > CTH_HI) return 0.0;
  return dalt_flux(p, x) * x;
}

// Muon rate crossing a horizontal plane, per unit area [Hz / m^2]:
//   R/A = 2*pi * integral_0^1 dx  integral_plo^phi dp  dalt_flux(p,x) * x
// Trapezoidal, log-spaced in p (the flux spans three decades and falls steeply).
inline double integrated_area_flux_hz() {
  constexpr int NP = 2000, NX = 1000;
  const double u_lo = std::log(PMU_LO), u_hi = std::log(PMU_HI);
  const double du = (u_hi - u_lo) / NP;
  const double dx = (CTH_HI - CTH_LO) / NX;

  double sum = 0.0;
  for (int ix = 0; ix <= NX; ++ix) {
    double x  = CTH_LO + ix * dx;
    double wx = (ix == 0 || ix == NX) ? 0.5 : 1.0;
    double inner = 0.0;
    for (int ip = 0; ip <= NP; ++ip) {
      double p  = std::exp(u_lo + ip * du);
      double wp = (ip == 0 || ip == NP) ? 0.5 : 1.0;
      inner += wp * horpdf(p, x) * p;      // dp = p du  (log spacing)
    }
    sum += wx * inner * du;
  }
  return 2.0 * PI * sum * dx;
}

// One Metropolis-Hastings step for the horizontal-surface flux.
// Proposal: p *= 10^(N(0,1)*DLGP_STEP) (log-space), x += N(0,1)*DCT_STEP.
// The log-space p proposal is symmetric in log10(p) but the target is in p,
// so detailed balance needs the Jacobian factor new_p/cur_p in the ratio.
// gauss() -> N(0,1), uniform() -> U(0,1).
template <class FG, class FU>
inline void mh_single(double &cur_p, double &cur_x, double &cur_f, FG &&gauss, FU &&uniform) {
  double new_p = cur_p * std::pow(10.0, gauss() * DLGP_STEP);
  double new_x = cur_x + gauss() * DCT_STEP;
  double new_f = horpdf(new_p, new_x);
  if (new_f <= 0.0) return;
  if (cur_f == 0.0) {
    cur_p = new_p;  cur_x = new_x;  cur_f = new_f;
  } else {
    double ratio = (new_f * new_p) / (cur_f * cur_p);   // includes Jacobian new_p/cur_p
    if (uniform() < std::min(1.0, ratio)) { cur_p = new_p;  cur_x = new_x;  cur_f = new_f; }
  }
}

template <class FG, class FU>
inline void mh_thinned(double &cur_p, double &cur_x, double &cur_f, FG &&gauss, FU &&uniform) {
  for (int i = 0; i < NTHIN; ++i) mh_single(cur_p, cur_x, cur_f, gauss, uniform);
}

}  // namespace ecomug
}  // namespace END

#endif
