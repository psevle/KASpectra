#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "kaspectra/constants.hpp"
#include "kaspectra/io/photon_field.hpp"
#include "kaspectra/io/proton_spectrum.hpp"
#include "kaspectra/math/integrate.hpp"

namespace kaspectra::bh {

    namespace detail {

        inline double blumenthal_T(double k, double p_minus, double cos_theta_minus) {
            return std::sqrt(k * k + p_minus * p_minus - 2.0 * k * p_minus * cos_theta_minus);
        }

        inline double blumenthal_Y(double k, double p_minus, double p_plus, double E_minus, double E_plus) {
            return (2.0 / (p_minus * p_minus)) * std::log((E_plus * E_minus + p_plus * p_minus + 1.0) / k);
        }

        inline double blumenthal_y_plus(double p_plus, double E_plus) {
            return std::log((E_plus + p_plus) / (E_plus - p_plus)) / p_plus;
        }

        inline double blumenthal_delta_plus_T(double T, double p_plus) {
            return std::log((T + p_plus) / (T - p_plus));
        }

    }   // namespace detail

    // W(k, E_minus, cos_theta_minus): dimensionless bracketed factor of
    // Blumenthal 1970 eq.10 -- dsigma/(dE_- dcos(theta_-))
    //     = (alpha*Z^2*r0^2*p_-*p_+/(2k^3)) * W(k, E_minus, cos_theta_minus)
    // Returns ONLY the bracket, not the prefactor.
    // Units: k, E_minus in m_e*c^2; E_plus = k - E_minus (energy conservation).
    // Precondition (unguarded, matches pgamma::x_pm convention): k>=2,
    // E_minus in [1,k-1]. Outside that range p_minus/p_plus are NaN.
    inline double W(double k, double E_minus, double cos_theta_minus) {
        const double E_plus = k - E_minus;
        const double p_minus = std::sqrt(E_minus * E_minus - 1.0);
        const double p_plus = std::sqrt(E_plus * E_plus - 1.0);

        const double Delta_minus = E_minus - p_minus * cos_theta_minus;
        const double Dm2 = Delta_minus * Delta_minus;
        const double Dm4 = Dm2 * Dm2;

        const double T            = detail::blumenthal_T(k, p_minus, cos_theta_minus);
        const double Y            = detail::blumenthal_Y(k, p_minus, p_plus, E_minus, E_plus);
        const double y_plus       = detail::blumenthal_y_plus(p_plus, E_plus);
        const double delta_plus_T = detail::blumenthal_delta_plus_T(T, p_plus);

        const double sin2_theta_minus = 1.0 - cos_theta_minus * cos_theta_minus;
        double const pm2 = p_minus * p_minus;

        const double term1 = -4.0 * sin2_theta_minus * (2.0 * E_minus * E_minus + 1.0) / (pm2 * Dm4);
        const double term2 = (5.0 * E_minus * E_minus - 2.0 * E_plus * E_minus + 3.0) / (pm2 * Dm2);
        const double term3 = (pm2 - k * k) / (T * T * Dm2);
        const double term4 = 2.0 * E_plus / (pm2 * Delta_minus);

        const double inner1 = 2.0 * E_minus * sin2_theta_minus * (3.0 * k + pm2 * E_plus) / Dm4;
        const double inner2 = (2.0 * E_minus * E_minus * (E_minus * E_minus + E_plus * E_plus)
                                - 7.0 * E_minus * E_minus - 3.0 * E_plus * E_minus - E_plus * E_plus + 1.0) / Dm2;
        const double inner3 = k * (E_minus * E_minus - E_minus * E_plus - 1.0) / Delta_minus;
        const double term5 = (Y / (p_minus * p_plus)) * (inner1 + inner2 + inner3);

        const double delta_group = (2.0 / Dm2) - (3.0 * k / Delta_minus) - (k * (pm2 - k * k) / (T * T * Delta_minus));
        const double term6 = -(delta_plus_T / (p_plus * T)) * delta_group;

        const double term7 = -2.0 * y_plus / Delta_minus;

        return term1 + term2 + term3 +term4 + term5 + term6 + term7;
    }

    namespace detail {
        // KA2008 eq.60 (V_p->1 ultrarelativistic limit).
        inline double xi_cos_theta_minus(double gamma_p, double e_e, double E_minus, double p_minus) {
            return (gamma_p * E_minus - e_e) / (gamma_p * p_minus);
        }
    }   // namespace detail

    // KA2008 eq.60's xi=+-1 boundary solved for E_minus. Always >= 1 (AM-GM on
    // gamma_p, e_e), matching the physical E_- >= 1 constraint -- no clamp needed.
    // e_e = E_e/m_e (dimensionless, Blumenthal's m_e*c^2 unit convention).
    inline double E_minus_lo_bound(double gamma_p, double e_e) {
        return (gamma_p * gamma_p + e_e * e_e) / (2.0 * gamma_p * e_e);
    }

    // omega_lo = E_minus_lo_bound + 1: smallest rest-frame photon energy for which
    // [E_minus_lo, omega-1] is non-empty. Always >= 2 (AM-GM), matching W's own
    // k>=2 precondition exactly -- no clamp needed.
    inline double omega_lo_bound(double gamma_p, double e_e) {
        const double s = gamma_p + e_e;
        return (s * s) / (2.0 * gamma_p * e_e);
    }

    // eps_lo in physical GeV: lab photon energy at which [omega_lo, 2*gamma_p*eps/m_e]
    // first becomes non-empty. Unlike the two bounds above, has NO floor as
    // gamma_p->infinity other than m_e^2/(4*E_e) -- the eps_lo>=epsilon_max guard
    // in dN_dEe is load-bearing, not just an optimization.
    inline double eps_lo_bound(double gamma_p, double E_e, double m_e) {
        const double num = gamma_p * m_e + E_e;
        return (num * num) / (4.0 * gamma_p * gamma_p * E_e);
    }

    namespace detail {

        // Innermost integral of eq.62/eq.67, shared by the general and the
        // Planckian-specialized paths:
        //     I(omega) = Int[E_minus_lo, omega-1] dE_- p_+ W(omega, E_-, xi)
        // (eq.62's dE_-/p_- against W_KA = alpha r0^2 (p_- p_+ / 2 omega^3) W
        // cancels the p_- COMPLETELY -- an earlier version kept a spurious
        // extra 1/p_-, which agreed accidentally at near-threshold windows
        // where the integral is dominated by E_- ~ 1.2-1.4 with 1/p_- ~ 1.1,
        // but suppressed wide-window/UHECR configurations by up to ~10x; see
        // the dN_dEe moment-consistency and Fig.10 tests that now pin this.)
        // E_plus = omega - E_minus -> 1 (p_plus -> 0) exactly at E_minus = omega-1.
        // W() has a removable 0/0 there (delta_plus_T/p_plus in its term6); the
        // quadrature below always samples its upper panel edge exactly, so
        // evaluating right at that edge returns NaN and poisons the whole integral
        // (silently, since NaN comparisons never trigger tolerance convergence --
        // the adaptive recursion then runs to max_depth on every affected panel,
        // at every level of nesting). The true limit there is finite (the p_plus
        // factor cancels W's 1/p_plus term analytically), so backing the bound
        // off by a relative epsilon drops only a measure-zero sliver, not real
        // integration accuracy.
        inline double E_minus_integral(double omega, double gamma_p, double e_e, double E_minus_lo,
                                        double abs_tol, double rel_tol, int max_depth, int panels) {
            // Mirror backoff at the LOWER edge: E_minus_lo == 1 exactly when
            // e_e == gamma_p (a query right at the spectral peak E_e =
            // gamma_p*m_e). p_minus = 0 there, so xi (division by p_minus) and
            // W's own 1/p_minus^2 terms are sampled as inf/NaN at the panel
            // edge and the adaptive recursion grinds to max_depth at every
            // nesting level -- measured as a multi-HOUR near-hang for a single
            // call (vs ~seconds one grid point away). The edge singularity is
            // integrable, so a relative backoff drops only a measure-zero
            // sliver, same as the upper edge's.
            const double E_minus_lo_eff = std::max(E_minus_lo, 1.0) * (1.0 + 1e-9);
            const double E_minus_hi = (omega - 1.0) * (1.0 - 1e-9);
            if (E_minus_lo_eff >= E_minus_hi) return 0.0;

            auto integrand = [&](double E_minus) {
                const double p_minus = std::sqrt(E_minus * E_minus - 1.0);
                const double E_plus   = omega - E_minus;
                const double p_plus   = std::sqrt(E_plus * E_plus - 1.0);
                const double xi = xi_cos_theta_minus(gamma_p, e_e, E_minus, p_minus);
                return p_plus * W(omega, E_minus, xi);
            };

            // Integrate in v = ln(E_- - 1) rather than ln(E_-): when
            // E_minus_lo == 1 (a query at exactly E_e = gamma_p*m_e) the
            // integrand carries a PHYSICAL integrable log singularity
            // p_+ W ~ 1/(E_- - 1) at the lower edge (W ~ 1/p_-^2 with
            // Delta_- pinned at e_e/gamma_p). In ln(E_-) coordinates the
            // adaptive recursion grinds against 1/x for minutes per call and
            // still truncates the edge at the backoff; in ln(E_- - 1)
            // coordinates that factor is FLAT (f * (E_- - 1) bounded), so the
            // quadrature is fast and captures the edge mass. Away from the
            // edge, v ~ ln E_- and the substitution behaves like the plain
            // log grid.
            auto g = [&](double v) {
                const double x = std::exp(v);
                return integrand(1.0 + x) * x;
            };
            return math::integrate(g, std::log(E_minus_lo_eff - 1.0), std::log(E_minus_hi - 1.0),
                                    abs_tol, rel_tol, max_depth, panels);
        }

    }   // namespace detail

    // Planckian-specialized single-proton spectrum: KA2008 eq.67, extended to a
    // finite epsilon_max. For a blackbody field the eps integral of eq.62 is
    // done analytically by parts (eq.66):
    //     deps f_ph(eps)/eps^2 = kT/(pi^2 hbar_c3) d ln(1 - e^{-eps/kT})
    // Swapping the (eps, omega) integration order (the region
    // {eps in [eps_lo, epsilon_max], omega in [omega_lo, 2 gamma_p eps/m_e]} is
    // exactly {omega in [omega_lo, 2 gamma_p epsilon_max/m_e],
    // eps in [omega m_e/(2 gamma_p), epsilon_max]} -- the inner eps bound is
    // >= eps_lo automatically since 2 gamma_p eps_lo/m_e == omega_lo) collapses
    // the triple integral to a double one:
    //
    //   dN/dE_e = alpha_r0sq_c m_e kT / (4 pi^2 hbar_c3 gamma_p^3)
    //             Int[omega_lo, omega_hi] domega/omega^2
    //               [ ln(1-e^{-epsilon_max/kT}) - ln(1-e^{-omega m_e/(2 gamma_p kT)}) ]
    //               I(omega)
    //
    // The bracket is the by-parts boundary term: it handles the finite
    // epsilon_max cutoff EXACTLY (no Wien-tail approximation), vanishes exactly
    // at omega_hi = 2 gamma_p epsilon_max/m_e, and reduces to eq.67's
    // -ln(1-e^{-omega/(2 gamma_p kT)}) weight as epsilon_max -> infinity. This
    // is an algebraic rearrangement of the same integrand the general path
    // evaluates, one adaptive-quadrature level cheaper -- results agree with
    // dN_dEe_general to within the quadrature tolerances (see
    // test/bh/spectrum_test.cpp's parity check).
    inline double dN_dEe_planck(double E_e, double gamma_p, double kT, double epsilon_max,
                                    double abs_tol = 1e-25, double rel_tol = 1e-4,
                                    int max_depth = 20, int panels = 24) {
        using namespace kaspectra::constants;

        const double e_e = E_e / m_e;
        const double eps_lo = eps_lo_bound(gamma_p, E_e, m_e);
        if (eps_lo >= epsilon_max) return 0.0;

        const double omega_lo = omega_lo_bound(gamma_p, e_e);
        const double omega_hi = 2.0 * gamma_p * epsilon_max / m_e;
        const double E_minus_lo = E_minus_lo_bound(gamma_p, e_e);

        // ln(1 - e^{-x}) via log1p for accuracy at both ends; exp underflow to
        // 0.0 (x >~ 745) gives log1p(-0.0) == 0.0, the correct limit.
        const double L_max = std::log1p(-std::exp(-epsilon_max / kT));

        auto omega_integrand = [&](double omega) {
            const double L_omega = std::log1p(-std::exp(-omega * m_e / (2.0 * gamma_p * kT)));
            const double weight = L_max - L_omega;
            if (weight <= 0.0) return 0.0;
            const double inner = detail::E_minus_integral(omega, gamma_p, e_e, E_minus_lo,
                                                            abs_tol, rel_tol, max_depth, panels);
            return weight * inner / (omega * omega);
        };

        return (alpha_r0sq_c * m_e * kT / (4.0 * pi * pi * hbar_c3 * gamma_p * gamma_p * gamma_p)) *
                math::integrate_log(omega_integrand, omega_lo, omega_hi, abs_tol, rel_tol, max_depth, panels);
    }

    // Single-proton lab-frame e+ (or e-, identical by symmetry -- KA2008 states this
    // explicitly after eq.67) differential spectrum, KA2008 eq.62.
    //
    //   dN/dE_e = 1/(2 gamma_p^3) Int[eps_lo,eps_max] deps f_ph(eps)/eps^2
    //             Int[omega_lo, 2 gamma_p eps] domega omega
    //             Int[E_lo, omega-1] dE_-/p_- W_KA(omega, E_-, xi)
    //
    // W_KA == d^2(sigma)/dE_- dcos(theta_-) is Blumenthal's FULL cross section
    // (eq.58) -- NOT this file's W(), which returns only the dimensionless
    // bracket (see its own doc comment). Relation:
    //
    //   W_KA(omega,E_-,cos) = alpha*r0^2 * (p_-*p_+ / (2*omega^3)) * W(omega,E_-,cos)
    //
    // Substituting: omega*domega combines with 1/omega^3 into domega/omega^2, and
    // the p_- in W_KA's prefactor cancels eq.62's dE_-/p_- COMPLETELY, leaving p_+:
    //
    //   dN/dE_e = (alpha*r0^2)/(4 gamma_p^3) Int deps f_ph(eps)/eps^2
    //             Int domega/omega^2 Int dE_- p_+ W(omega,E_-,xi)
    //
    // eq.62 uses Blumenthal's c=hbar=m_e=1 convention; restoring physical units
    // (E_e, eps in GeV, gamma_p dimensionless) needs e_e=E_e/m_e substituted
    // throughout, plus one m_e from the eps-Jacobian, plus one c_light (folded
    // into alpha_r0sq_c). Net prefactor: alpha_r0sq_c * m_e / (4 gamma_p^3).
    // Verified by the moment-consistency check in test/bh/spectrum_test.cpp
    // (integral of dN_dEe over E_e reproduces interaction_rate; energy-weighted
    // integral reproduces energy_loss_rate).
    //
    // Default tolerances are deliberately looser than bh/source.hpp's single-
    // integral functions (1e-10/1e-8/50/32 there): this is a 3-level nested
    // integral, so every additional bit of relative precision demanded at the
    // outer level multiplies the cost of the (already expensive) inner two
    // integrals at every recursion step. Measured directly: at UHECR-scale
    // gamma_p (~1e11) with a CMB target field, rel_tol=1e-6 does not finish in
    // any practical time (adaptive recursion compounds across all 3 levels),
    // while rel_tol=1e-4 with max_depth=20 completes in ~6s and already agrees
    // with the moment-consistency check (see test file) to within a few
    // percent -- ample for the ~0.1-1% intrinsic accuracy of the underlying
    // W()/eq.62 physics itself. abs_tol is set far below any physically
    // meaningful scale so it never dominates over rel_tol.
    inline double dN_dEe_general(double E_e, double gamma_p, const io::PhotonField& f_ph, double epsilon_max,
                                    double abs_tol = 1e-25, double rel_tol = 1e-4,
                                    int max_depth = 20, int panels = 24) {
        using namespace kaspectra::constants;

        const double e_e = E_e / m_e;
        const double eps_lo = eps_lo_bound(gamma_p, E_e, m_e);
        if (eps_lo >= epsilon_max) return 0.0;

        const double omega_lo = omega_lo_bound(gamma_p, e_e);
        const double E_minus_lo = E_minus_lo_bound(gamma_p, e_e);

        auto eps_integrand = [&](double eps) {
            // Cheap short-circuit: at UHECR-scale gamma_p, most of [eps_lo, epsilon_max]
            // sits far out on the photon field's own suppressed tail (e.g. a blackbody's
            // Wien tail, where f_ph underflows to exact 0.0 well before epsilon_max is
            // reached), while the omega/E_minus domain size below GROWS with eps -- so
            // without this check, the most expensive nested-integral evaluations land
            // exactly where the physical contribution is already zero. Always safe:
            // f_ph(eps)==0.0 makes the whole integrand identically zero regardless of
            // the (expensive) omega integral's value.
            const double f_ph_eps = f_ph(eps);
            if (f_ph_eps == 0.0) return 0.0;

            const double omega_hi = 2.0 * gamma_p * eps / m_e;
            if (omega_lo >= omega_hi) return 0.0;

            auto omega_integrand = [&](double omega) {
                return detail::E_minus_integral(omega, gamma_p, e_e, E_minus_lo,
                                                abs_tol, rel_tol, max_depth, panels) / (omega * omega);
            };

            return (f_ph_eps / (eps * eps)) *
                    math::integrate_log(omega_integrand, omega_lo, omega_hi, abs_tol, rel_tol, max_depth, panels);
        };

        return (alpha_r0sq_c * m_e / (4.0 * gamma_p * gamma_p * gamma_p)) *
                math::integrate_log(eps_integrand, eps_lo, epsilon_max, abs_tol, rel_tol, max_depth, panels);
    }

    // Field-agnostic fast path: the SAME integration-order swap that powers
    // dN_dEe_planck, but with the eps cumulative
    //     G(x) = Int[x, epsilon_max] f_ph(eps)/eps^2 deps
    // evaluated numerically instead of analytically (for a Planckian field G
    // has the closed form eq.66 gives; for any other field it doesn't, but
    // nothing about the swap itself required the closed form):
    //
    //   dN/dE_e = alpha_r0sq_c m_e / (4 gamma_p^3)
    //             Int[omega_lo, omega_hi] domega/omega^2
    //                 G(omega m_e/(2 gamma_p)) I(omega)
    //
    // G is set up once per call as suffix sums of per-segment adaptive
    // integrals on a log grid over [eps_lo, epsilon_max] (f_ph evaluations
    // only -- cheap), and each query point is then completed EXACTLY by one
    // small adaptive integral over the remaining sub-segment (never wider
    // than one grid step), so no interpolation error enters anywhere: the
    // result differs from dN_dEe_general only by quadrature tolerances, at
    // one adaptive-nesting level less cost. See the parity tests in
    // test/bh/spectrum_test.cpp (power-law and tabulated fields).
    inline double dN_dEe_fast(double E_e, double gamma_p, const io::PhotonField& f_ph, double epsilon_max,
                                double abs_tol = 1e-25, double rel_tol = 1e-4,
                                int max_depth = 20, int panels = 24) {
        using namespace kaspectra::constants;

        const double e_e = E_e / m_e;
        const double eps_lo = eps_lo_bound(gamma_p, E_e, m_e);
        if (eps_lo >= epsilon_max) return 0.0;

        const double omega_lo = omega_lo_bound(gamma_p, e_e);
        const double omega_hi = 2.0 * gamma_p * epsilon_max / m_e;
        const double E_minus_lo = E_minus_lo_bound(gamma_p, e_e);

        auto f_over_eps2 = [&](double eps) { return f_ph(eps) / (eps * eps); };

        // Log grid over [eps_lo, epsilon_max], ~24 points/decade (min 4
        // segments); suffix[i] = Int[grid[i], epsilon_max].
        const double decades = std::log10(epsilon_max / eps_lo);
        const int n_seg = std::max(4, static_cast<int>(std::ceil(decades * 24.0)));
        std::vector<double> grid(static_cast<std::size_t>(n_seg) + 1);
        const double log_lo = std::log(eps_lo), log_hi = std::log(epsilon_max);
        for (int i = 0; i <= n_seg; ++i) {
            grid[static_cast<std::size_t>(i)] =
                std::exp(log_lo + (log_hi - log_lo) * static_cast<double>(i) / n_seg);
        }
        grid.front() = eps_lo;
        grid.back() = epsilon_max;

        std::vector<double> suffix(grid.size(), 0.0);
        for (int i = n_seg - 1; i >= 0; --i) {
            const auto ui = static_cast<std::size_t>(i);
            suffix[ui] = suffix[ui + 1] +
                math::integrate_log(f_over_eps2, grid[ui], grid[ui + 1], abs_tol, rel_tol, max_depth, 4);
        }

        auto G = [&](double x) {
            if (x >= epsilon_max) return 0.0;
            if (x <= eps_lo) return suffix.front();
            const auto it = std::upper_bound(grid.begin(), grid.end(), x);
            const auto hi = static_cast<std::size_t>(it - grid.begin());
            return suffix[hi] + math::integrate_log(f_over_eps2, x, grid[hi], abs_tol, rel_tol, max_depth, 4);
        };

        auto omega_integrand = [&](double omega) {
            const double g = G(omega * m_e / (2.0 * gamma_p));
            if (g <= 0.0) return 0.0;
            const double inner = detail::E_minus_integral(omega, gamma_p, e_e, E_minus_lo,
                                                            abs_tol, rel_tol, max_depth, panels);
            return g * inner / (omega * omega);
        };

        return (alpha_r0sq_c * m_e / (4.0 * gamma_p * gamma_p * gamma_p)) *
                math::integrate_log(omega_integrand, omega_lo, omega_hi, abs_tol, rel_tol, max_depth, panels);
    }

    // Public entry point: dispatches to the Planckian-specialized eq.67 path
    // when f_ph is a BlackbodyPhotonField, otherwise to the field-agnostic
    // dN_dEe_fast -- both are exact algebraic rearrangements of the eq.62
    // triple integral, one adaptive-quadrature level cheaper, not
    // approximations. dN_dEe_general stays public as the straightforward
    // eq.62 reference path.
    inline double dN_dEe(double E_e, double gamma_p, const io::PhotonField& f_ph, double epsilon_max,
                            double abs_tol = 1e-25, double rel_tol = 1e-4,
                            int max_depth = 20, int panels = 24) {
        if (const auto* bb = dynamic_cast<const io::BlackbodyPhotonField*>(&f_ph)) {
            return dN_dEe_planck(E_e, gamma_p, bb->kT(), epsilon_max, abs_tol, rel_tol, max_depth, panels);
        }
        return dN_dEe_fast(E_e, gamma_p, f_ph, epsilon_max, abs_tol, rel_tol, max_depth, panels);
    }

    // Minimum proton energy for which lab electron energy E_e is kinematically
    // reachable at all, given photon field cutoff epsilon_max. Solves
    // eps_lo_bound(gamma_p,E_e,m_e) == epsilon_max for gamma_p: quadratic
    // A*g^2+B*g+C=0 with A=m_e^2-4*epsilon_max*E_e, B=2*m_e*E_e, C=E_e^2.
    // Returns +infinity if no gamma_p reaches eps_lo<=epsilon_max for this E_e
    // (i.e. epsilon_max*E_e <= m_e^2/4).
    inline double E_p_min(double E_e, double epsilon_max) {
        using namespace kaspectra::constants;

        const double A = m_e * m_e - 4.0 * epsilon_max * E_e;
        if (A >= 0.0) return std::numeric_limits<double>::infinity();

        const double B = 2.0 * m_e * E_e;
        const double C = E_e * E_e;
        const double disc = B * B - 4.0 * A * C;
        const double q = -(B + std::sqrt(disc)) / 2.0;   // stable form, B>0 always here

        return (q / A) * m_p;
    }

    // Population-level e+ (or e-) lab-frame spectrum: outer integral of dN_dEe
    // over the proton spectrum J_p(E_p), mirroring q_pair_rate/q_pair_energy_loss.
    //
    // abs_tol/rel_tol/max_depth/panels here govern ONLY this outer E_p integral,
    // not dN_dEe's own internal (eps,omega,E_minus) integrals -- dN_dEe is called
    // with its own defaults instead of threading these through, because a single
    // dN_dEe evaluation already costs seconds at UHECR scale (measured directly:
    // ~6s at gamma_p~1e11 against a CMB field); reusing dN_dEe's own tuned inner
    // defaults here and keeping the OUTER integral's own panel count small keeps
    // the total evaluation count (panels * dN_dEe calls) tractable. J_p(E_p) is
    // expected to be smooth (a power law or similar), so a loose outer rel_tol
    // and few panels are sufficient -- no sharp feature is expected in E_p beyond
    // the domain edges already guarded by E_p_lo/E_p_max.
    inline double q_pair_spectrum(double E_e, const io::ProtonSpectrum& J_p, const io::PhotonField& f_ph,
                                    double E_p_max, double epsilon_max,
                                    double abs_tol = 1e-25, double rel_tol = 1e-3,
                                    int max_depth = 10, int panels = 6) {
        using namespace kaspectra::constants;

        const double E_p_lo = E_p_min(E_e, epsilon_max);
        if (E_p_lo >= E_p_max) return 0.0;

        auto integrand = [&](double E_p) {
            return J_p(E_p) * dN_dEe(E_e, E_p / m_p, f_ph, epsilon_max);
        };

        return math::integrate_log(integrand, E_p_lo, E_p_max, abs_tol, rel_tol, max_depth, panels);
    }

}   // namespace kaspectra::bh