#pragma once

#include <algorithm>
#include <vector>
#include "kaspectra/math/interp1d.hpp"

namespace kaspectra::pgamma{

    // (B, s, delta) fit paramters at a given eta/eta_0
    // B is in cm^3/s; s and delta are dimless
    struct PhiTableParams{
        double B;
        double s;
        double delta;
    };

    namespace detail {

        //Table I (gamma) eta/eta_0 grid
        inline const std::vector<double>& table_I_eta() {
            static const std::vector<double> eta = {
                1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                20.0, 30.0, 40.0, 100.0
            };
            return eta;
        }

        // Table II (e+, numubar, numu, nue) eta/eta_0 grid
        inline const std::vector<double>& table_II_eta() {
            static const std::vector<double> eta = {
                1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                30.0, 100.0
            };
            return eta;
        }

        // Table III (e-, nuebar) eta/eta_0 grid
        inline const std::vector<double>& table_III_eta() {
            static const std::vector<double> eta = {
                3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 30.0, 100.0
            };
            return eta;
        }

    } // namespace detail

    // -----------------------------------------------------------------------
    // Table I - gamma-ray (KA2008 eq. 27, Table I, p.4)
    //
    // Edge behavior: the paper states B_gamma -> 0 exactly at eta/eta_0 = 1
    // (text below eq. 29) but tabulates nothing between 1.0 and the first row
    // at 1.1. B is therefore blended linearly from 0 at rho=1 to the row-1.1
    // value (continuous at 1.1 by construction) instead of clamping high --
    // the previous clamp overestimated the near-threshold amplitude across
    // ~10% of the fit's domain. s/delta keep the clamp: they only shape a
    // spectrum whose amplitude B already vanishes at threshold, and the paper
    // gives no threshold statement for them.
    // Above eta/eta_0 = 100 (the last row), Interpolator1D clamps as before.
    // -----------------------------------------------------------------------
    inline PhiTableParams lookup_table_I(double eta_over_eta0) {
        static const math::Interpolator1D s_interp(
            detail::table_I_eta(),
            { 0.0768, 0.106, 0.182, 0.201, 0.219, 0.216, 0.233, 0.233, 0.248, 0.244,
            0.188, 0.131, 0.120, 0.107, 0.102, 0.0932, 0.0838, 0.0761, 0.107, 0.0928, 0.0772, 0.0479 },
            math::InterpMode::Linear);
        static const math::Interpolator1D delta_interp(
            detail::table_I_eta(),
            { 0.544, 0.540, 0.750, 0.791, 0.788, 0.831, 0.839, 0.825, 0.805, 0.779,
            1.23, 1.82, 2.05, 2.19, 2.23, 2.29, 2.37, 2.43,
            2.27, 2.33, 2.42, 2.59 },
            math::InterpMode::Linear);
        static const math::Interpolator1D B_interp(
            detail::table_I_eta(),
            { 2.86e-19, 2.24e-18, 5.61e-18, 1.02e-17, 1.60e-17, 2.23e-17, 3.10e-17,
            4.07e-17, 5.30e-17, 6.74e-17, 1.51e-16, 1.24e-16, 1.37e-16, 1.62e-16,
            1.71e-16, 1.78e-16, 1.84e-16, 1.93e-16, 4.74e-16, 7.70e-16, 1.06e-15, 2.73e-15 },
            math::InterpMode::Linear);

        double B = B_interp(eta_over_eta0);
        if (eta_over_eta0 < 1.1) {
            B *= std::max(0.0, (eta_over_eta0 - 1.0) / 0.1);
        }
        return { B, s_interp(eta_over_eta0), delta_interp(eta_over_eta0) };
    }

    // -----------------------------------------------------------------------
    // Table II - e+, numubar, numu, nue (KA2008 eq. 31/34/35/36/37, Table II, p.6)
    //
    // One lookup function per species: each maps 1:1 to the eventual per-species Phi_*
    // function in secondary_spectra.hpp. Edge behavior is the same
    // clamp-at-row-1.1/row-100 policy as Table I
    // -----------------------------------------------------------------------
    inline PhiTableParams lookup_table_II_e_plus(double eta_over_eta0) {
        static const math::Interpolator1D s_interp(
            detail::table_II_eta(),
            { 0.367, 0.282, 0.260, 0.239, 0.224, 0.207, 0.198, 0.193, 0.187, 0.181,
            0.122, 0.106, 0.0983, 0.0875, 0.0830, 0.0783, 0.0735, 0.0644,
            0.0333, 0.0224 },
            math::InterpMode::Linear);
        static const math::Interpolator1D delta_interp(
            detail::table_II_eta(),
            { 3.12, 2.96, 2.83, 2.76, 2.69, 2.66, 2.62, 2.56, 2.52, 2.49,
            2.48, 2.50, 2.46, 2.46, 2.44, 2.44, 2.45, 2.50, 2.77, 2.86 },
            math::InterpMode::Linear);
        static const math::Interpolator1D B_interp(
            detail::table_II_eta(),
            { 8.09e-19, 7.70e-18, 2.05e-17, 3.66e-17, 5.48e-17, 7.39e-17, 9.52e-17,
            1.20e-16, 1.47e-16, 1.75e-16, 3.31e-16, 4.16e-16, 5.57e-16, 6.78e-16,
            7.65e-16, 8.52e-16, 9.17e-16, 9.57e-16, 3.07e-15, 1.58e-14 },
            math::InterpMode::Linear);

        return { B_interp(eta_over_eta0), s_interp(eta_over_eta0), delta_interp(eta_over_eta0) };
    }

    inline PhiTableParams lookup_table_II_numubar(double eta_over_eta0) {
        static const math::Interpolator1D s_interp(
            detail::table_II_eta(),
            { 0.365, 0.287, 0.250, 0.238, 0.220, 0.206, 0.197, 0.193, 0.187, 0.178,
            0.123, 0.106, 0.0944, 0.0829, 0.0801, 0.0752, 0.0680, 0.0615, 0.0361, 0.0228 },
            math::InterpMode::Linear);
        static const math::Interpolator1D delta_interp(
            detail::table_II_eta(),
            { 3.09, 2.96, 2.89, 2.76, 2.71, 2.67, 2.62, 2.56, 2.52, 2.51,
            2.48, 2.50, 2.57, 2.58, 2.54, 2.53, 2.56, 2.60, 2.78, 2.88 },
            math::InterpMode::Linear);
        static const math::Interpolator1D B_interp(
            detail::table_II_eta(),
            { 8.09e-19, 7.70e-18, 1.99e-17, 3.62e-17, 5.39e-17, 7.39e-17, 9.48e-17,
            1.20e-16, 1.47e-16, 1.74e-16, 3.38e-16, 5.17e-16, 7.61e-16, 9.57e-16,
            1.11e-15, 1.25e-15, 1.36e-15, 1.46e-15, 5.87e-15, 3.10e-14 },
            math::InterpMode::Linear);

        return { B_interp(eta_over_eta0), s_interp(eta_over_eta0), delta_interp(eta_over_eta0) };
    }

    // Note (paper text, p.5): for eta < 2.14*eta_0 only pi+ is produced, so
    // e+ and numubar come from the same mu+ decay and their table entries
    // should be identical up to simulation noise - visible above in the
    // near-equal e+/numubar columns for eta/eta_0 <= 2.0. Above 2.14*eta_0,
    // pi- production opens a new numubar channel and the columns diverge
    // significantly.
    inline PhiTableParams lookup_table_II_numu(double eta_over_eta0) {
        static const math::Interpolator1D s_interp(
            detail::table_II_eta(),
            { 0.0, 0.0778, 0.242, 0.377, 0.440, 0.450, 0.461, 0.451, 0.464, 0.446,
            0.366, 0.249, 0.204, 0.174, 0.156, 0.140, 0.121, 0.107, 0.0705, 0.0463 },
            math::InterpMode::Linear);
        static const math::Interpolator1D delta_interp(
            detail::table_II_eta(),
            { 0.0, 0.306, 0.792, 1.09, 1.06, 0.953, 0.956, 0.922, 0.912, 0.940,
            1.49, 2.03, 2.18, 2.24, 2.28, 2.32, 2.39, 2.46, 2.53, 2.62 },
            math::InterpMode::Linear);
        static const math::Interpolator1D B_interp(
            detail::table_II_eta(),
            { 1.08e-18, 9.91e-18, 2.47e-17, 4.43e-17, 6.70e-17, 9.04e-17, 1.18e-16,
            1.32e-16, 1.77e-16, 2.11e-16, 3.83e-16, 5.09e-16, 7.26e-16, 9.26e-16,
            1.07e-15, 1.19e-15, 1.29e-15, 1.40e-15, 5.65e-15, 3.01e-14 },
            math::InterpMode::Linear);

        return { B_interp(eta_over_eta0), s_interp(eta_over_eta0), delta_interp(eta_over_eta0) };
    }

    inline PhiTableParams lookup_table_II_nue(double eta_over_eta0) {
        static const math::Interpolator1D s_interp(
            detail::table_II_eta(),
            { 0.768, 0.569, 0.491, 0.395, 0.31, 0.323, 0.305, 0.285, 0.270, 0.259,
            0.158, 0.129, 0.113, 0.0996, 0.0921, 0.0861, 0.0800, 0.0723, 0.0411, 0.0283 },
            math::InterpMode::Linear);
        static const math::Interpolator1D delta_interp(
            detail::table_II_eta(),
            { 2.49, 2.35, 2.41, 2.45, 2.45, 2.43, 2.40, 2.39, 2.37, 2.35,
            2.42, 2.46, 2.45, 2.46, 2.46, 2.45, 2.47, 2.51, 2.70, 2.77 },
            math::InterpMode::Linear);
        static const math::Interpolator1D B_interp(
            detail::table_II_eta(),
            { 9.43e-19, 9.22e-18, 2.35e-17, 4.20e-17, 6.26e-17, 8.57e-17, 1.13e-16,
            1.39e-16, 1.70e-16, 2.05e-16, 3.81e-16, 4.74e-16, 6.30e-16, 7.65e-16,
            8.61e-16, 9.61e-16, 1.03e-15, 1.10e-15, 3.55e-15, 1.86e-14 },
            math::InterpMode::Linear);

        return { B_interp(eta_over_eta0), s_interp(eta_over_eta0), delta_interp(eta_over_eta0) };
    }

    // -----------------------------------------------------------------------
    // Table III - e-, nuebar (KA2008 eq. 31/40/41, Table III, p.8)
    //
    // IMPORTANT domain note: this table's grid starts at eta/eta_0 = 3.0.
    // The two-pion production threshold (eq. 38) is eta/eta_0 = 2.14, and
    // eq. 41's psi Heaviside step Theta(rho-4) turns on at rho = 4 - neither
    // of those thresholds coincides with the table's own data floor of 3.0.
    // Querying below 3.0 with Interpolator1D's clamp policy returns the
    // row-3.0 values rather than zero or an extrapolated fit
    // This table intentionally does NOT special-case eta/eta_0 < 3.0
    // that decision belongs one layer up.
    // -----------------------------------------------------------------------
    inline PhiTableParams lookup_table_III_e_minus(double eta_over_eta0) {
        static const math::Interpolator1D s_interp(
            detail::table_III_eta(),
            { 0.658, 0.348, 0.256, 0.220, 0.217, 0.220, 0.217, 0.192, 0.125, 0.0507 },
            math::InterpMode::Linear);
        static const math::Interpolator1D delta_interp(
            detail::table_III_eta(),
            { 3.09, 2.81, 2.39, 2.27, 2.13, 2.20, 2.13, 2.19, 2.27, 2.63 },
            math::InterpMode::Linear);
        static const math::Interpolator1D B_interp(
            detail::table_III_eta(),
            { 6.43e-19, 9.91e-18, 1.24e-16, 2.67e-16, 3.50e-16,
            4.03e-16, 4.48e-16, 4.78e-16, 1.64e-15, 4.52e-15 },
            math::InterpMode::Linear);

        return { B_interp(eta_over_eta0), s_interp(eta_over_eta0), delta_interp(eta_over_eta0) };
    }

    inline PhiTableParams lookup_table_III_nuebar(double eta_over_eta0) {
        static const math::Interpolator1D s_interp(
            detail::table_III_eta(),
            { 0.985, 0.378, 0.31, 0.327, 0.308, 0.292, 0.260, 0.233, 0.135, 0.0770 },
            math::InterpMode::Linear);
        static const math::Interpolator1D delta_interp(
            detail::table_III_eta(),
            { 2.63, 2.98, 2.31, 2.11, 2.03, 1.98, 2.02, 2.07, 2.24, 2.40 },
            math::InterpMode::Linear);
        static const math::Interpolator1D B_interp(
            detail::table_III_eta(),
            { 6.61e-19, 9.74e-18, 1.34e-16, 2.91e-16, 3.81e-16,
            4.48e-16, 4.83e-16, 5.13e-16, 1.75e-15, 5.48e-15 },
            math::InterpMode::Linear);

        return { B_interp(eta_over_eta0), s_interp(eta_over_eta0), delta_interp(eta_over_eta0) };
    }

} // namespace kaspectra::pgamma