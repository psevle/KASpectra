#pragma once

#include <cmath>
#include <memory>
#include <vector>

#include "kaspectra/math/interp1d.hpp"

namespace kaspectra::io {

    // J_p(E_p): differentail proton flux/number density vs energy [GeV]
    // KASpectra does not enforce or convert flux units (e.g. GeV^-1 cm^-2 s^-1 sr^-1
    // vs cm^-3) - normalisation and units are entirely the users responsibility
    class ProtonSpectrum {
        public:
            virtual ~ProtonSpectrum() = default;
            virtual double operator()(double E_p) const = 0;
    };

    class PowerLawSpectrum final : public ProtonSpectrum {
        public:
            PowerLawSpectrum(double norm, double index, double E_ref = 1.0)
                : norm_(norm), index_(index), E_ref_(E_ref) {}
            
            double operator()(double E_p) const override {
                return norm_ * std::pow(E_p / E_ref_, -index_);
            }
        
        private:
            double norm_, index_, E_ref_;
    };

    class ExpCutoffPowerLawSpectrum final : public ProtonSpectrum {
        public:
            ExpCutoffPowerLawSpectrum(double norm, double index, double E_cutoff, double E_ref = 1.0)
                : norm_(norm), index_(index), E_cutoff_(E_cutoff), E_ref_(E_ref) {}
            
            double operator()(double E_p) const override {
                return norm_ * std::pow(E_p / E_ref_, -index_) * std::exp(-E_p / E_cutoff_);
            }
        
        private:
            double norm_, index_, E_cutoff_, E_ref_;
    };

    // norm is defined at E_ref (must lie below E_break); the high-energy segmet's
    // norm is derived so the two power laws agree at E_break - continuous by
    // construction rather than by a fitted second parameter
    class BrokenPowerLawSpectrum final : public ProtonSpectrum {
        public:
            BrokenPowerLawSpectrum(double norm, double index_lo, double index_hi,
                                    double E_break, double E_ref = 1.0)
                : norm_(norm), index_lo_(index_lo), index_hi_(index_hi), E_break_(E_break), E_ref_(E_ref),
                    norm_hi_(norm * std::pow(E_break / E_ref, -index_lo)) {}
            
            double operator()(double E_p) const override {
                if (E_p < E_break_) {
                    return norm_ * std::pow(E_p / E_ref_, -index_lo_);
                }
                return norm_hi_ * std::pow(E_p / E_break_, -index_hi_);
            }
        
            private:
                double norm_, index_lo_, index_hi_, E_break_, E_ref_, norm_hi_;
    };

    class TabulatedSpectrum final : public ProtonSpectrum {
        public:
            TabulatedSpectrum(const std::vector<double>& E, const std::vector<double>& J)
                : interp_(E, J, math::InterpMode::LogLog) {}
            
            double operator()(double E_p) const override { return interp_(E_p); }
        
        private:
            math::Interpolator1D interp_;
    };

    // Usage:
    //   PowerLawSpectrum J_p(1.0, 2.0);
    //   double flux = integrate_log([&](double E){ return J_p(E) / E; }, 1.0, 1e6);

} // namespace kaspectra::io