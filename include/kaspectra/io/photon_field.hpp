#pragma once

#include <cmath>
#include <vector>

#include "kaspectra/constants.hpp"
#include "kaspectra/math/interp1d.hpp"

namespace kaspectra::io {

    // f_ph(epsilon): differential photon number density vs photon energy [GeV].
    // f_ph(epsilon) dEpsilon is the number of photons per 1 cm^3 in the energy
    // interval [epsilon, epsilon + dEpsilon]
    class PhotonField {
        public:
            virtual ~PhotonField() = default;
            virtual double operator()(double epsilon) const = 0;
    };

    // Isotropic blackbody (Planckian) photon field.
    // Eq.(65) (natural units, hbar=c=1):
    //     f_ph(epsilon) = (1/pi^2) * epsilon^2 / (exp(epsilon/kT) - 1)
    // Restored to GeV/cm units via constants::hbar_c3 (see constants.hpp).
    // T is supplied in Kelvin; kT is converted to GeV via constants::k_boltzmann.
    class BlackbodyPhotonField final :  public PhotonField {
        public:
            explicit BlackbodyPhotonField(double T_kelvin)
                : kT_(T_kelvin * constants::k_boltzmann) {}
            
            double operator()(double epsilon) const override {
                return (epsilon * epsilon) / (constants::pi * constants::pi * constants::hbar_c3 * std::expm1(epsilon / kT_));
            }

            // kT in GeV. Exposed so Planckian-specialized fast paths (bh::dN_dEe's
            // eq.67 branch) can recover the temperature they need to collapse the
            // epsilon integral analytically.
            double kT() const { return kT_; }

            private:
                double kT_;
    };

    // Power-law photon field with exponential cutoff
    // Same shape as ExpCutoffPowerLawSpectrum for protons
    class PowerLawPhotonField final : public PhotonField {
        public:
            PowerLawPhotonField(double norm, double index, double epsilon_cutoff, double epsilon_ref = 1.0)
                : norm_(norm), index_(index), epsilon_cutoff_(epsilon_cutoff), epsilon_ref_(epsilon_ref) {}
            
            double operator()(double epsilon) const override {
                return norm_ * std::pow(epsilon / epsilon_ref_, -index_) * std::exp(-epsilon / epsilon_cutoff_);
            }
        
        private:
            double norm_, index_, epsilon_cutoff_, epsilon_ref_;
    };

    // Tabulated (interpolated) photon field, mirroring TabulatedSpectrum
    class TabulatedPhotonField final : public PhotonField {
        public:
            TabulatedPhotonField(const std::vector<double>& epsilon, const std::vector<double>& f_ph)
                : interp_(epsilon, f_ph, math::InterpMode::LogLog) {}
            
            double operator()(double epsilon) const override {
                return interp_(epsilon);
            }
        
        private:
            math::Interpolator1D interp_;
    };

    // Sum of photon fields -- the standard astrophysical target is a
    // superposition (CMB + infrared + starlight), and every routine in this
    // library takes a single PhotonField&. NON-OWNING: stores pointers to the
    // added fields, which the caller must keep alive for the composite's
    // lifetime (the Python binding enforces this with keep_alive). Note that
    // a composite is not a BlackbodyPhotonField even if every component is,
    // so bh::dN_dEe routes it through the field-agnostic fast path rather
    // than the Planckian one.
    class CompositePhotonField final : public PhotonField {
        public:
            CompositePhotonField() = default;

            void add(const PhotonField& field) { fields_.push_back(&field); }

            double operator()(double epsilon) const override {
                double sum = 0.0;
                for (const PhotonField* f : fields_) sum += (*f)(epsilon);
                return sum;
            }

            std::size_t size() const { return fields_.size(); }

        private:
            std::vector<const PhotonField*> fields_;
    };

    // Usage:
    //   kaspectra::io::BlackbodyPhotonField cmb(2.725); // 2.725 K CMB
    //   double n = cmb(1e-13); // dn/depsilon near CMB peak [cm^-3 GeV^-1]

} // namespace kaspectra::io