#pragma once

#include "kaspectra/bh/spectrum.hpp"
#include "kaspectra/io/photon_field.hpp"
#include "kaspectra/io/proton_spectrum.hpp"
#include "kaspectra/pgamma/source.hpp"
#include "kaspectra/pp/source.hpp"
#include "kaspectra/species.hpp"

namespace kaspectra {

    // Combined secondary spectrum from every implemented channel, per shared
    // Species [particles GeV^-1 cm^-3 s^-1]:
    //
    //   pp      (KAB2006)          -- all 7 species, needs the target gas
    //                                 density n_H [cm^-3]
    //   pgamma  (KA2008 photomeson)-- all 7 species, needs the photon field
    //   bh      (Bethe-Heitler)    -- e+/e- only (pair production makes no
    //                                 photons or neutrinos)
    //
    // Each channel is evaluated with its OWN tuned default tolerances (they
    // differ by orders of magnitude for good reasons documented at each
    // definition -- bh's triple integral cannot afford pp's 1e-8 rel_tol);
    // callers needing custom tolerances should sum the channels themselves.
    // Zero-cost opt-outs: n_H == 0 skips pp, and an epsilon_max below every
    // channel's threshold effectively skips the photon-target channels.
    inline double q_total_species(Species s, double E,
                                    const io::ProtonSpectrum& J_p, double n_H,
                                    const io::PhotonField& f_ph,
                                    double E_p_max, double E_pi_max, double epsilon_max) {
        double q = 0.0;
        if (n_H > 0.0) {
            q += pp::q_species(s, E, J_p, n_H, E_p_max, E_pi_max);
        }
        q += pgamma::q_species(s, E, J_p, f_ph, E_p_max, epsilon_max);
        if (s == Species::Positron || s == Species::Electron) {
            q += bh::q_pair_spectrum(E, J_p, f_ph, E_p_max, epsilon_max);
        }
        return q;
    }

}   // namespace kaspectra
