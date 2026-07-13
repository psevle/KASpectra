#include <pybind11/pybind11.h>

#include "kaspectra/species.hpp"
#include "kaspectra/constants.hpp"
#include "kaspectra/math/integrate.hpp"
#include "kaspectra/total.hpp"

namespace py = pybind11;

void bind_io(py::module_&);
void bind_pp(py::module_&);
void bind_pgamma(py::module_&);
void bind_bh(py::module_&);

PYBIND11_MODULE(_kaspectra, m) {
    m.doc() = "KASpectra: cosmic-ray secondary spectra (pp, pgamma, Bethe-Heitler). "
              "Low-level typed extension -- import kaspectra instead of this module directly.";

    py::enum_<kaspectra::Species>(m, "Species",
        "Secondary particle species, shared by the pp and pgamma channels.")
        .value("Gamma", kaspectra::Species::Gamma)
        .value("Positron", kaspectra::Species::Positron)
        .value("Electron", kaspectra::Species::Electron)
        .value("NuMu", kaspectra::Species::NuMu)
        .value("NuMuBar", kaspectra::Species::NuMuBar)
        .value("NuE", kaspectra::Species::NuE)
        .value("NuEBar", kaspectra::Species::NuEBar);

    // GeV / natural units throughout, matching constants.hpp.
    auto constants = m.def_submodule("constants", "Physical constants (see kaspectra/constants.hpp)");
    constants.attr("m_p")                = kaspectra::constants::m_p;
    constants.attr("m_pi_charged")       = kaspectra::constants::m_pi_charged;
    constants.attr("m_mu")               = kaspectra::constants::m_mu;
    constants.attr("m_e")                = kaspectra::constants::m_e;
    constants.attr("E_threshold_pp")     = kaspectra::constants::E_threshold_pp;
    constants.attr("c_light")            = kaspectra::constants::c_light;
    constants.attr("E_delta_approx_max") = kaspectra::constants::E_delta_approx_max;
    constants.attr("pi")                 = kaspectra::constants::pi;
    constants.attr("k_boltzmann")        = kaspectra::constants::k_boltzmann;
    constants.attr("hbar_c3")            = kaspectra::constants::hbar_c3;
    constants.attr("alpha_r0sq_c")       = kaspectra::constants::alpha_r0sq_c;
    constants.attr("r")                  = kaspectra::constants::r;
    constants.attr("eta_0")              = kaspectra::constants::eta_0;
    constants.attr("cm_per_Mpc")         = kaspectra::constants::cm_per_Mpc;
    constants.attr("seconds_per_year")   = kaspectra::constants::seconds_per_year;

    // Quadrature diagnostics (per-thread counters; see math/integrate.hpp).
    m.def("max_depth_hits", &kaspectra::math::max_depth_hits,
          "Panels that exhausted max_depth WITHOUT meeting tolerance since the "
          "last reset -- the signature of a possibly under-converged integral.");
    m.def("reset_max_depth_hits", &kaspectra::math::reset_max_depth_hits,
          "Zero the max_depth_hits counter for this thread.");

    bind_io(m);
    bind_pp(m);
    bind_pgamma(m);
    bind_bh(m);

    // Combined pp + pgamma (+ bh for e+/e-) spectrum; the pure-Python
    // kaspectra.q_total_species wrapper adds species-string coercion.
    // Registered AFTER bind_io so the signature renders the Python
    // ProtonSpectrum/PhotonField types instead of raw C++ names (matters for
    // stub generation).
    m.def("_q_total_species", &kaspectra::q_total_species,
          py::arg("species"), py::arg("E"), py::arg("J_p"), py::arg("n_H"),
          py::arg("f_ph"), py::arg("E_p_max"), py::arg("E_pi_max"), py::arg("epsilon_max"),
          "Sum of all implemented channels for one species; each channel uses "
          "its own tuned default tolerances.");
}
