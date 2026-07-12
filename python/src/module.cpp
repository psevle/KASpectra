#include <pybind11/pybind11.h>

#include "kaspectra/species.hpp"
#include "kaspectra/constants.hpp"

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

    bind_io(m);
    bind_pp(m);
    bind_pgamma(m);
    bind_bh(m);
}
