#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "kaspectra/species.hpp"
#include "kaspectra/io/proton_spectrum.hpp"
#include "kaspectra/pp/source.hpp"

namespace py = pybind11;
using kaspectra::Species;
using kaspectra::io::ProtonSpectrum;

namespace {

    // Vectorized q_species: loops in C++ with the GIL released. Safe because
    // nothing in the loop body touches a Python object -- J_p's virtual call
    // and the adaptive quadrature are pure C++/numeric.
    py::array_t<double> q_species_array(Species s, py::array_t<double, py::array::c_style | py::array::forcecast> E,
                                         const ProtonSpectrum& J_p, double n_H,
                                         double E_p_max, double E_pi_max, double K_pi,
                                         double abs_tol, double rel_tol, int max_depth) {
        auto E_buf = E.unchecked<1>();
        py::array_t<double> out(E_buf.shape(0));
        auto out_buf = out.mutable_unchecked<1>();
        {
            py::gil_scoped_release release;
            for (py::ssize_t i = 0; i < E_buf.shape(0); ++i) {
                out_buf(i) = kaspectra::pp::q_species(s, E_buf(i), J_p, n_H, E_p_max, E_pi_max,
                                                       K_pi, abs_tol, rel_tol, max_depth);
            }
        }
        return out;
    }

} // namespace

void bind_pp(py::module_& m) {
    auto pp = m.def_submodule("pp", "Proton-proton secondary spectra (arXiv:astro-ph/0606058)");

    pp.def("q_species", &kaspectra::pp::q_species,
           py::arg("species"), py::arg("E"), py::arg("J_p"), py::arg("n_H"),
           py::arg("E_p_max"), py::arg("E_pi_max"),
           py::arg("K_pi") = kaspectra::pp::K_pi_default,
           py::arg("abs_tol") = 1e-10, py::arg("rel_tol") = 1e-8, py::arg("max_depth") = 50,
           "Population-level pp secondary differential spectrum q(E). Below "
           "E_delta_approx_max (100 GeV), Gamma is stitched to a delta-function "
           "approximation; other species are q_accurate extrapolated past its "
           "fitted validity range. See kaspectra/pp/source.hpp.");

    pp.def("q_species", &q_species_array,
           py::arg("species"), py::arg("E"), py::arg("J_p"), py::arg("n_H"),
           py::arg("E_p_max"), py::arg("E_pi_max"),
           py::arg("K_pi") = kaspectra::pp::K_pi_default,
           py::arg("abs_tol") = 1e-10, py::arg("rel_tol") = 1e-8, py::arg("max_depth") = 50,
           "Vectorized form of q_species: E may be a numpy array.");
}
