#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "kaspectra/species.hpp"
#include "kaspectra/io/proton_spectrum.hpp"
#include "kaspectra/io/photon_field.hpp"
#include "kaspectra/pgamma/source.hpp"

namespace py = pybind11;
using kaspectra::Species;
using kaspectra::io::ProtonSpectrum;
using kaspectra::io::PhotonField;

namespace {

    py::array_t<double> q_species_array(Species s, py::array_t<double, py::array::c_style | py::array::forcecast> E,
                                         const ProtonSpectrum& J_p, const PhotonField& f_ph,
                                         double E_p_max, double epsilon_max,
                                         double abs_tol, double rel_tol, int max_depth) {
        auto E_buf = E.unchecked<1>();
        py::array_t<double> out(E_buf.shape(0));
        auto out_buf = out.mutable_unchecked<1>();
        {
            py::gil_scoped_release release;
            for (py::ssize_t i = 0; i < E_buf.shape(0); ++i) {
                out_buf(i) = kaspectra::pgamma::q_species(s, E_buf(i), J_p, f_ph, E_p_max, epsilon_max,
                                                           abs_tol, rel_tol, max_depth);
            }
        }
        return out;
    }

} // namespace

void bind_pgamma(py::module_& m) {
    auto pgamma = m.def_submodule("pgamma", "Photomeson (pgamma) secondary spectra (KA2008)");

    pgamma.def("q_species", &kaspectra::pgamma::q_species,
               py::arg("species"), py::arg("E"), py::arg("J_p"), py::arg("f_ph"),
               py::arg("E_p_max"), py::arg("epsilon_max"),
               py::arg("abs_tol") = 1e-10, py::arg("rel_tol") = 1e-8, py::arg("max_depth") = 50,
               "Population-level pgamma secondary differential spectrum q(E). "
               "NOTE: same name as kaspectra.pp.q_species but a different signature "
               "(photon field instead of n_H) -- kept in separate submodules, mirroring "
               "the two C++ namespaces, so there is no collision.");

    pgamma.def("q_species", &q_species_array,
               py::arg("species"), py::arg("E"), py::arg("J_p"), py::arg("f_ph"),
               py::arg("E_p_max"), py::arg("epsilon_max"),
               py::arg("abs_tol") = 1e-10, py::arg("rel_tol") = 1e-8, py::arg("max_depth") = 50,
               "Vectorized form of q_species: E may be a numpy array.");
}
