#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "kaspectra/io/proton_spectrum.hpp"
#include "kaspectra/io/photon_field.hpp"

namespace py = pybind11;
using namespace kaspectra::io;

void bind_io(py::module_& m) {
    auto io = m.def_submodule("io", "Proton spectrum and photon field input types");

    // Abstract bases: not constructible from Python (no py::init registered).
    // __call__ binds the pure-virtual member; dynamic dispatch to whichever
    // concrete subclass is actually held works exactly as in C++.
    py::class_<ProtonSpectrum>(io, "ProtonSpectrum",
        "J_p(E_p): differential proton flux/number density [GeV]. Construct one of "
        "the concrete subclasses below (or use TabulatedSpectrum for arbitrary data). "
        "No unit/normalization convention is enforced -- caller's responsibility.")
        .def("__call__", &ProtonSpectrum::operator(), py::arg("E_p"));

    py::class_<PowerLawSpectrum, ProtonSpectrum>(io, "PowerLawSpectrum",
        "J_p(E_p) = norm * (E_p/E_ref)^-index")
        .def(py::init<double, double, double>(),
             py::arg("norm"), py::arg("index"), py::arg("E_ref") = 1.0);

    py::class_<ExpCutoffPowerLawSpectrum, ProtonSpectrum>(io, "ExpCutoffPowerLawSpectrum",
        "J_p(E_p) = norm * (E_p/E_ref)^-index * exp(-E_p/E_cutoff)")
        .def(py::init<double, double, double, double>(),
             py::arg("norm"), py::arg("index"), py::arg("E_cutoff"), py::arg("E_ref") = 1.0);

    py::class_<BrokenPowerLawSpectrum, ProtonSpectrum>(io, "BrokenPowerLawSpectrum",
        "Two power laws joined continuously at E_break; norm is defined at E_ref < E_break.")
        .def(py::init<double, double, double, double, double>(),
             py::arg("norm"), py::arg("index_lo"), py::arg("index_hi"),
             py::arg("E_break"), py::arg("E_ref") = 1.0);

    py::class_<TabulatedSpectrum, ProtonSpectrum>(io, "TabulatedSpectrum",
        "Log-log interpolated proton spectrum over tabulated (E, J) points.")
        .def(py::init<const std::vector<double>&, const std::vector<double>&>(),
             py::arg("E"), py::arg("J"));

    py::class_<PhotonField>(io, "PhotonField",
        "f_ph(epsilon): differential photon number density [GeV^-1 cm^-3].")
        .def("__call__", &PhotonField::operator(), py::arg("epsilon"));

    py::class_<BlackbodyPhotonField, PhotonField>(io, "BlackbodyPhotonField",
        "Isotropic Planckian photon field at temperature T_kelvin [K] (e.g. 2.725 for the CMB).")
        .def(py::init<double>(), py::arg("T_kelvin"));

    py::class_<PowerLawPhotonField, PhotonField>(io, "PowerLawPhotonField",
        "f_ph(epsilon) = norm * (epsilon/epsilon_ref)^-index * exp(-epsilon/epsilon_cutoff)")
        .def(py::init<double, double, double, double>(),
             py::arg("norm"), py::arg("index"), py::arg("epsilon_cutoff"), py::arg("epsilon_ref") = 1.0);

    py::class_<TabulatedPhotonField, PhotonField>(io, "TabulatedPhotonField",
        "Log-log interpolated photon field over tabulated (epsilon, f_ph) points.")
        .def(py::init<const std::vector<double>&, const std::vector<double>&>(),
             py::arg("epsilon"), py::arg("f_ph"));
}
