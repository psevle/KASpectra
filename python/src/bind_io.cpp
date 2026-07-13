#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "kaspectra/io/proton_spectrum.hpp"
#include "kaspectra/io/photon_field.hpp"

namespace py = pybind11;
using namespace kaspectra::io;

namespace {

    // Trampolines: allow Python subclasses of the abstract bases (a custom
    // J_p(E_p) or f_ph(epsilon) written directly in Python). PYBIND11_OVERRIDE
    // acquires the GIL before calling into Python, so a Python-defined field
    // is safe to evaluate from the GIL-released hot loops (map_array,
    // DNdEeTable's std::async build) -- PROVIDED the thread that launched the
    // work has itself released the GIL first; the expensive entry points in
    // bind_bh.cpp use py::call_guard<py::gil_scoped_release> for exactly that
    // reason (otherwise the launcher would hold the GIL while worker threads
    // block on acquiring it: deadlock). Python-side overrides are, of course,
    // orders of magnitude slower per call than the C++ implementations.
    class PyProtonSpectrum : public ProtonSpectrum {
        public:
            using ProtonSpectrum::ProtonSpectrum;
            double operator()(double E_p) const override {
                PYBIND11_OVERRIDE_PURE_NAME(double, ProtonSpectrum, "__call__", operator(), E_p);
            }
    };

    class PyPhotonField : public PhotonField {
        public:
            using PhotonField::PhotonField;
            double operator()(double epsilon) const override {
                PYBIND11_OVERRIDE_PURE_NAME(double, PhotonField, "__call__", operator(), epsilon);
            }
    };

}   // namespace

void bind_io(py::module_& m) {
    auto io = m.def_submodule("io", "Proton spectrum and photon field input types");

    // Abstract bases: constructible only via Python subclassing (the
    // trampoline); C++ concrete subclasses below need no Python-side ctor.
    py::class_<ProtonSpectrum, PyProtonSpectrum>(io, "ProtonSpectrum",
        "J_p(E_p): differential proton flux/number density [GeV]. Construct one of "
        "the concrete subclasses below, use TabulatedSpectrum for arbitrary data, or "
        "subclass this in Python and implement __call__(self, E_p) -> float. "
        "No unit/normalization convention is enforced -- caller's responsibility.")
        .def(py::init<>())
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

    py::class_<PhotonField, PyPhotonField>(io, "PhotonField",
        "f_ph(epsilon): differential photon number density [GeV^-1 cm^-3]. "
        "Subclassable from Python: implement __call__(self, epsilon) -> float.")
        .def(py::init<>())
        .def("__call__", &PhotonField::operator(), py::arg("epsilon"));

    py::class_<BlackbodyPhotonField, PhotonField>(io, "BlackbodyPhotonField",
        "Isotropic Planckian photon field at temperature T_kelvin [K] (e.g. 2.725 for the CMB).")
        .def(py::init<double>(), py::arg("T_kelvin"))
        .def("kT", &BlackbodyPhotonField::kT, "kT in GeV.");

    py::class_<PowerLawPhotonField, PhotonField>(io, "PowerLawPhotonField",
        "f_ph(epsilon) = norm * (epsilon/epsilon_ref)^-index * exp(-epsilon/epsilon_cutoff)")
        .def(py::init<double, double, double, double>(),
             py::arg("norm"), py::arg("index"), py::arg("epsilon_cutoff"), py::arg("epsilon_ref") = 1.0);

    py::class_<TabulatedPhotonField, PhotonField>(io, "TabulatedPhotonField",
        "Log-log interpolated photon field over tabulated (epsilon, f_ph) points.")
        .def(py::init<const std::vector<double>&, const std::vector<double>&>(),
             py::arg("epsilon"), py::arg("f_ph"));

    py::class_<CompositePhotonField, PhotonField>(io, "CompositePhotonField",
        "Sum of photon fields (e.g. CMB + IR + starlight). add() stores a "
        "reference; keep_alive ties each added field's lifetime to the composite.")
        .def(py::init<>())
        .def("add", &CompositePhotonField::add, py::arg("field"),
             py::keep_alive<1, 2>())
        .def("__len__", &CompositePhotonField::size);
}
