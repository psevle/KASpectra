#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <limits>

#include "kaspectra/io/proton_spectrum.hpp"
#include "kaspectra/io/photon_field.hpp"
#include "kaspectra/bh/source.hpp"
#include "kaspectra/bh/spectrum.hpp"
#include "kaspectra/bh/spectrum_cache.hpp"
#include "kaspectra/math/integrate.hpp"
#include "kaspectra/constants.hpp"

namespace py = pybind11;
using kaspectra::io::ProtonSpectrum;
using kaspectra::io::PhotonField;
using namespace kaspectra::bh;

namespace {

    template <typename F>
    py::array_t<double> map_array(py::array_t<double, py::array::c_style | py::array::forcecast> x, F&& f) {
        auto x_buf = x.unchecked<1>();
        py::array_t<double> out(x_buf.shape(0));
        auto out_buf = out.mutable_unchecked<1>();
        py::gil_scoped_release release;
        for (py::ssize_t i = 0; i < x_buf.shape(0); ++i) out_buf(i) = f(x_buf(i));
        return out;
    }

} // namespace

void bind_bh(py::module_& m) {
    auto bh = m.def_submodule("bh", "Bethe-Heitler pair-production rates and spectra");

    // Default tolerance args below mirror the C++ headers' own defaults
    // verbatim (bh/source.hpp, bh/spectrum.hpp, bh/spectrum_cache.hpp);
    // if a header default changes, update here too.
    bh.def("interaction_rate", &interaction_rate,
           py::arg("E_p"), py::arg("f_ph"), py::arg("epsilon_max"),
           py::arg("abs_tol") = 1e-10, py::arg("rel_tol") = 1e-8,
           py::arg("max_depth") = 20, py::arg("panels") = kaspectra::math::kDefaultPanels,
           "Single-proton BH pair-production interaction rate [s^-1].");

    bh.def("energy_loss_rate", &energy_loss_rate,
           py::arg("E_p"), py::arg("f_ph"), py::arg("epsilon_max"),
           py::arg("abs_tol") = 1e-10, py::arg("rel_tol") = 1e-8,
           py::arg("max_depth") = 20, py::arg("panels") = kaspectra::math::kDefaultPanels,
           "Single-proton BH energy-loss rate -dE_p/dt [GeV/s].");

    bh.def("q_pair_rate", &q_pair_rate,
           py::arg("J_p"), py::arg("f_ph"), py::arg("E_p_max"), py::arg("epsilon_max"),
           py::arg("abs_tol") = 1e-10, py::arg("rel_tol") = 1e-8,
           py::arg("max_depth") = 20, py::arg("panels") = kaspectra::math::kDefaultPanels,
           "Population-level BH interaction rate, integrated over J_p.");

    bh.def("q_pair_energy_loss", &q_pair_energy_loss,
           py::arg("J_p"), py::arg("f_ph"), py::arg("E_p_max"), py::arg("epsilon_max"),
           py::arg("abs_tol") = 1e-10, py::arg("rel_tol") = 1e-8,
           py::arg("max_depth") = 20, py::arg("panels") = kaspectra::math::kDefaultPanels,
           "Population-level BH energy-loss rate, integrated over J_p.");

    bh.def("dN_dEe", &dN_dEe,
           py::arg("E_e"), py::arg("gamma_p"), py::arg("f_ph"), py::arg("epsilon_max"),
           py::arg("abs_tol") = 1e-25, py::arg("rel_tol") = 1e-4,
           py::arg("max_depth") = 20, py::arg("panels") = 24,
           "Single-proton lab-frame e+/e- spectrum dN/dE_e [GeV^-1]. Expensive: a "
           "3-level nested adaptive integral, up to several seconds per call at "
           "UHECR gamma_p against a CMB-like field. Prefer DNdEeTable for repeated "
           "evaluations at the same gamma_p or E_e.");

    bh.def("dN_dEe_array",
           [](py::array_t<double, py::array::c_style | py::array::forcecast> E_e, double gamma_p,
              const PhotonField& f_ph, double epsilon_max,
              double abs_tol, double rel_tol, int max_depth, int panels) {
               return map_array(E_e, [&](double e) {
                   return dN_dEe(e, gamma_p, f_ph, epsilon_max, abs_tol, rel_tol, max_depth, panels);
               });
           },
           py::arg("E_e"), py::arg("gamma_p"), py::arg("f_ph"), py::arg("epsilon_max"),
           py::arg("abs_tol") = 1e-25, py::arg("rel_tol") = 1e-4,
           py::arg("max_depth") = 20, py::arg("panels") = 24,
           "Vectorized dN_dEe over an array of E_e -- still one expensive nested "
           "integral PER element; prefer DNdEeTable for a dense sweep.");

    bh.def("E_p_min", &E_p_min, py::arg("E_e"), py::arg("epsilon_max"),
           "Minimum proton energy for which E_e is kinematically reachable "
           "(returns +inf if never reachable at any E_p).");

    bh.def("q_pair_spectrum", &q_pair_spectrum,
           py::arg("E_e"), py::arg("J_p"), py::arg("f_ph"),
           py::arg("E_p_max"), py::arg("epsilon_max"),
           py::arg("abs_tol") = 1e-25, py::arg("rel_tol") = 1e-3,
           py::arg("max_depth") = 10, py::arg("panels") = 6,
           "Population-level BH pair spectrum -- slow but exact (calls dN_dEe at "
           "every outer quadrature node). See q_pair_spectrum_cached for repeated calls.");

    py::enum_<CacheAxis>(bh, "CacheAxis")
        .value("GammaP", CacheAxis::GammaP)
        .value("ElectronEnergy", CacheAxis::ElectronEnergy);

    py::class_<EeWindow>(bh, "EeWindow",
        "Kinematically reachable E_e window at fixed gamma_p (see E_e_window).")
        .def_readonly("lo", &EeWindow::lo)
        .def_readonly("hi", &EeWindow::hi)
        .def_readonly("reachable", &EeWindow::reachable);

    // The C++ E_e_window takes m_e explicitly (no default); the Python-side
    // default supplies the library's own constants::m_e, which is what every
    // C++ call site passes anyway.
    bh.def("E_e_window", &E_e_window,
           py::arg("gamma_p"), py::arg("epsilon_max"), py::arg("m_e") = kaspectra::constants::m_e,
           "Closed-form kinematically reachable E_e window at fixed gamma_p.");

    bh.def("suggested_n_points", &suggested_n_points,
           py::arg("x_lo"), py::arg("x_hi"), py::arg("points_per_decade") = 24,
           "Grid-density heuristic for DNdEeTable build points (default 24/decade).");

    py::class_<DNdEeTable>(bh, "DNdEeTable",
        "Cached/interpolated stand-in for dN_dEe. Build once via over_gamma_p or "
        "over_E_e, then call repeatedly -- cheap log-log interpolation instead of "
        "the underlying nested integral. Immutable after construction.")
        .def(py::init<>(), "Empty table: reachable()==False, __call__ always returns 0.0.")
        .def_static("over_gamma_p", &DNdEeTable::over_gamma_p,
                    py::arg("E_e"), py::arg("f_ph"), py::arg("epsilon_max"), py::arg("gamma_p_max"),
                    py::arg("n_points") = 0,
                    py::arg("abs_tol") = 1e-25, py::arg("rel_tol") = 1e-4,
                    py::arg("max_depth") = 20, py::arg("panels") = 24,
                    "Build over gamma_p at fixed E_e (feeds q_pair_spectrum_cached's outer integral).")
        .def_static("over_E_e", &DNdEeTable::over_E_e,
                    py::arg("gamma_p"), py::arg("f_ph"), py::arg("epsilon_max"),
                    py::arg("E_e_min") = 0.0,
                    py::arg("E_e_max") = std::numeric_limits<double>::infinity(),
                    py::arg("n_points") = 0,
                    py::arg("abs_tol") = 1e-25, py::arg("rel_tol") = 1e-4,
                    py::arg("max_depth") = 20, py::arg("panels") = 24,
                    "Build over E_e at fixed gamma_p (feeds a dense E_e sweep/plot).")
        .def("__call__", &DNdEeTable::operator(), py::arg("x"))
        .def("__call__", [](const DNdEeTable& t, py::array_t<double, py::array::c_style | py::array::forcecast> x) {
                 return map_array(x, [&](double xi) { return t(xi); });
             }, py::arg("x"))
        .def("reachable", &DNdEeTable::reachable)
        .def("axis", &DNdEeTable::axis)
        .def("fixed_value", &DNdEeTable::fixed_value)
        .def("domain_lo", &DNdEeTable::domain_lo)
        .def("domain_hi", &DNdEeTable::domain_hi)
        .def("__len__", &DNdEeTable::size);

    bh.def("make_table_for_q_pair_spectrum", &make_table_for_q_pair_spectrum,
           py::arg("E_e"), py::arg("f_ph"), py::arg("epsilon_max"), py::arg("E_p_max"),
           py::arg("n_points") = 0,
           "Convenience: builds a DNdEeTable via over_gamma_p with gamma_p_max "
           "derived from E_p_max, matching what q_pair_spectrum_cached expects.");

    bh.def("q_pair_spectrum_cached", &q_pair_spectrum_cached,
           py::arg("E_e"), py::arg("J_p"), py::arg("table"),
           py::arg("E_p_max"), py::arg("epsilon_max"),
           py::arg("abs_tol") = 1e-25, py::arg("rel_tol") = 1e-3,
           py::arg("max_depth") = 10, py::arg("panels") = 6,
           "Like q_pair_spectrum, but consults a pre-built DNdEeTable instead of "
           "calling dN_dEe at every outer quadrature node. Raises ValueError "
           "(pybind11's default translation of std::invalid_argument) if the table "
           "doesn't match this E_e/axis, or doesn't cover E_p_max.");
}
