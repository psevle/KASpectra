#include <kaspectra/bh/spectrum.hpp>
#include <kaspectra/bh/spectrum_cache.hpp>
#include <kaspectra/io/photon_field.hpp>
#include <kaspectra/constants.hpp>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

using kaspectra::io::BlackbodyPhotonField;
using kaspectra::bh::dN_dEe;
using kaspectra::bh::DNdEeTable;

namespace {

    void print_usage(const char* prog) {
        std::cerr <<
            "Usage: " << prog << " [options]\n"
            "  --Ep <v>      single proton energy E_p [GeV]                 (default 1e9)\n"
            "  --T <v>       photon field blackbody temperature [K]        (default 2.725, CMB)\n"
            "  --epsmax <v>  maximum photon energy epsilon_max [GeV]       (default 1e-6)\n"
            "  --Emin <v>    minimum electron energy E_e [GeV]              (default 1e2)\n"
            "  --Emax <v>    maximum electron energy E_e [GeV]              (default 1e8)\n"
            "  --n <v>       number of log-spaced energy points             (default 20)\n"
            "  --direct      call dN_dEe per point instead of a cached table\n"
            "  -h, --help    print this message and exit\n"
            "\n"
            "Dense sweeps (--n above the table's own build size) go through a\n"
            "DNdEeTable::over_E_e cache by default: the table costs a fixed number\n"
            "of dN_dEe evaluations (built in parallel), after which every output\n"
            "row is a cheap interpolation. Blackbody fields additionally take the\n"
            "Planckian fast path inside dN_dEe itself.\n";
    }

    double next_double(int argc, char** argv, int& i) {
        if (i + 1 >= argc) {
            std::cerr << "error: missing value for " << argv[i] << "\n";
            std::exit(1);
        }
        return std::stod(argv[++i]);
    }
} // namespace

int main(int argc, char** argv) {
    double E_p = 1e9;
    double T_kelvin = 2.725;
    double epsilon_max = 1e-6;
    double E_min = 1e2;
    double E_max = 1e8;
    int n_points = 20;
    bool direct = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--Ep")          E_p = next_double(argc, argv, i);
        else if (arg == "--T")      T_kelvin = next_double(argc, argv, i);
        else if (arg == "--epsmax") epsilon_max = next_double(argc, argv, i);
        else if (arg == "--Emin")   E_min = next_double(argc, argv, i);
        else if (arg == "--Emax")   E_max = next_double(argc, argv, i);
        else if (arg == "--n")      n_points = static_cast<int>(next_double(argc, argv, i));
        else if (arg == "--direct") direct = true;
        else if (arg == "-h" || arg == "--help") { print_usage(argv[0]); return 0; }
        else {
            std::cerr << "error: unknown argument '" << arg << "'\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (n_points < 2 || E_min <= 0.0 || E_max <= E_min || E_p <= 0.0) {
        std::cerr << "error: require n >= 2, 0 < Emin < Emax, and Ep > 0\n";
        return 1;
    }

    BlackbodyPhotonField f_ph(T_kelvin);
    const double gamma_p = E_p / kaspectra::constants::m_p;

    // Cache pays off once the sweep is denser than the table's own build
    // grid; below that, per-point direct evaluation is both cheaper and free
    // of interpolation error.
    DNdEeTable table;
    const bool use_table = !direct &&
        n_points > kaspectra::bh::suggested_n_points(E_min, E_max);
    if (use_table) {
        table = DNdEeTable::over_E_e(gamma_p, f_ph, epsilon_max, E_min, E_max);
    }

    std::cout << "E_e_GeV,dN_dEe_per_GeV\n";
    std::cout << std::scientific << std::setprecision(6);

    const double log_min = std::log10(E_min);
    const double log_max = std::log10(E_max);

    for (int i = 0; i < n_points; ++i) {
        const double frac = static_cast<double>(i) / (n_points - 1);
        const double E_e = std::pow(10.0, log_min + frac * (log_max - log_min));
        const double v = use_table ? table(E_e) : dN_dEe(E_e, gamma_p, f_ph, epsilon_max);
        std::cout << E_e << "," << v << "\n";
    }
    return 0;
}
