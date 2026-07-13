// Population-level Bethe-Heitler pair SED: q_pair_spectrum over a dense E_e
// grid, amortized through ONE DNdEeTable2D build (the per-point exact
// q_pair_spectrum would re-evaluate dN_dEe at every outer quadrature node of
// every point). Demonstrates the 2D cache's intended use.
#include <kaspectra/bh/spectrum_cache.hpp>
#include <kaspectra/io/photon_field.hpp>
#include <kaspectra/io/proton_spectrum.hpp>
#include <kaspectra/constants.hpp>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

using namespace kaspectra;

namespace {

    void print_usage(const char* prog) {
        std::cerr <<
            "Usage: " << prog << " [options]\n"
            "  --norm <v>    proton power-law normalization                 (default 1.0)\n"
            "  --index <v>   proton power-law index                          (default 2.0)\n"
            "  --T <v>       photon field blackbody temperature [K]         (default 2.725, CMB)\n"
            "  --Epmax <v>   proton budget E_p_max [GeV]                     (default 1e12)\n"
            "  --epsmax <v>  maximum photon energy epsilon_max [GeV]        (default 1e-6)\n"
            "  --Emin <v>    minimum electron energy E_e [GeV]               (default 1e4)\n"
            "  --Emax <v>    maximum electron energy E_e [GeV]               (default 1e8)\n"
            "  --n <v>       number of log-spaced E_e points                 (default 50)\n"
            "  -h, --help    print this message and exit\n";
    }

    double next_double(int argc, char** argv, int& i) {
        if (i + 1 >= argc) {
            std::cerr << "error: missing value for " << argv[i] << "\n";
            std::exit(1);
        }
        return std::stod(argv[++i]);
    }

}   // namespace

int main(int argc, char** argv) {
    double norm = 1.0, index = 2.0;
    double T_kelvin = 2.725;
    double E_p_max = 1e12;
    double epsilon_max = 1e-6;
    double E_min = 1e4, E_max = 1e8;
    int n_points = 50;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--norm")        norm = next_double(argc, argv, i);
        else if (arg == "--index")  index = next_double(argc, argv, i);
        else if (arg == "--T")      T_kelvin = next_double(argc, argv, i);
        else if (arg == "--Epmax")  E_p_max = next_double(argc, argv, i);
        else if (arg == "--epsmax") epsilon_max = next_double(argc, argv, i);
        else if (arg == "--Emin")   E_min = next_double(argc, argv, i);
        else if (arg == "--Emax")   E_max = next_double(argc, argv, i);
        else if (arg == "--n")      n_points = static_cast<int>(next_double(argc, argv, i));
        else if (arg == "-h" || arg == "--help") { print_usage(argv[0]); return 0; }
        else {
            std::cerr << "error: unknown argument '" << arg << "'\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (n_points < 2 || E_min <= 0.0 || E_max <= E_min || E_p_max <= 0.0) {
        std::cerr << "error: require n >= 2, 0 < Emin < Emax, and Epmax > 0\n";
        return 1;
    }

    io::PowerLawSpectrum J_p(norm, index);
    io::BlackbodyPhotonField f_ph(T_kelvin);

    auto table = bh::DNdEeTable2D::build(f_ph, epsilon_max, E_p_max, E_min, E_max);
    if (!table.built()) {
        std::cerr << "error: table build failed (degenerate E_e band?)\n";
        return 1;
    }

    std::cout << "E_e_GeV,q_pair_per_GeV_cm3_s\n";
    std::cout << std::scientific << std::setprecision(6);

    const double log_min = std::log10(E_min);
    const double log_max = std::log10(E_max);
    for (int i = 0; i < n_points; ++i) {
        const double frac = static_cast<double>(i) / (n_points - 1);
        const double E_e = std::pow(10.0, log_min + frac * (log_max - log_min));
        std::cout << E_e << ","
                  << bh::q_pair_spectrum_cached(E_e, J_p, table, E_p_max, epsilon_max) << "\n";
    }
    return 0;
}
