#include <kaspectra/bh/source.hpp>
#include <kaspectra/io/proton_spectrum.hpp>
#include <kaspectra/io/photon_field.hpp>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

using kaspectra::io::PowerLawSpectrum;
using kaspectra::io::BlackbodyPhotonField;
using kaspectra::bh::interaction_rate;
using kaspectra::bh::energy_loss_rate;

namespace {

    void print_usage(const char* prog) {
        std::cerr <<
            "Usage: " << prog << " [options]\n"
            "  --norm <v>    proton spectrum normalization at E_ref=1 GeV   (default 1.0)\n"
            "  --index <v>   proton spectral index (J_p ~ E_p^-index)      (default 2.0)\n"
            "  --T <v>       photon field blackbody temperature [K]        (default 2.725, CMB)\n"
            "  --epsmax <v>  maximum photon energy epsilon_max [GeV]       (default 1e-6)\n"
            "  --Emin <v>    minimum proton energy E_p [GeV]                (default 1e9)\n"
            "  --Emax <v>    maximum proton energy E_p [GeV]                (default 1e13)\n"
            "  --n <v>       number of log-spaced energy points             (default 50)\n"
            "  -h, --help    print this message and exit\n";
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
    double norm = 1.0;
    double index = 2.0;
    double const E_ref = 1.0;
    double T_kelvin = 2.725;
    double epsilon_max = 1e-6;
    double E_min = 1e9;
    double E_max = 1e13;
    int n_points = 50;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--norm")        norm = next_double(argc, argv, i);
        else if (arg == "--index")  index = next_double(argc, argv, i);
        else if (arg == "--T")      T_kelvin = next_double(argc, argv, i);
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

    if (n_points < 2 || E_min <= 0.0 || E_max <= E_min) {
        std::cerr << "error: require n >= 2 and 0 < Emin < Emax\n";
        return 1;
    }

    PowerLawSpectrum J_p(norm, index, E_ref);
    BlackbodyPhotonField f_ph(T_kelvin);

    std::cout << "E_p_GeV,interaction_rate_s,energy_loss_rate_GeVs\n";
    std::cout << std::scientific << std::setprecision(6);

    const double log_min = std::log10(E_min);
    const double log_max = std::log10(E_max);

    for (int i = 0; i < n_points; ++i) {
        const double frac = static_cast<double>(i) / (n_points - 1);
        const double E_p = std::pow(10.0, log_min + frac * (log_max - log_min));
        std::cout << E_p << ","
                   << interaction_rate(E_p, f_ph, epsilon_max) << ","
                   << energy_loss_rate(E_p, f_ph, epsilon_max) << "\n";
    }
    return 0;
}
