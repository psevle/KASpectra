add_executable(kaspectra_tests
    harness_sanity_test.cpp
    species_test.cpp
    math/interpid1d_test.cpp
    math/integrate_test.cpp
    io/proton_spectrum_test.cpp
    io/photon_field_test.cpp
    pp/cross_section_test.cpp
    pp/secondary_spectra_test.cpp
    pp/delta_approx_test.cpp
    pp/source_test.cpp
    pgamma/kinematics_test.cpp
    pgamma/phi_shape_test.cpp
    pgamma/phi_tables_test.cpp
    pgamma/secondary_spectra_test.cpp
    pgamma/source_test.cpp
    pgamma/multiplicity_test.cpp
    bh/fits_test.cpp
    bh/source_test.cpp
    bh/spectrum_test.cpp
    bh/spectrum_cache_test.cpp
)
target_link_libraries(kaspectra_tests PRIVATE kaspectra Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(kaspectra_tests)
