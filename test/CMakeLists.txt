add_executable(kaspectra_tests
    harness_sanity_test.cpp
    math/interpid1d_test.cpp
    math/integrate_test.cpp
    io/proton_spectrum_test.cpp
    pp/cross_section_test.cpp
    pp/secondary_spectra_test.cpp
    pp/delta_approx_test.cpp
    pp/source_test.cpp
)
target_link_libraries(kaspectra_tests PRIVATE kaspectra Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(kaspectra_tests)
