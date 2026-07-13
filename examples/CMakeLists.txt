add_executable(pgamma_spectrum pgamma_spectrum.cpp)
target_link_libraries(pgamma_spectrum PRIVATE kaspectra)

add_executable(pp_spectrum pp_spectrum.cpp)
target_link_libraries(pp_spectrum PRIVATE kaspectra)

add_executable(bh_rates bh_rates.cpp)
target_link_libraries(bh_rates PRIVATE kaspectra)

add_executable(bh_spectrum bh_spectrum.cpp)
target_link_libraries(bh_spectrum PRIVATE kaspectra)

add_executable(bh_pair_sed bh_pair_sed.cpp)
target_link_libraries(bh_pair_sed PRIVATE kaspectra)
