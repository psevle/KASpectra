add_executable(pgamma_spectrum pgamma_spectrum.cpp)
target_link_libraries(pgamma_spectrum PRIVATE kaspectra)

add_executable(pp_spectrum pp_spectrum.cpp)
target_link_libraries(pp_spectrum PRIVATE kaspectra)
