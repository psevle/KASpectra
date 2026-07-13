cmake_minimum_required(VERSION 3.16)
project(KASpectra VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# kaspectra is header-only for now; becomes a compiled lib once src/ gains .cpp files.
add_library(kaspectra INTERFACE)
add_library(KASpectra::kaspectra ALIAS kaspectra)
target_include_directories(kaspectra INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(kaspectra INTERFACE cxx_std_17)

# install + find_package(KASpectra) support for C++ consumers:
#   cmake --install build --prefix <prefix>
#   find_package(KASpectra REQUIRED)
#   target_link_libraries(app PRIVATE KASpectra::kaspectra)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS kaspectra EXPORT KASpectraTargets)
install(EXPORT KASpectraTargets
        NAMESPACE KASpectra::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/KASpectra)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/KASpectraConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/KASpectraConfig.cmake"
    "include(\"\${CMAKE_CURRENT_LIST_DIR}/KASpectraTargets.cmake\")\n")
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/KASpectraConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/KASpectraConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/KASpectra)

option(KASPECTRA_BUILD_TESTS "Build KASpectra unit tests" ON)

if(KASPECTRA_BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.7.1
    )
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

    enable_testing()
    add_subdirectory(test)
endif()

option(KASPECTRA_BUILD_EXAMPLES "Build KASpectra example CLI programs" ON)

if(KASPECTRA_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# OFF by default, unlike the two options above: pybind11 needs Python dev
# headers and a matching interpreter -- a real new dependency -- and a plain
# C++-only build should keep configuring without any Python environment.
option(KASPECTRA_BUILD_PYTHON "Build KASpectra Python bindings (pybind11)" OFF)

if(KASPECTRA_BUILD_PYTHON)
    add_subdirectory(python)
endif()
