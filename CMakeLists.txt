cmake_minimum_required(VERSION 3.16)
project(KASpectra LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# kaspectra is header-only for now; becomes a compiled lib once src/ gains .cpp files.
add_library(kaspectra INTERFACE)
target_include_directories(kaspectra INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(kaspectra INTERFACE cxx_std_17)

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
