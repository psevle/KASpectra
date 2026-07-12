include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(_kaspectra
    src/module.cpp
    src/bind_io.cpp
    src/bind_pp.cpp
    src/bind_pgamma.cpp
    src/bind_bh.cpp
)
target_link_libraries(_kaspectra PRIVATE kaspectra)

# Land the compiled extension inside the pure-Python package directory so
# `import kaspectra` works straight out of the build tree, matching what
# scikit-build-core's packaging expects too.
set_target_properties(_kaspectra PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/kaspectra"
)

# When driven by pip/scikit-build-core (which defines SKBUILD), the wheel's
# pure-Python half comes from pyproject.toml's wheel.packages; the compiled
# extension must be install()ed into that same package directory or the
# wheel ships without it.
if(SKBUILD)
    install(TARGETS _kaspectra LIBRARY DESTINATION kaspectra)
endif()
