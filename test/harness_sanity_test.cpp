#include <catch2/catch_test_macros.hpp>

// Confirms the CMake + Catch2 harness builds/links/runs before any physics
// code exists. Delete once real tests (math/, pp/) land.
TEST_CASE("build harness is wired up", "[sanity]") {
    REQUIRE(1 + 1 == 2);
}
