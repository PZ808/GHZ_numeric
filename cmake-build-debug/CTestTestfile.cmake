# CMake generated Testfile for 
# Source directory: /Users/antares/Projects/physics_codes/GHZ_numeric
# Build directory: /Users/antares/Projects/physics_codes/GHZ_numeric/cmake-build-debug
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[spectral_unit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/cmake-build-debug/spectral_unit_tests")
set_tests_properties([=[spectral_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;281;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[spectral_physical_unit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/cmake-build-debug/spectral_physical_unit_tests")
set_tests_properties([=[spectral_physical_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;282;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[spectral_solve_unit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/cmake-build-debug/spectral_solve_unit_tests")
set_tests_properties([=[spectral_solve_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;283;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[cubic_splines_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/cmake-build-debug/cubic_splines_tests")
set_tests_properties([=[cubic_splines_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;284;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[kinnersley_held_operators_unit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/cmake-build-debug/kinnersley_held_operators_unit_tests")
set_tests_properties([=[kinnersley_held_operators_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;285;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[effective_source_unit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/cmake-build-debug/effective_source_unit_tests")
set_tests_properties([=[effective_source_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;286;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[teukolsky_source_unit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/cmake-build-debug/teukolsky_source_unit_tests")
set_tests_properties([=[teukolsky_source_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;287;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[bound_orbit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/cmake-build-debug/bound_orbit_tests")
set_tests_properties([=[bound_orbit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;288;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[corrector_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/cmake-build-debug/corrector_tests")
set_tests_properties([=[corrector_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;289;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[proxy_solver_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/cmake-build-debug/proxy_solver_tests")
set_tests_properties([=[proxy_solver_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;290;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
subdirs("_deps/catch2-build")
