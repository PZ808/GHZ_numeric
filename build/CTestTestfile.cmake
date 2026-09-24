# CMake generated Testfile for 
# Source directory: /Users/antares/Projects/physics_codes/GHZ_numeric
# Build directory: /Users/antares/Projects/physics_codes/GHZ_numeric/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[spectral_unit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/spectral_unit_tests")
set_tests_properties([=[spectral_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;343;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[spectral_physical_unit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/spectral_physical_unit_tests")
set_tests_properties([=[spectral_physical_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;344;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[spectral_solve_unit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/spectral_solve_unit_tests")
set_tests_properties([=[spectral_solve_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;345;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[cubic_splines_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/cubic_splines_tests")
set_tests_properties([=[cubic_splines_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;346;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[kinnersley_held_operators_unit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/kinnersley_held_operators_unit_tests")
set_tests_properties([=[kinnersley_held_operators_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;347;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[effective_source_unit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/effective_source_unit_tests")
set_tests_properties([=[effective_source_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;348;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[teukolsky_source_unit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/teukolsky_source_unit_tests")
set_tests_properties([=[teukolsky_source_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;349;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[bound_orbit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/bound_orbit_tests")
set_tests_properties([=[bound_orbit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;350;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[corrector_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/corrector_tests")
set_tests_properties([=[corrector_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;351;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[proxy_solver_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/proxy_solver_tests")
set_tests_properties([=[proxy_solver_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;352;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[bondi_schwarzschild_data_collocation_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/bondi_schwarzschild_data_collocation_tests" "/Users/antares/Projects/physics_codes/GHZ_numeric/tests/data/psi0_schwarzschild_r0_10_lmax20_mostly_minus.csv")
set_tests_properties([=[bondi_schwarzschild_data_collocation_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;353;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[bondi_zeta_circular_orbit_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/bondi_zeta_circular_orbit_tests" "/Users/antares/Projects/physics_codes/GHZ_numeric/tests/data/psi0_schwarzschild_r0_10_lmax20_mostly_minus.csv")
set_tests_properties([=[bondi_zeta_circular_orbit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;356;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
add_test([=[bondi_metric_falloff_tests]=] "/Users/antares/Projects/physics_codes/GHZ_numeric/build/bondi_metric_falloff_tests" "/Users/antares/Projects/physics_codes/GHZ_numeric/tests/data/psi0_schwarzschild_r0_10_lmax20_mostly_minus.csv")
set_tests_properties([=[bondi_metric_falloff_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;359;add_test;/Users/antares/Projects/physics_codes/GHZ_numeric/CMakeLists.txt;0;")
subdirs("_deps/catch2-build")
