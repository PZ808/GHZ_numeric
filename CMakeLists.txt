cmake_minimum_required(VERSION 3.20)

# Force Homebrew LLVM before project() so CMake does not enter a compiler-cache
# reconfigure loop.
set(CMAKE_CXX_COMPILER /opt/homebrew/opt/llvm/bin/clang++)
set(CMAKE_C_COMPILER   /opt/homebrew/opt/llvm/bin/clang)

project(GHZ_numeric LANGUAGES CXX)

# --------------------------
# Compiler settings
# --------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Force Homebrew LLVM
#set(CMAKE_C_COMPILER   "/opt/homebrew/opt/llvm@21/bin/clang")
#set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm@21/bin/clang++")

include(FetchContent)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.14.0
)

add_compile_definitions(TEUK_PRECISION=2)

# --------------------------
# Homebrew include and link paths
# --------------------------
include_directories(
        /opt/homebrew/include
        /opt/homebrew/include/eigen3
        /opt/homebrew/opt/libomp/include  # OpenMP headers
)
link_directories(
        /opt/homebrew/lib
        /opt/homebrew/opt/libomp/lib      # OpenMP library
)

# --------------------------
# OpenMP (Homebrew LLVM)
# --------------------------
set(OpenMP_CXX_FLAGS "-fopenmp -I/opt/homebrew/opt/libomp/include")
set(OpenMP_C_FLAGS   "-fopenmp -I/opt/homebrew/opt/libomp/include")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_C_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP: ${OpenMP_CXX_FLAGS}")
endif()


# --------------------------
# External dependencies
# --------------------------
find_package(Eigen3 REQUIRED)
find_package(FFTW3 REQUIRED)

include(FetchContent)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.14.0
)

# --------------------------
# Executable
# --------------------------
add_executable(GHZ_numeric
        main.cpp

        include/ghz/core/GhzTypes.hpp
        include/ghz/core/EigenTeukNumTraits.hpp
        include/ghz/core/VectorsGHZ.hpp
        include/ghz/core/MathMacros.hpp

        include/ghz/geom/KerrParams.hpp src/ghz/geom/KerrParams.cpp
        include/ghz/geom/Metric.hpp src/ghz/geom/KerrMetric.cpp
        include/ghz/geom/KerrMetric.hpp
        include/ghz/geom/KerrMetricBL.hpp src/ghz/geom/KerrMetricBL.cpp
        include/ghz/geom/Coords.hpp src/ghz/geom/Coords.cpp
        include/ghz/geom/KerrMetricIngoing.hpp src/ghz/geom/KerrMetricIngoing.cpp
        include/ghz/geom/KerrMetricOutgoing.hpp src/ghz/geom/KerrMetricOutgoing.cpp

        include/ghz/geom/DataDomain.hpp
        include/ghz/geom/Domain.hpp

        include/ghz/geom/Tetrads.hpp src/ghz/geom/KinnersleyTetrad.cpp
        include/ghz/geom/KinnersleyTetrad.hpp

        include/ghz/ghp/SpinCoeffsNP.hpp src/ghz/ghp/SpinCoeffNP.cpp
        include/ghz/ghp/WeylScalars.hpp src/ghz/ghp/WeylScalars.cpp
        include/ghz/ghp/GHPScalars.hpp src/ghz/ghp/GHPScalars.cpp


        include/ghz/special/EllipticIntegrals.hpp
        include/ghz/orbit/KerrOrbit.hpp src/ghz/orbit/KerrOrbit.cpp

        include/ghz/spectral/SpectralDiffer.hpp src/ghz/spectral/SpectralDiffer.cpp
        include/ghz/ghp/HeldScalars.hpp src/ghz/ghp/HeldScalars.cpp
        include/ghz/geom/KerrCharts.hpp

        include/ghz/spectral/SpectralGHPFieldVectorized.hpp
        src/ghz/spectral/KinnersleySpectralHeldOperators.cpp

        include/ghz/orbit/Splines.hpp

        include/ghz/transport/Corrector.hpp
        src/ghz/transport/Corrector.cpp
        #include/ghz/transport/RK/ODE.hpp
        #include/ghz/transport/RK/Builders.hpp
        #include/ghz/transport/RK/ZSliceSolver.hpp

        src/ghz/transport/SourceBuilders.cpp
        include/ghz/source/EffectiveSource.hpp

        include/ghz/source/BinaryEffectiveSourceArchive.hpp
        src/ghz/source/BinaryEffectiveSourceArchive.cpp

        include/ghz/source/ConditionSource.hpp
        src/ghz/source/ConditionSource.cpp
        include/ghz/source/TeukolskySource.hpp
        src/ghz/source/TeukolskySource.cpp
        include/ghz/source/AngularProjection.hpp
        src/ghz/source/AngularProjection.cpp

        include/ghz/asymptotic/SchwarzschildPsi0Modes.hpp
        src/ghz/asymptotic/SchwarzschildPsi0Modes.cpp
        include/ghz/asymptotic/BondiHeldSeedSolve.hpp
        src/ghz/asymptotic/BondiHeldSeedSolve.cpp
        include/ghz/asymptotic/SchwarzschildBondiGauge.hpp
        src/ghz/asymptotic/SchwarzschildBondiGauge.cpp

        include/ghz/transport/collocation/SpectralSliceSolve.hpp
        src/ghz/transport/collocation/SpectralSliceSolve.cpp

        include/ghz/source/ConditionBoundaryData.hpp

        src/ghz/source/BinaryBoundaryDataArchive.cpp
        src/ghz/source/ConditionBoundaryData.cpp
        src/ghz/transport/collocation/XmmbarWorldtubeSolve.cpp
        src/ghz/spectral/PoleFactorizedSpectralHeldOperators.cpp include/ghz/core/Floating.hpp)


# --------------------------
# Target include directories
# --------------------------
target_include_directories(GHZ_numeric PRIVATE include)

# --------------------------
# Link libraries
# --------------------------
target_link_libraries(GHZ_numeric PRIVATE
        OpenMP::OpenMP_CXX
        fftw3 fftw3l m
        )
add_executable(spectral_unit_tests
        tests/unit/spectral_unit_tests.cpp
        src/ghz/spectral/SpectralDiffer.cpp
        src/ghz/ghp/GHPScalars.cpp
        src/ghz/ghp/SpinCoeffNP.cpp
        )
add_executable(spectral_physical_unit_tests
        tests/unit/spectral_physical_unit_tests.cpp
        src/ghz/spectral/SpectralDiffer.cpp
        src/ghz/ghp/GHPScalars.cpp
        src/ghz/ghp/SpinCoeffNP.cpp
        )
add_executable(spectral_solve_unit_tests
        tests/unit/spectral_solve_unit_tests.cpp
        src/ghz/transport/collocation/SpectralSliceSolve.cpp
        src/ghz/spectral/SpectralDiffer.cpp
        src/ghz/ghp/GHPScalars.cpp
        src/ghz/ghp/SpinCoeffNP.cpp
        )
add_executable(cubic_splines_tests
        tests/unit/cubic_splines_tests.cpp
        )
add_executable(kinnersley_held_operators_unit_tests
        tests/unit/kinnersley_held_operators_unit_tests.cpp
        src/ghz/spectral/KinnersleySpectralHeldOperators.cpp
        src/ghz/spectral/SpectralDiffer.cpp
        src/ghz/ghp/GHPScalars.cpp
        src/ghz/ghp/SpinCoeffNP.cpp
        src/ghz/geom/KerrParams.cpp
        src/ghz/geom/KerrMetric.cpp
        src/ghz/geom/KerrMetricBL.cpp
        src/ghz/geom/Coords.cpp
        src/ghz/geom/KinnersleyTetrad.cpp
        src/ghz/ghp/WeylScalars.cpp
        src/ghz/ghp/HeldScalars.cpp
        src/ghz/spectral/PoleFactorizedSpectralHeldOperators.cpp
        )
add_executable(effective_source_unit_tests
        tests/unit/effective_source_unit_tests.cpp
        src/ghz/source/BinaryEffectiveSourceArchive.cpp src/ghz/source/ConditionSource.cpp
        )
add_executable(teukolsky_source_unit_tests
        tests/unit/teukolsky_source_unit_tests.cpp
        src/ghz/source/TeukolskySource.cpp
        src/ghz/source/AngularProjection.cpp
        src/ghz/geom/KerrParams.cpp
        src/ghz/geom/KerrMetric.cpp
        src/ghz/geom/Coords.cpp
        src/ghz/geom/KinnersleyTetrad.cpp
        src/ghz/ghp/SpinCoeffNP.cpp
        src/ghz/ghp/GHPScalars.cpp
        src/ghz/ghp/WeylScalars.cpp
        src/ghz/ghp/HeldScalars.cpp
        )
add_executable(bound_orbit_tests
        tests/unit/bound_orbit_tests.cpp
        src/ghz/orbit/KerrOrbit.cpp
        src/ghz/geom/KerrMetric.cpp
        src/ghz/geom/KerrParams.cpp
        src/ghz/spectral/PoleFactorizedSpectralHeldOperators.cpp
        src/ghz/spectral/KinnersleySpectralHeldOperators.cpp
        src/ghz/spectral/SpectralDiffer.cpp
        src/ghz/ghp/GHPScalars.cpp
        src/ghz/ghp/SpinCoeffNP.cpp
        src/ghz/geom/Coords.cpp
        src/ghz/geom/KinnersleyTetrad.cpp
        src/ghz/ghp/WeylScalars.cpp
        src/ghz/ghp/HeldScalars.cpp)
add_executable(corrector_tests
        tests/unit/corrector_tests.cpp
        src/ghz/transport/collocation/XmmbarWorldtubeSolve.cpp
        src/ghz/transport/collocation/SpectralSliceSolve.cpp
        src/ghz/source/ConditionSource.cpp
        src/ghz/source/BinaryBoundaryDataArchive.cpp
        src/ghz/source/ConditionBoundaryData.cpp
        src/ghz/transport/Corrector.cpp
        src/ghz/spectral/SpectralDiffer.cpp
        src/ghz/ghp/GHPScalars.cpp
        src/ghz/geom/KerrMetric.cpp
        src/ghz/geom/KerrParams.cpp
        src/ghz/ghp/SpinCoeffNP.cpp
        src/ghz/geom/KerrMetricOutgoing.cpp
        )
add_executable(proxy_solver_tests
        tests/unit/proxy_solver_tests.cpp
        src/ghz/transport/Corrector.cpp
        src/ghz/transport/collocation/SpectralSliceSolve.cpp
        src/ghz/spectral/SpectralDiffer.cpp
        src/ghz/ghp/GHPScalars.cpp
        src/ghz/ghp/SpinCoeffNP.cpp
        )
add_executable(bondi_schwarzschild_data_collocation_tests
        tests/unit/bondi_schwarzschild_data_collocation_tests.cpp
        src/ghz/asymptotic/SchwarzschildPsi0Modes.cpp
        src/ghz/asymptotic/BondiHeldSeedSolve.cpp
        src/ghz/asymptotic/BondiHeldMetricGrid.cpp
        src/ghz/asymptotic/SchwarzschildBondiMetric.cpp
        src/ghz/spectral/KinnersleySpectralHeldOperators.cpp
        src/ghz/spectral/PoleFactorizedSpectralHeldOperators.cpp
        src/ghz/spectral/SpectralDiffer.cpp
        src/ghz/ghp/GHPScalars.cpp
        src/ghz/ghp/SpinCoeffNP.cpp
        src/ghz/geom/KerrParams.cpp
        src/ghz/geom/KerrMetric.cpp
        src/ghz/geom/Coords.cpp
        src/ghz/geom/KinnersleyTetrad.cpp
        src/ghz/ghp/WeylScalars.cpp
        src/ghz/ghp/HeldScalars.cpp)
add_executable(bondi_zeta_circular_orbit_tests
        tests/unit/bondi_zeta_circular_orbit_tests.cpp
        src/ghz/asymptotic/SchwarzschildPsi0Modes.cpp
        src/ghz/asymptotic/SchwarzschildBondiGauge.cpp
        src/ghz/orbit/KerrOrbit.cpp
        src/ghz/geom/KerrMetric.cpp
        src/ghz/geom/KerrParams.cpp)
add_executable(bondi_plot_data
        tools/bondi_plot_data.cpp
        src/ghz/asymptotic/SchwarzschildPsi0Modes.cpp
        src/ghz/asymptotic/SchwarzschildBondiGauge.cpp
        src/ghz/asymptotic/BondiHeldSeedSolve.cpp
        src/ghz/orbit/KerrOrbit.cpp
        src/ghz/spectral/KinnersleySpectralHeldOperators.cpp
        src/ghz/spectral/PoleFactorizedSpectralHeldOperators.cpp
        src/ghz/spectral/SpectralDiffer.cpp
        src/ghz/ghp/GHPScalars.cpp
        src/ghz/ghp/SpinCoeffNP.cpp
        src/ghz/geom/KerrParams.cpp
        src/ghz/geom/KerrMetric.cpp
        src/ghz/geom/Coords.cpp
        src/ghz/geom/KinnersleyTetrad.cpp
        src/ghz/ghp/WeylScalars.cpp
        src/ghz/ghp/HeldScalars.cpp)
add_executable(bondi_metric_falloff_tests
        tests/unit/bondi_metric_falloff_tests.cpp
        src/ghz/asymptotic/SchwarzschildPsi0Modes.cpp
        src/ghz/asymptotic/SchwarzschildBondiMetric.cpp)


FetchContent_MakeAvailable(Catch2)


include(CTest)
#include(Catch)


target_include_directories(spectral_physical_unit_tests PRIVATE include)
target_include_directories(spectral_solve_unit_tests PRIVATE include)
target_include_directories(cubic_splines_tests PRIVATE include)
target_include_directories(corrector_tests PRIVATE include)
target_include_directories(proxy_solver_tests PRIVATE include)

target_include_directories(spectral_unit_tests PRIVATE include)
target_include_directories(effective_source_unit_tests PRIVATE include)
target_include_directories(teukolsky_source_unit_tests PRIVATE include)
target_include_directories(kinnersley_held_operators_unit_tests PRIVATE include)
target_include_directories(bound_orbit_tests PRIVATE include)
target_include_directories(bondi_schwarzschild_data_collocation_tests PRIVATE include)
target_include_directories(bondi_zeta_circular_orbit_tests PRIVATE include)
target_include_directories(bondi_plot_data PRIVATE include)
target_include_directories(bondi_metric_falloff_tests PRIVATE include)

target_link_libraries(bound_orbit_tests PRIVATE  OpenMP::OpenMP_CXX
        fftw3l)

target_link_libraries(corrector_tests PRIVATE  OpenMP::OpenMP_CXX
        fftw3l)
target_link_libraries(proxy_solver_tests PRIVATE OpenMP::OpenMP_CXX fftw3l)
target_link_libraries(kinnersley_held_operators_unit_tests PRIVATE
        OpenMP::OpenMP_CXX  Catch2::Catch2WithMain
        fftw3)
target_link_libraries(bondi_schwarzschild_data_collocation_tests PRIVATE
        OpenMP::OpenMP_CXX fftw3)
target_link_libraries(bondi_zeta_circular_orbit_tests PRIVATE
        OpenMP::OpenMP_CXX fftw3l)
target_link_libraries(bondi_plot_data PRIVATE OpenMP::OpenMP_CXX fftw3l)
target_link_libraries(bondi_metric_falloff_tests PRIVATE OpenMP::OpenMP_CXX fftw3l)
#catch_discover_tests(kinnersley_held_operators_unit_tests)

add_test(NAME spectral_unit_tests COMMAND spectral_unit_tests)
add_test(NAME spectral_physical_unit_tests COMMAND spectral_physical_unit_tests)
add_test(NAME spectral_solve_unit_tests COMMAND spectral_solve_unit_tests)
add_test(NAME cubic_splines_tests COMMAND cubic_splines_tests)
add_test(NAME kinnersley_held_operators_unit_tests COMMAND kinnersley_held_operators_unit_tests)
add_test(NAME effective_source_unit_tests COMMAND effective_source_unit_tests)
add_test(NAME teukolsky_source_unit_tests COMMAND teukolsky_source_unit_tests)
add_test(NAME bound_orbit_tests COMMAND bound_orbit_tests)
add_test(NAME corrector_tests COMMAND corrector_tests)
add_test(NAME proxy_solver_tests COMMAND proxy_solver_tests)
add_test(NAME bondi_schwarzschild_data_collocation_tests
        COMMAND bondi_schwarzschild_data_collocation_tests
        ${CMAKE_SOURCE_DIR}/tests/data/psi0_schwarzschild_r0_10_lmax20_mostly_minus.csv)
add_test(NAME bondi_zeta_circular_orbit_tests
        COMMAND bondi_zeta_circular_orbit_tests
        ${CMAKE_SOURCE_DIR}/tests/data/psi0_schwarzschild_r0_10_lmax20_mostly_minus.csv)
add_test(NAME bondi_metric_falloff_tests
        COMMAND bondi_metric_falloff_tests
        ${CMAKE_SOURCE_DIR}/tests/data/psi0_schwarzschild_r0_10_lmax20_mostly_minus.csv)
