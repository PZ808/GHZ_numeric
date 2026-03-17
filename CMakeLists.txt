cmake_minimum_required(VERSION 3.20)
project(GHZ_numeric LANGUAGES CXX)

# --------------------------
# Compiler settings
# --------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Force Homebrew LLVM
set(CMAKE_C_COMPILER   "/opt/homebrew/opt/llvm@20/bin/clang")
set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm@20/bin/clang++")

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
        include/ghz/orbit/KerrOrbit.hpp src/ghz/source/KerrOrbit.cpp

        include/ghz/spectral/SpectralDiffer.hpp src/ghz/spectral/SpectralDiffer.cpp
        include/ghz/ghp/HeldScalars.hpp src/ghz/ghp/HeldScalars.cpp
        include/ghz/geom/KerrCharts.hpp

        include/ghz/spectral/SpectralGHPFieldVectorized.hpp
        src/ghz/spectral/KinnersleySpectralHeldOperators.cpp

        include/ghz/orbit/Splines.hpp

        include/ghz/transport/ODE.hpp
        include/ghz/transport/Corrector.hpp
        src/ghz/transport/Corrector.cpp
        include/ghz/transport/Builders.hpp
        include/ghz/transport/ZSliceSolver.hpp
        #include/ghz/transport/Operators.hpp
        # src/ghz/transport/Operators.cpp

        # src/ghz/transport/Operators.cpp
        src/ghz/transport/SourceBuilders.cpp
        include/ghz/source/EffectiveSource.hpp
        include/ghz/source/BinaryEffectiveSourceArchive.hpp
        src/ghz/source/BinaryEffectiveSourceArchive.cpp
        )

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
add_executable(kinnersley_held_operators_unit_tests
        tests/unit/kinnersley_held_operators_unit_tests.cpp
        src/ghz/spectral/KinnersleySpectralHeldOperators.cpp
        src/ghz/spectral/SpectralDiffer.cpp
        src/ghz/ghp/GHPScalars.cpp
        src/ghz/ghp/SpinCoeffNP.cpp
        src/ghz/geom/KerrParams.cpp
        src/ghz/geom/KerrMetric.cpp
        src/ghz/geom/Coords.cpp
        src/ghz/geom/KinnersleyTetrad.cpp
        src/ghz/ghp/WeylScalars.cpp
        src/ghz/ghp/HeldScalars.cpp
        )
add_executable(effective_source_unit_tests
        tests/unit/effective_source_unit_tests.cpp
        src/ghz/source/BinaryEffectiveSourceArchive.cpp)
target_include_directories(spectral_unit_tests PRIVATE include)
target_include_directories(effective_source_unit_tests PRIVATE include)
target_include_directories(kinnersley_held_operators_unit_tests PRIVATE include)

target_link_libraries(kinnersley_held_operators_unit_tests PRIVATE
        OpenMP::OpenMP_CXX
        fftw3)